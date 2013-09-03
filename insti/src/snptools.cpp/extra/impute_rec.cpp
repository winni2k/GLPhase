/**
@file	impute.cpp
@brief	Haplotype imputation by cFDSL distribution
@author	Yi Wang
@date	04/01/2011
*/
#include	<gsl/gsl_rng.h>
#include	<sys/time.h>
#include	<algorithm>
#include	<iostream>
#include	<stdint.h>
#include	<unistd.h>
#include	<cstdlib>
#include	<cstring>
#include	<sstream>
#include	<tabix.h>
#include	<cfloat>
#include	<cstdio>
#include	<vector>
#include	<zlib.h>
#include	<cmath>
#include	<ctime>
#include	<omp.h>
#define	WordShift	6
#define	WordMod	63
#define	MaskHA	0x80
#define	MaskHB	0x40
#define	MaskG	0x3f
using	namespace	std;
typedef	float	fast;
typedef	double	real;
typedef	unsigned	uint;
typedef uint64_t	word;
const	fast	norm=powf(FLT_MIN,	2.0f/3.0f);

struct	Site{
	string	chr,	all;
	uint	pos;
};

struct	Rate{
	uint	pos;
	real	dis;
	bool	operator()(Rate	X,	Rate	Y){	return	X.pos<Y.pos;	}
};

class	Impute{
private:
	gsl_rng	*rng;
	uint	hn,	pn,	en,	wn;
	vector<word>	haps,	hnew;
	vector<uint>	hsum;
	vector<fast>	tran,	emit;
	vector<Rate>	rate;
	
	fast	pc[4][4];	
	void	set1(word	*P,	uint	I){	P[I>>WordShift]|=static_cast<word>(1)<<(I&WordMod);	}
	void	set0(word	*P,	uint	I){	P[I>>WordShift]&=~(static_cast<word>(1)<<(I&WordMod));	}
	word	test(word	*P,	uint	I){	return	(P[I>>WordShift]>>(I&WordMod))&static_cast<word>(1);	}
	void	pos2dis(void);
	fast	hmm_like(uint	I,	uint	*P);
	fast	hmm(uint	I,	uint	*P,	fast	S);
	void	hmm_work(uint	I,	uint	*P,	fast	S);
	fast	solve(uint	I,	uint	N,	fast	S);
	void	replace(uint	I);
	void	result(void);
public:
	uint	in,	mn,	bn,	sn,	nn;
	real	density;
	vector<string>	vcf_file;
	vector<Site>	site;
	vector<string>	name;
	vector<real>	posi;
	vector<fast>	prob;
	
	Impute(){	rng=gsl_rng_alloc(gsl_rng_default);	gsl_rng_set(rng,	time(NULL));	}
	~Impute(){	gsl_rng_free(rng);	}
	bool	load_rec(const	char	*F);
	bool	load_bin(const	char	*F);
	uint	load_vcf(const	char	*F);
	void	initialize(void);
	void	estimate(void);
	void	save_vcf(const	char	*F);
	void	document(void);
};

bool	Impute::load_rec(const	char	*F){
	Rate	r;	
	FILE	*f=fopen(F,	"rt");
	if(f==NULL||fscanf(f,	"%*s%*s%*s%*s"))	return	false;
	for(fscanf(f,	"%*s%u%*s%lf",	&r.pos,	&r.dis);	!feof(f);	fscanf(f,	"%*s%u%*s%lf",	&r.pos,	&r.dis))	rate.push_back(r);
	fclose(f);
	sort(rate.begin(),	rate.end(),	Rate());
	return	true;
}
	
bool	Impute::load_bin(const	char	*F){
	name.clear();	site.clear();	prob.clear();	posi.clear();
	gzFile	f=gzopen(F,	"rt");
	if(f==Z_NULL)	return	false;
	char	*buffer=new	char[1024*1024];
	{
		gzgets(f,	buffer,	1024*1024);
		string	s;	istringstream	si(buffer);
		for(si>>s>>s>>s>>s;	!si.eof();	si>>s)	name.push_back(s);
		in=name.size();
	}
	Site	temp;	fast	rawp;
	while(gzgets(f,	buffer,	1024*1024)!=NULL){
		istringstream	si(buffer);
		si>>temp.chr>>temp.pos>>temp.all;
		site.push_back(temp);	posi.push_back(temp.pos);
		for(uint	i=0;	i<2*in;	i++){	si>>rawp;	prob.push_back(rawp);	}
	}
	mn=site.size();
	delete	[]	buffer;
	gzclose(f);
	cerr<<F<<endl;
	cerr<<"sites\t"<<mn<<endl;
	cerr<<"sample\t"<<in<<endl;
	return	true;
}

uint	Impute::load_vcf(const	char	*F){
	uint	known=0;
	if(!vcf_file.size())	return	0;
	tabix_t *t;	ti_iter_t	iter;	const	char	*s;	int	len;
	if((t=ti_open(F, 0))==0)	return	0;
	vector<uint16_t>	vid;
	uint	vin;
	{
		vector<string>	vname;
		iter=ti_query(t, 0, 0, 0);
		while((s=ti_read(t,	iter,	&len))!=0){
			if(*s!='#') break;
			else	if(strstr(s,	"#CHROM")==NULL)	continue;
			istringstream	si(s);
			string	st;
			si>>st>>st>>st>>st>>st>>st>>st>>st>>st;
			while(!si.eof()){	si>>st;	vname.push_back(st);	}
		}
		ti_iter_destroy(iter);
		vin=vname.size();	vid.resize(vin);
		for(uint	i=0;	i<vin;	i++){
			vid[i]=0xffff;
			for(uint	j=0;	j<in;	j++)	if(vname[i]==name[j]){	vid[i]=j;	break;	}
		}
	}
	iter=ti_query(t,	site[0].chr.c_str(),	site[0].pos,	site[mn-1].pos);
	while((s=ti_read(t,	iter,	&len))!=0){
		istringstream	si(s);	string	st,	chr,	all;	uint	pos;
		si>>chr>>pos>>st>>all>>st;	all+=st;
		si>>st>>st>>st>>st;
		for(uint	m=0;	m<mn;	m++)	if(site[m].chr==chr&&site[m].pos==pos&&site[m].all==all){
			for(uint	i=0;	i<vin;	i++){
				si>>st;	if(vid[i]==0xffff)	continue;
				fast	*p=&prob[m*in*2+vid[i]*2],	pa=fmaxf(1-p[0]-p[1],0);
				if(st[0]=='0'&&st[2]=='0'){	pa*=0.998;	p[0]*=1e-3;	p[1]*=1e-3;	}
				else	if((st[0]=='0'&&st[2]=='1')||(st[0]=='1'&&st[2]=='0')){	pa*=1e-3;	p[0]*=0.998;	p[1]*=1e-3;	}
				else	if(st[0]=='1'&&st[2]=='1'){	pa*=1e-3;	p[0]*=1e-3;	p[1]*=0.998;	}
				fast	sum=1/(pa+p[0]+p[1]);
				p[0]*=sum;	p[1]*=sum;	known++;
			}
			break;
		}
	}	
	ti_iter_destroy(iter);
	ti_close(t);
	return	known;
}

void	Impute::pos2dis(void){
	Rate	r;
	for(uint	m=0;	m<mn;	m++){
		r.pos=(uint)(posi[m]+0.5);
		vector<Rate>::iterator	vi=lower_bound(rate.begin(),	rate.end(),	r,	Rate());
		if(vi==rate.begin()){	posi[m]=vi->dis+1e-6*(posi[m]-vi->pos);	}
		else	if(vi==rate.end()){	vi--;	posi[m]=vi->dis+1e-6*(posi[m]-vi->pos);	}
		else{
			vector<Rate>::iterator	ui=vi;	ui--;
			posi[m]=ui->dis+(posi[m]-ui->pos)/(vi->pos-ui->pos)*(vi->dis-ui->dis);
		}
	}	
}

void	Impute::initialize(void){
	if(rate.size())	pos2dis();
	wn=(mn&WordMod)?(mn>>WordShift)+1:(mn>>WordShift);
	hn=in*2;	haps.resize(hn*wn);	hnew.resize(hn*wn);
	hsum.assign(hn*mn,	0);
	pn=3*mn;	tran.resize(pn);	vector<fast>	temp(in*pn);	
	en=4*mn;	emit.resize(in*en);	
	
	fast	mu=0,	rho;
	for(uint	i=1;	i<hn;	i++)	mu+=1.0/i;
	if(posi.size()!=mn){	posi.resize(mn);	for(uint	m=0;	m<mn;	m++)	posi[m]=m;	}
	mu=1/mu;	rho=0.5*mu*(mn-1)/(posi[mn-1]-posi[0])/density;	mu=mu/(hn+mu);
	for(uint	m=mn-1;	m;	m--){
		posi[m]=(posi[m]-posi[m-1])*rho;	fast	r=posi[m]/(posi[m]+hn);
		tran[m*3]=(1-r)*(1-r);	tran[m*3+1]=r*(1-r);	tran[m*3+2]=r*r;
	}
	pc[0][0]=pc[1][1]=pc[2][2]=pc[3][3]=(1-mu)*(1-mu);	
	pc[0][1]=pc[0][2]=pc[1][0]=pc[1][3]=pc[2][0]=pc[2][3]=pc[3][1]=pc[3][2]=mu*(1-mu);
	pc[0][3]=pc[1][2]=pc[2][1]=pc[3][0]=mu*mu;
	
	for(uint	i=0;	i<in;	i++){
		word	*ha=&haps[i*2*wn],	*hb=ha+wn;
		fast	*t=&temp[i*pn],	*e=&emit[i*en],	*p=&prob[i*2];
		for(uint	m=0;	m<mn;	m++,	t+=3,	e+=4,	p+=hn){
			if(gsl_rng_get(rng)&1)	set1(ha,	m);
			if(gsl_rng_get(rng)&1)	set1(hb,	m);
			t[0]=fmaxf(1-p[0]-p[1],0.0f);	t[1]=fmaxf(p[0],0.0f);	t[2]=fmaxf(p[1],0.0f);
			for(uint	j=0;	j<4;	j++)
				e[j]=pc[j][0]*t[0]+pc[j][1]*t[1]+pc[j][2]*t[1]+pc[j][3]*t[2];
			
		}
	}
	swap(temp,	prob);
}

fast	Impute::hmm_like(uint	I,	uint	*P){
	word	*f0=&haps[P[0]*wn],	*f1=&haps[P[1]*wn],	*m0=&haps[P[2]*wn],	*m1=&haps[P[3]*wn];	
	fast	*e=&emit[I*en],	*t=&tran[0],	sum,	score=0;
	fast	l00=0.25f*e[(test(f0,0)<<1)|test(m0,0)],	l01=0.25f*e[(test(f0,0)<<1)|test(m1,0)];
	fast	l10=0.25f*e[(test(f1,0)<<1)|test(m0,0)],	l11=0.25f*e[(test(f1,0)<<1)|test(m1,0)];
	fast	b00,	b01,	b10,	b11;
	e+=4;	t+=3;
	for(uint	m=1;	m<mn;	m++,	e+=4,	t+=3){
		b00=l00*t[0]+(l01+l10)*t[1]+l11*t[2];	b01=l01*t[0]+(l00+l11)*t[1]+l10*t[2];
		b10=l10*t[0]+(l00+l11)*t[1]+l01*t[2];	b11=l11*t[0]+(l01+l10)*t[1]+l00*t[2];
		l00=b00*e[(test(f0,m)<<1)|test(m0,m)];	l01=b01*e[(test(f0,m)<<1)|test(m1,m)];
		l10=b10*e[(test(f1,m)<<1)|test(m0,m)];	l11=b11*e[(test(f1,m)<<1)|test(m1,m)];
		if((sum=l00+l01+l10+l11)<norm){
			sum=1.0f/sum;	score-=logf(sum);
			l00*=sum;	l01*=sum;	l10*=sum;	l11*=sum;
		}
	}
	return	score+logf(l00+l01+l10+l11);
}

void	Impute::hmm_work(uint	I,	uint	*P,	fast	S){
	word	*f0=&haps[P[0]*wn],	*f1=&haps[P[1]*wn],	*m0=&haps[P[2]*wn],	*m1=&haps[P[3]*wn];	
	//	backward recursion
	vector<fast>	beta(mn*4);
	fast	*e=&emit[(I+1)*en-4],	*t=&tran[(mn-1)*3],	sum,	*b=&beta[(mn-1)*4];
	fast	l00=0,	l01=0,	l10=0,	l11=0;
	fast	b00,	b01,	b10,	b11;
	
	b[0]=b[1]=b[2]=b[3]=1;
	for(uint	m=mn-1;	m;	m--,	e-=4,	t-=3){
		b00=b[0]*e[(test(f0,m)<<1)|test(m0,m)];	b01=b[1]*e[(test(f0,m)<<1)|test(m1,m)];
		b10=b[2]*e[(test(f1,m)<<1)|test(m0,m)];	b11=b[3]*e[(test(f1,m)<<1)|test(m1,m)];
		b-=4;
		b[0]=b00*t[0]+(b01+b10)*t[1]+b11*t[2];	b[1]=b01*t[0]+(b00+b11)*t[1]+b10*t[2];
		b[2]=b10*t[0]+(b00+b11)*t[1]+b01*t[2];	b[3]=b11*t[0]+(b01+b10)*t[1]+b00*t[2];
		sum=1.0f/(b[0]+b[1]+b[2]+b[3]);	b[0]*=sum;	b[1]*=sum;	b[2]*=sum;	b[3]*=sum;
	}
	
	//	forward sampling	
	word	*ha=&hnew[I*2*wn],	*hb=ha+wn;
	fast	*p=&prob[I*pn];	e=&emit[I*en];	t=&tran[0];	b=&beta[0];
	uint	s00,	s01,	s10,	s11;
	for(uint	m=0;	m<mn;	m++,	e+=4,	t+=3,	p+=3,	b+=4){
		s00=(test(f0,m)<<1)|test(m0,m);	s01=(test(f0,m)<<1)|test(m1,m);
		s10=(test(f1,m)<<1)|test(m0,m);	s11=(test(f1,m)<<1)|test(m1,m);
		if(m){
			b00=l00*t[0]+l01*t[1]+l10*t[1]+l11*t[2],	b01=l00*t[1]+l01*t[0]+l10*t[2]+l11*t[1];
			b10=l00*t[1]+l01*t[2]+l10*t[0]+l11*t[1],	b11=l00*t[2]+l01*t[1]+l10*t[1]+l11*t[0];
			l00=b00*e[s00];	l01=b01*e[s01];	l10=b10*e[s10];	l11=b11*e[s11];
		}
		else{	l00=0.25f*e[s00];	l01=0.25f*e[s01];	l10=0.25f*e[s10];	l11=0.25f*e[s11];	}	
		sum=1.0f/(l00+l01+l10+l11);	l00*=sum;	l01*=sum;	l10*=sum;	l11*=sum;
		fast	p00=l00*b[0],	p01=l01*b[1],	p10=l10*b[2],	p11=l11*b[3];
		fast	c00=powf(p[0]*(p00*pc[s00][0]+p01*pc[s01][0]+p10*pc[s10][0]+p11*pc[s11][0]),	S);
		fast	c01=powf(p[1]*(p00*pc[s00][1]+p01*pc[s01][1]+p10*pc[s10][1]+p11*pc[s11][1]),	S);
		fast	c10=powf(p[1]*(p00*pc[s00][2]+p01*pc[s01][2]+p10*pc[s10][2]+p11*pc[s11][2]),	S);
		fast	c11=powf(p[2]*(p00*pc[s00][3]+p01*pc[s01][3]+p10*pc[s10][3]+p11*pc[s11][3]),	S);
		sum=gsl_rng_uniform(rng)*(c00+c01+c10+c11);
		if(sum<c00){	set0(ha,m);	set0(hb,m);	}
		else	if(sum<c00+c01){	set0(ha,m);	set1(hb,m);	}
		else	if(sum<c00+c01+c10){	set1(ha,m);	set0(hb,m);	}
		else{	set1(ha,m);	set1(hb,m);	}
	}
}

fast	Impute::solve(uint	I,	uint	N,	fast	S){
	uint	p[4];
	for(uint	j=0;	j<4;	j++){	
		do	p[j]=gsl_rng_get(rng)%hn;	while(p[j]/2==I);
	}
	fast	curr=hmm_like(I,	p);
		
	for(uint	n=0;	n<N;	n++){
		uint	rp=gsl_rng_get(rng)&3,	oh=p[rp];
		do	p[rp]=gsl_rng_get(rng)%hn;	while(p[rp]/2==I);
		fast	prop=hmm_like(I,	p);
		if(prop>curr||gsl_rng_uniform(rng)<expf((prop-curr)*S))	curr=prop;
		else	p[rp]=oh;
	}
	hmm_work(I,	p,	S);
	return	curr;	
}

void	Impute::estimate(void){
	cerr.setf(ios::fixed);	cerr.precision(3);
	for(uint	n=0;	n<bn+sn;	n++){
		fast	sum=0,	pen=fminf(2*(n+1.0f)/bn,	1);	pen*=pen;
		#pragma omp parallel for
		for(uint	i=0;	i<in;	i++)	sum+=solve(i,	nn*in,	pen);
		swap(hnew,	haps);
		if(n>=bn){
			#pragma omp parallel for
			for(uint	i=0;	i<in;	i++)	replace(i);
		}
		cerr<<n<<'\t'<<pen<<'\t'<<sum/in/mn<<'\r';
	}
	cerr<<endl;
	result();	
}

void	Impute::replace(uint	I){
	word	*oa=&haps[I*2*wn],	*ob=oa+wn;
	uint	*ha=&hsum[I*2*mn],	*hb=ha+mn;
	uint	sis=0,	tra=0;
	for(uint	m=0;	m<mn;	m++){	
		if(test(oa,m)){	sis+=ha[m];	tra+=hb[m];	}
		if(test(ob,m)){	sis+=hb[m];	tra+=ha[m];	}
	}
	if(sis>tra)	for(uint	m=0;	m<mn;	m++){	ha[m]+=test(oa,m);	hb[m]+=test(ob,m);	}
	else	for(uint	m=0;	m<mn;	m++){	ha[m]+=test(ob,m);	hb[m]+=test(oa,m);	}	
}

void	Impute::result(void){
	prob.resize(mn*hn);
	fast	norm=1.0/sn;
	for(uint	m=0;	m<mn;	m++){
		fast	*p=&prob[m*hn];
		for(uint	i=0;	i<in;	i++,	p+=2){
			p[0]=hsum[i*2*mn+m]*norm;	p[1]=hsum[(i*2+1)*mn+m]*norm;
		}
	}
}

void	Impute::save_vcf(const	char	*F){
	string	temp=F;	temp+=".vcf.gz";
	gzFile	f=gzopen(temp.c_str(),	"wt");
	gzprintf(f,	"##fileformat=VCFv4.0\n");
	gzprintf(f,	"##source=BCM:SNPTools:impute\n");
	gzprintf(f,	"##reference=1000Genomes-NCBI37\n");
	gzprintf(f,	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	gzprintf(f,	"##FORMAT=<ID=AP,Number=2,Type=Numeric,Description=\"Allelic Probability, P(Allele=1|Haplotype)\">\n");
	gzprintf(f,	"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
	for(uint	i=0;	i<in;	i++)	gzprintf(f,	"\t%s",	name[i].c_str());
	for(uint	m=0;	m<mn;	m++){
		gzprintf(f,	"\n%s\t%u\t.\t%c\t%c\t100\tPASS\t.\tGT:AP",	
			site[m].chr.c_str(),	site[m].pos,	site[m].all[0],	site[m].all[1]);
		fast	*p=&prob[m*hn];	
		for(uint	i=0;	i<in;	i++,	p+=2){
			fast	prr=(1-p[0])*(1-p[1]),	pra=(1-p[0])*p[1]+p[0]*(1-p[1]),	paa=p[0]*p[1];
			if(prr>=pra&&prr>=paa)	gzprintf(f,	"\t0|0:%.3f,%.3f",	p[0],	p[1]);
			else	if(pra>=prr&&pra>=paa){
				if(p[0]>p[1])	gzprintf(f,	"\t1|0:%.3f,%.3f",	p[0],	p[1]);
				else	gzprintf(f,	"\t0|1:%.3f,%.3f",	p[0],	p[1]);
			}
			else	gzprintf(f,	"\t1|1:%.3f,%.3f",	p[0],	p[1]);
		}
	}
	gzprintf(f,	"\n");
	gzclose(f);
}

void	Impute::document(void){
	cerr<<"\nimpute";
	cerr<<"\nhaplotype imputation by cFDSL distribution";
	cerr<<"\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage	impute [options] 1.bin 2.bin ...";
	cerr<<"\n	-d <density>	relative SNP density to Sanger sequencing (1)";
	cerr<<"\n	-b <burn>	burn-in generations (56)";
	cerr<<"\n	-l <file>	list of input files";
	cerr<<"\n	-m <mcmc>	sampling generations (200)";
	cerr<<"\n	-n <fold>	sample size*fold of nested MH sampler iteration (3)";
	cerr<<"\n	-t <thread>	number of threads (0=MAX)";
	cerr<<"\n	-v <vcf>	integrate known genotype in VCF format";
	cerr<<"\n	-r <file>	use recombination rate";	
	cerr<<"\n\n";
	exit(0);
}

int	main(int	ac,	char	**av){
	Impute	lp;
	lp.bn=56;	lp.sn=200;	lp.nn=3;	lp.density=1.0;
	uint	threads=0;
	vector<string>	file;
	
	char	opt;
	while((opt=getopt(ac,	av,	"d:b:l:m:n:t:v:r:"))>=0){
		switch(opt) {
		case	'd':	lp.density=atof(optarg);	break;
		case	'b':	lp.bn=atoi(optarg);	break;
		case	'm':	lp.sn=atoi(optarg);	break;
		case	'n':	lp.nn=atoi(optarg);	break;
		case	't':	threads=atoi(optarg);	break;
		case	'v':	lp.vcf_file.push_back(optarg);	break;
		case	'r':{	
			if(!lp.load_rec(optarg))	cerr<<"fail to load "<<optarg<<endl;	
			break;
		}	
		case	'l':{
			char	temp[256];
			FILE	*f=fopen(optarg,	"rt");
			while(fscanf(f,	"%s",	temp)!=EOF)	file.push_back(temp);
			fclose(f);
			}
		break;
	
		default:	lp.document();
		}
	}
	if(threads)	omp_set_num_threads(threads);
	for(int	i=optind;	i<ac;	i++)	file.push_back(av[i]);
	sort(file.begin(),	file.end());
	uint	fn=unique(file.begin(),	file.end())-file.begin();
	if(!fn)	lp.document();
	
	for(uint	i=0;	i<fn;	i++){
		timeval	sta,	end;
		gettimeofday(&sta,	NULL);

		lp.load_bin(file[i].c_str());
		for(uint	j=0;	j<lp.vcf_file.size();	j++){
			cerr<<lp.vcf_file[j]<<'\t'<<lp.load_vcf(lp.vcf_file[j].c_str())<<endl;
		}	
		lp.initialize();
		lp.estimate();
		lp.save_vcf(file[i].c_str());
		char	temp[256];
		sprintf(temp,	"mv %s %s.ok",	file[i].c_str(),	file[i].c_str());
		gettimeofday(&end,	NULL);
		cerr<<"rename\t"<<system(temp)<<endl;
		cerr<<"time\t"<<end.tv_sec-sta.tv_sec+1e-6*(end.tv_usec-sta.tv_usec)<<endl<<endl;
	}
	return	0;
}

