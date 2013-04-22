#include	<algorithm>
#include	<iostream>
#include	<dirent.h>
#include	<fstream>
#include	<sstream>
#include	"ebd.cpp"
#include	<cmath>
#include	<map>

struct	Result{
	uint8_t	ref,	alt,	ind,	val;
	float	score;
};

class	varisite{
private:
	uint	in,	gn;
	vector<EBD>	ebd;
	vector<uint32_t>	site;
	vector<Result>	result;
	vector<uint8_t>	ref;
	vector<uint32_t>	group;
	
	void	evaluate(uint	P,	uint	I);
public:
	bool	indel;
	float	cutoff,	stable;
	string	chr,	vcf;
	uint32_t	mask;
	static	void	document(void);
	bool	load_file(const	char	*F,	const	char	*C);
	bool	load_site(const	char	*F,	const	char	*C);
	bool	load_ref(const	char	*F);
	void	estimate(void);
};

void	varisite::document(void){
	cerr<<"\nvarisite";
	cerr<<"\nallele selection and site scoring of variants";
	cerr<<"\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage	varisite [options] <ebd.list chr chr.fa>";
	cerr<<"\n	-i	consider indel allele (false)";	
	cerr<<"\n	-s <FLT>	significance cutoff (1.5)";
	cerr<<"\n	-v <site.vcf>	only scoring on give sites";
	cerr<<"\n	-g <bit_mask>	group files by fields specified with bit_mask (0)";
	cerr<<"\n	-l <FLT>	variance stabilizer (1)";
	cerr<<"\n\n";
	exit(0);
}

bool	varisite::load_file(const	char	*F,	const	char	*C){
	ifstream	fi(F);
	if(!fi)	return	false;
	EBD	e;
	for(fi>>e.file;	!fi.eof();	fi>>e.file)	ebd.push_back(e);
	fi.close();
	cerr<<"binding\t"<<(in=ebd.size())<<endl;
	for(uint	i=0;	i<in;	i++)	ebd[i].bind(C);
	if(mask==0){	gn=1;	group.assign(in,	0);	}
	else{
		map<string,	uint32_t>	m;
		vector<string>	key(in);
		for(uint	i=0;	i<in;	i++){
			size_t	j=ebd[i].file.find_last_of('/');
			string	buff,	temp;
			if(j==string::npos)	buff=ebd[i].file;
			else	buff.assign(ebd[i].file,	j+1,	ebd[i].file.size()-j-1);	
			istringstream	si(buff);	j=0;
			while(!si.eof()){	getline(si,	temp,	'.');	if(mask&(1<<j))	key[i]+=temp;	j++;	}
			m[key[i]]++;
		}
		gn=0;	for(map<string,	uint32_t>::iterator	mi=m.begin();	mi!=m.end();	++mi,	gn++)	
			mi->second=gn;
		group.resize(in);
		for(uint	i=0;	i<in;	i++)	group[i]=m[key[i]];
	}
	cerr<<"groups\t"<<gn<<endl;
	return	true;
}

bool	varisite::load_site(const	char	*F,	const	char	*C){
	site.clear();
	FILE	*f=fopen(F,	"rt");
	if(f==NULL)	return	false;
	char	buff[1024],	temp[1024];
	uint32_t	pos;
	while(fgets(buff,	1024,	f)!=NULL){
		temp[0]=0;
		sscanf(buff,	"%s%u",	temp,	&pos);
		if(!strcmp(temp,	C))	site.push_back(pos);
	}
	fclose(f);
	sort(site.begin(),	site.end());
	cout<<"sites\t"<<site.size()<<endl;
	return	true;
}

bool	varisite::load_ref(const	char	*F){
	uint8_t	code[256];
	for(uint	i=0;	i<256;	i++)	code[i]=4;
	code['A']=code['a']=0;	code['C']=code['c']=1;
	code['G']=code['g']=2;	code['T']=code['t']=3;
		
	ref.clear();
	FILE	*f=fopen(F,	"rt");
	if(f==NULL)	return	false;
	cerr<<"loading\t"<<F<<endl;	
	char	buff[1024];
	while(fgets(buff,	1024,	f)!=NULL)	if(buff[0]!='>'){
		for(char	*p=buff;	*p>=33;	p++)	ref.push_back(*(code+*p));
	}	
	fclose(f);
	return	true;
}

void	varisite::evaluate(uint	P,	uint	I){
	if(site.size()){
		vector<uint32_t>::iterator	vi=lower_bound(site.begin(),	site.end(),	P+I+1);
		if(*vi!=P+I+1){	result[I].val=0;	return;	}
	}
	uint8_t	ref_all=ref[P+I],	alt_all=ref_all?0:1;
	if(ref_all==4){	result[I].val=0;	return;	}
	//	calculate ebd2
	float	ebd2[6]={0,0,0,0,0,0};
	uint	p=I*6;
	for(uint	i=0;	i<in;	i++){
		float	x;
		for(uint	j=0;	j<6;	j++){	x=0.1f*ebd[i].cnt[p+j];	ebd2[j]+=x*x;	}
	}
	//	choose allele
	float	best=-1;
	for(uint	j=0;	j<(indel?6:4);	j++)	if(ebd2[j]>best&&j!=ref_all){	best=ebd2[j];	alt_all=j;	}
	if(best<1){	result[I].val=0;	return;	}

	vector<float>	na(gn),	nb(gn),	f(gn);
	float	v0=0,	v1=0,	v2=0;
	for(uint	i=0;	i<in;	i++){
		na[group[i]]+=0.1f*ebd[i].cnt[p+ref_all];
		nb[group[i]]+=0.1f*ebd[i].cnt[p+alt_all];
	}
	for(uint	i=0;	i<gn;	i++)	if(na[i]+nb[i]){
		f[i]=nb[i]/(na[i]+nb[i]);	v0+=na[i]*nb[i]/(na[i]+nb[i]);
	}	
	for(uint	i=0;	i<in;	i++){
		float	a=0.1f*ebd[i].cnt[p+ref_all],	b=0.1f*ebd[i].cnt[p+alt_all],	e=f[group[i]]*(a+b);
		v1+=(b-e)*(b-e);
		v2+=fminf(fminf(a*a,	b*b),	0.25*(a-b)*(a-b));
	}
	result[I].ref=ref_all;	result[I].alt=alt_all;	result[I].ind=alt_all>=4?1:0;	result[I].val=1;
	result[I].score=(v1-v0+stable)/(v2+stable);
}

void	varisite::estimate(void){
	cerr.precision(1);	cerr.setf(ios::fixed);
	result.resize(65536);	string	s=chr+".sites.vcf";
	
	ofstream	fo(s.c_str());	s="ACGT+-";
	for(uint	p=0;	p<ref.size();	p+=65536){
		uint	len=p+65536<ref.size()?65536:ref.size()-p;
		#pragma omp parallel for
		for(uint	i=0;	i<in;	i++)
			if(!ebd[i].load(p>>16))	cerr<<"fail to read "<<ebd[i].file<<endl;
		#pragma omp parallel for	
		for(uint	i=0;	i<len;	i++)	evaluate(p,	i);
		for(uint	i=0;	i<len;	i++)	if(result[i].val){
			if(result[i].score<cutoff)	continue;
			fo<<chr<<'\t'<<p+i+1<<"\t.\t"<<s[result[i].ref]<<'\t'<<s[result[i].alt];
			fo<<'\t'<<result[i].score<<"\tPASS\t"<<(result[i].ind?"INDEL\n":"SNP\n");
		}
		cerr<<1e-6*p<<"M\r";	
	}
	fo.close();
}

int	main(int	ac,	char	**av){
	varisite	vs;	vs.cutoff=1.5;	vs.indel=false;	vs.mask=0;	vs.stable=1;

	int	opt;
	while((opt=getopt(ac,	av,	"is:v:g:l:"))>=0){
		switch(opt) {
		case	'i':	vs.indel=true;	break;
		case	's':	vs.cutoff=atof(optarg);	break;
		case	'v':	vs.vcf=optarg;	break;
		case	'g':	vs.mask=atoi(optarg);	break;
		case	'l':	vs.stable=atof(optarg);	break;
		default:	vs.document();
		}
	}
	if(ac<optind+3)	vs.document();	
	vs.chr=av[optind+1];
	if(vs.vcf.size())	vs.load_site(vs.vcf.c_str(),	av[optind+1]);
	vs.load_file(av[optind],	av[optind+1]);
	vs.load_ref(av[optind+2]);
	vs.estimate();	
	return	0;
}

