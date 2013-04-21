/**
@file	bamodel.cpp
@brief	Three-component binomial mixture modeling of BAMs on one same individual
@author	Yi Wang
@date	03/29/2011
*/
#include	<algorithm>
#include	<iostream>
#include	<stdint.h>
#include	<cstring>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<cmath>
#include	<sam.h>
using	namespace	std;
typedef	uint32_t	uint;

struct	Site{
	char	chr[10],	all[2];
	uint	pos;
	bool	operator()(Site	X,	Site	Y){
		int	cmp=strcmp(X.chr,Y.chr);
		return	cmp?cmp<0:X.pos<Y.pos;
	}
};

struct	Site_Equal{
	bool	operator()(Site	X,	Site	Y){	return	X.pos==Y.pos&&!strcmp(X.chr,Y.chr);	}
};

struct	bam_file{
	string	file;	//	bam file name
	Site	*site;	//	pointer to site list
	char	**name;	//	chromosome names
	uint32_t	sn;	//	number of sites
	vector<float>	ebd;//	effective base depth
	
	static	float	phred[256];	//	phred[i]=1-P(Error[i])
	static int pileup_func(uint tid,	uint pos,	int n,	const	bam_pileup1_t	*pl,	void	*data);
	bool	scan_bam(void);
};

float	bam_file::phred[256];

int	bam_file::pileup_func(uint tid,	uint	pos,	int n,	const	bam_pileup1_t	*pl,	void	*data){
	bam_file	*bf=(bam_file*)data;
	Site	s;	memset(&s,	0,	sizeof(Site));
	strncpy(s.chr,	bf->name[tid],	7);	s.pos=pos+1;
	Site	*vi=lower_bound(bf->site,	bf->site+bf->sn,	s,	Site());
	if(s.pos!=vi->pos||strcmp(s.chr,vi->chr))	return	0;
	float	*pd=&(bf->ebd[(vi-bf->site)*2]);
	const	bam_pileup1_t	*p=pl;
	uint8_t	ref=*(bam_nt16_table+vi->all[0]),	alt=*(bam_nt16_table+vi->all[1]);
	
	if(vi->all[1]=='+'){
		for(int	i=0;	i<n;	i++,	p++){
			if(p->indel>0)	pd[1]+=phred[p->b->core.qual];
			else	if(!p->is_del&&bam1_seqi(bam1_seq(p->b),p->qpos)==ref)	pd[0]+=phred[bam1_qual(p->b)[p->qpos]]*phred[p->b->core.qual];
		}
	}
	else	if(vi->all[1]=='-'){
		for(int	i=0;	i<n;	i++,	p++){
			if(p->indel<0)	pd[1]+=phred[p->b->core.qual];
			else	if(!p->is_del&&bam1_seqi(bam1_seq(p->b),p->qpos)==ref)	pd[0]+=phred[bam1_qual(p->b)[p->qpos]]*phred[p->b->core.qual];
		}	
	}
	else{	
		for(int	i=0;	i<n;	i++,	p++)	if(!p->is_del&&p->b->core.qual){
			uint8_t	c=bam1_seqi(bam1_seq(p->b),	p->qpos);
			if(c==ref)	pd[0]+=phred[bam1_qual(p->b)[p->qpos]]*phred[p->b->core.qual];
			else	if(c==alt)	pd[1]+=phred[bam1_qual(p->b)[p->qpos]]*phred[p->b->core.qual];
		}
	}	
	return	0;	
}

bool	bam_file::scan_bam(void){
	samfile_t *bam;
	bam=samopen(file.c_str(), "rb", 0);
	if(bam==0)	return	false;
	cerr<<"scaning\t"<<file<<endl;
	name=bam->header->target_name;
	ebd.assign(sn*2,	0);
	sampileup(bam,	-1,	pileup_func,	this);
	samclose(bam);
	return	true;
}

class	bamodel{
private:
	uint	bn,	sn;
	vector<Site>	site;
	vector<uint16_t>	data;
public:
	static	double	prior_a,	prior_h,	prior_r,	precision;	//	priors for emission probabilities
	static	void	document(void);

	vector<bam_file>	bam;
	bool	load_site(const	char	*F);
	bool	scan_bams(void);
	void	estimate(const	char	*F,	bool	INDEL);
	void	save(const	char	*F);
};

double	bamodel::prior_a;
double	bamodel::prior_h;
double	bamodel::prior_r;
double	bamodel::precision;

bool	bamodel::load_site(const	char	*F){
	string	temp,	a,	b;
	ifstream	fi(F);	if(!fi)	return	false;
	cerr<<"loading\t"<<F<<endl;
	Site	s;
	for(getline(fi,	temp);	!fi.eof();	getline(fi,	temp))	if(temp.size()&&temp[0]!='#'){
		istringstream	si(temp);
		memset(&s,	0,	sizeof(Site));
		si>>a;	strncpy(s.chr,	a.c_str(),	9);
		si>>s.pos>>a>>a>>b;	s.all[0]=a[0];	s.all[1]=b[0];
		site.push_back(s);
	}
	fi.close();
	sort(site.begin(),	site.end(),	Site());
	sn=unique(site.begin(),	site.end(),	Site_Equal())-site.begin();
	data.resize(sn*2);
	cerr<<"unique\t"<<sn<<endl;
	return	true;
}

bool	bamodel::scan_bams(void){
	bn=bam.size();
	for(uint	i=0;	i<bn;	i++){
		bam[i].site=&site[0];	bam[i].sn=sn;	bam[i].ebd.resize(2*sn);
		if(!bam[i].scan_bam()){	cerr<<"fail to scan "<<bam[i].file<<endl;	return	false;	}
	}
	return	true;
}

void	bamodel::estimate(const	char	*F,	bool	INDEL){
	double	na=0.91,	nh=0.05,	nb=0.04,	oa,	ob,	oh,	delta;
	vector<double>	ea(bn),	eh(bn),	eb(bn);
	vector<double>	lap(bn),	laq(bn),	lhp(bn),	lhq(bn),	lbp(bn),	lbq(bn);
	vector<double>	sap(bn),	saq(bn),	shp(bn),	shq(bn),	sbp(bn),	sbq(bn);
	
	for(uint	i=0;	i<bn;	i++){	ea[i]=0.99;	eh[i]=0.5;	eb[i]=0.02;	}
	cout.precision(4);	cout.setf(ios::fixed);
	do{
		double	la=log(na),	lh=log(nh),	lb=log(nb);
		oa=na;	oh=nh;	ob=nb;	na=nh=nb=0;
		for(uint	i=0;	i<bn;	i++){
			lap[i]=log(ea[i]);	laq[i]=log(1-ea[i]);
			lhp[i]=log(eh[i]);	lhq[i]=log(1-eh[i]);
			lbp[i]=log(eb[i]);	lbq[i]=log(1-eb[i]);
			sap[i]=saq[i]=shp[i]=shq[i]=sbp[i]=sbq[i]=0;
		}
		for(uint	m=0;	m<sn;	m++)	if((site[m].all[1]=='+'||site[m].all[1]=='-')==INDEL){
			uint	p=m*2;
			double	wa=la,	wh=lh,	wb=lb;
			for(uint	i=0;	i<bn;	i++){
				double	ta=bam[i].ebd[p],	tb=bam[i].ebd[p+1];
				wa+=ta*lap[i]+tb*laq[i];	wh+=ta*lhp[i]+tb*lhq[i];	wb+=ta*lbp[i]+tb*lbq[i];
			}
			if(wa>=wh&&wa>=wb){	wh=exp(wh-wa);	wb=exp(wb-wa);	wa=1;	}
			else	if(wh>=wa&&wh>=wb){	wa=exp(wa-wh);	wb=exp(wb-wh);	wh=1;	}
			else{	wa=exp(wa-wb);	wh=exp(wh-wb);	wb=1;	}
			double	sum=1.0/(wa+wh+wb);	wa*=sum;	wh*=sum;	wb*=sum;
			na+=wa;	nh+=wh;	nb+=wb;
			for(uint	i=0;	i<bn;	i++){
				double	ta=bam[i].ebd[p],	tb=bam[i].ebd[p+1];
				sap[i]+=ta*wa;	saq[i]+=tb*wa;	shp[i]+=ta*wh;	shq[i]+=tb*wh;	sbp[i]+=ta*wb;	sbq[i]+=tb*wb;
			}
		}
		double	sum=1.0/(na+nh+nb+1);
		na=(na+0.91)*sum;	nh=(nh+0.05)*sum;	nb=(nb+0.04)*sum;
		cout<<"Ref/Ref\tRef/Alt\tAlt/Alt\tFile\n";
		for(uint	i=0;	i<bn;	i++){
			ea[i]=(sap[i]+prior_r*precision)/(sap[i]+saq[i]+precision);
			eh[i]=(shp[i]+prior_h*precision)/(shp[i]+shq[i]+precision);
			eb[i]=(sbp[i]+prior_a*precision)/(sbp[i]+sbq[i]+precision);
			cout<<ea[i]<<'\t'<<eh[i]<<'\t'<<eb[i]<<'\t'<<bam[i].file<<endl;
		}
		cout<<na<<'\t'<<nh<<'\t'<<nb<<'\t'<<F<<endl<<endl;
		delta=sqrt((na-oa)*(na-oa)+(nh-oh)*(nh-oh)+(nb-ob)*(nb-ob));
	}while(delta>1e-6);
	
	string	s=F;	s+=".emit";
	FILE	*f=fopen(s.c_str(),	INDEL?"at":"wt");
	for(uint	i=0;	i<bn;	i++)	fprintf(f,	"%s\t%g\t%g\t%g\n",	bam[i].file.c_str(),	ea[i],	eh[i],	eb[i]);
	fclose(f);
	
	s=F;	s+=".weig";
	f=fopen(s.c_str(),	INDEL?"at":"wt");
	fprintf(f,	"%s\t%g\t%g\t%g\n",	F,	na,	nh,	nb);
	fclose(f);
	
	for(uint	m=0;	m<sn;	m++)	if((site[m].all[1]=='+'||site[m].all[1]=='-')==INDEL){
		uint	p=m*2;
		double	wa=0,	wh=0,	wb=0;
		for(uint	i=0;	i<bn;	i++){
			double	ta=bam[i].ebd[p],	tb=bam[i].ebd[p+1];
			wa+=ta*lap[i]+tb*laq[i];	wh+=ta*lhp[i]+tb*lhq[i];	wb+=ta*lbp[i]+tb*lbq[i];
		}
		if(wa>=wh&&wa>=wb){	wh=exp(wh-wa);	wb=exp(wb-wa);	wa=1;	}
		else	if(wh>=wa&&wh>=wb){	wa=exp(wa-wh);	wb=exp(wb-wh);	wh=1;	}
		else{	wa=exp(wa-wb);	wh=exp(wh-wb);	wb=1;	}
		double	sum=50000/(wa+wh+wb);
		data[p]=(uint16_t)(wh*sum+0.5);	data[p+1]=(uint16_t)(wb*sum+0.5);
	}
}

void	bamodel::save(const	char	*F){
	string	s=F;	s+=".raw";	
	FILE	*f=fopen(s.c_str(),	"wb");
	fwrite(&data[0],	4*sn,	1,	f);
	fclose(f);
}

void	bamodel::document(void){
	cerr<<"\nbamodel";
	cerr<<"\nthree-component binomial mixture modeling of BAMs on the SAME individual";
	cerr<<"\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage	bamodel [options] <out site.vcf in1.bam> [in2.bam in3.bam...]";
	cerr<<"\n	-a <FLOAT>	P(read=ref|geno=alt/alt) (0.010)";
	cerr<<"\n	-h <FLOAT>	P(read=ref|geno=ref/alt) (0.500)";
	cerr<<"\n	-r <FLOAT>	P(read=ref|geno=ref/ref) (0.995)";
	cerr<<"\n	-p <FLOAT>	precision of beta distribution prior (100)";
	cerr<<"\n\n";
	exit(0);	
}

int	main(int	ac,	char	**av){
	for(uint	i=0;	i<256;	i++)	bam_file::phred[i]=1-exp(-0.1*i);
	bamodel::prior_a=0.010;	bamodel::prior_h=0.500;	bamodel::prior_r=0.995;	bamodel::precision=100;
	
	int	opt;
	while((opt=getopt(ac,	av,	"a:h:r:p:"))>=0){
		switch(opt){
		case	'a':	bamodel::prior_a=atof(optarg);	break;
		case	'h':	bamodel::prior_h=atof(optarg);	break;
		case	'r':	bamodel::prior_r=atof(optarg);	break;
		case	'p':	bamodel::precision=atof(optarg);break;
		default:	bamodel::document();
		}
	}
	if(ac<optind+3)	bamodel::document();
	
	bamodel	bm;	bm.bam.resize(ac-optind-2);
	for(int	i=optind+2;	i<ac;	i++)	bm.bam[i-optind-2].file=av[i];
	if(!bm.load_site(av[optind+1])){	cerr<<"fail to load\t"<<av[optind]<<endl;	return	0;	}
	if(!bm.scan_bams())	return	0;
	bm.estimate(av[optind],	false);
	bm.estimate(av[optind],	true);
	bm.save(av[optind]);
	return	0;
}
