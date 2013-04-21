/**
@file	realhap.cpp
@brief	query single/paired read haplotype on given sites
@author	Yi Wang
@date	04/10/2011
*/
#include	<algorithm>
#include	<iostream>
#include	<stdint.h>
#include	<cstring>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<zlib.h>
#include	<cmath>
#include	<sam.h>
using	namespace	std;
typedef	uint32_t	uint;
typedef	__uint128_t	quad;

struct	Site{
	char	all[4];
	uint	pos;
	char	chr[8];
	bool	operator()(Site	X,	Site	Y){	return	*((quad*)&X)<*((quad*)&Y);	}
}__attribute__	((aligned(16)));

struct	Haplo{
	char	rid[28];
	uint32_t	sid;
	bool	operator()(Haplo	X,	Haplo	Y){
		int	cmp=strcmp(X.rid,	Y.rid);
		return	cmp?cmp<0:X.sid<Y.sid;
	}
};

struct	Equal{
	bool	operator()(Haplo	X,	Haplo	Y){
		return	(X.sid>>1)==(Y.sid>>1)&&!strcmp(X.rid,	Y.rid);
	}
};

class	realhap{
public:
	static	float	cutoff;
	static	float	phred[256];	//	phred[i]=1-P(Error[i])
	static	void	document(void);
	static int pileup_func(uint tid,	uint pos,	int n,	const	bam_pileup1_t	*pl,	void	*data);
		
	uint	sn;
	vector<Site>	site;
	char	**name;
	vector<Haplo>	haplo;
		
	bool	load_site(const	char	*F);
	bool	scan_bam(const	char	*F);
	void	save_hap(const	char	*F);
};

float	realhap::cutoff;
float	realhap::phred[256];

void	realhap::document(void){
	cerr<<"\nrealhap";
	cerr<<"\nquery single/paired read haplotype on given sites";
	cerr<<"\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage	realhap [options] <site.vcf in.bam>";
	cerr<<"\n	-c <FLT>	recalibrated base quality cutoff (0.9)";
	cerr<<"\n\n";
	exit(0);	
}

bool	realhap::load_site(const	char	*F){
	string	temp,	a,	b;
	ifstream	fi(F);	if(!fi)	return	false;
	cerr<<"loading\t"<<F<<endl;
	Site	s;
	for(getline(fi,	temp);	!fi.eof();	getline(fi,	temp))	if(temp.size()&&temp[0]!='#'){
		istringstream	si(temp);
		memset(&s,	0,	sizeof(Site));
		si>>a;	strncpy(s.chr,	a.c_str(),	7);
		si>>s.pos>>a>>a>>b;	s.all[0]=a[0];	s.all[1]=b[0];
		site.push_back(s);
	}
	fi.close();
	cerr<<"sorting "<<(sn=site.size())<<"\tsites\n";
	sort(site.begin(),	site.end(),	Site());
	return	true;
}

bool	realhap::scan_bam(const	char	*F){
	samfile_t *bam;
	bam=samopen(F, "rb", 0);
	if(bam==0)	return	false;
	cerr<<"scaning\t"<<F<<endl;
	name=bam->header->target_name;
	sampileup(bam,	-1,	pileup_func,	this);
	samclose(bam);
	return	true;
}

int	realhap::pileup_func(uint tid,	uint	pos,	int n,	const	bam_pileup1_t	*pl,	void	*data){
	realhap	*rh=(realhap*)data;
	Site	s;	memset(&s,	0,	sizeof(Site));
	strncpy(s.chr,	rh->name[tid],	7);	s.pos=pos+1;
	vector<Site>::iterator	vi=lower_bound(rh->site.begin(),	rh->site.end(),	s,	Site());
	if(s.pos!=vi->pos||strcmp(s.chr,vi->chr))	return	0;
	Haplo	h;	h.sid=(vi-rh->site.begin())*2;
	uint8_t	ref=*(bam_nt16_table+vi->all[0]),	alt=*(bam_nt16_table+vi->all[1]);
	const	bam_pileup1_t	*p=pl;
	for(int	i=0;	i<n;	i++,	p++)	if(!p->is_del&&p->b->core.qual){
		uint8_t	c=bam1_seqi(bam1_seq(p->b),	p->qpos);
		if(phred[bam1_qual(p->b)[p->qpos]]*phred[p->b->core.qual]<cutoff)	continue;
		if(c==ref){
			strcpy(h.rid,	bam1_qname(p->b));	h.sid&=0xFFFFFFFE;	rh->haplo.push_back(h);
		}	
		else	if(c==alt){
			strcpy(h.rid,	bam1_qname(p->b));	h.sid|=0x1;	rh->haplo.push_back(h);
		}	
	}
	return	0;	
}

void	realhap::save_hap(const	char	*F){
	string	temp=F;	temp+=".hap";
	cerr<<"saving\t"<<temp<<endl;
	sort(haplo.begin(),	haplo.end(),	Haplo());
	uint	un=unique(haplo.begin(),	haplo.end(),	Equal())-haplo.begin();	
	gzFile	f=gzopen(temp.c_str(),	"wt");
	for(uint	i=0;	i<un;	){
		uint	j;
		for(j=i+1;	j<un&&!strcmp(haplo[i].rid,haplo[j].rid);	j++);
		if(j-i>1){
			uint	s=haplo[i].sid>>1;
			gzprintf(f,	"%s\t%s\n",	haplo[i].rid,	site[s].chr);
			for(uint	k=i;	k<j;	k++){
				uint	t=haplo[k].sid>>1;
				if(strcmp(site[s].chr,	site[t].chr))	continue;
				gzprintf(f,	"\t%u\t%c\n",	site[t].pos,	site[t].all[haplo[k].sid&1]);
			}
		}
		i=j;
	}
	gzclose(f);
}

int	main(int	ac,	char	**av){
	for(uint	i=0;	i<256;	i++)	realhap::phred[i]=1-exp(-0.1*i);
	realhap::cutoff=0.9;
	
	char	opt;
	while((opt=getopt(ac,	av,	"c:"))>=0){
		switch(opt){
		case	'c':	realhap::cutoff=atof(optarg);	break;
		default:	realhap::document();
		}
	}
	if(ac<optind+2)	realhap::document();
	
	realhap	rh;
	if(!rh.load_site(av[optind])){	cerr<<"fail to load\t"<<av[optind]<<endl;	return	0;	}
	if(!rh.scan_bam(av[optind+1])){	cerr<<"fail to load\t"<<av[optind+1]<<endl;	return	0;	}
	rh.save_hap(av[optind+1]);
	return	0;
}
