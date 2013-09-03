/**
@file	hapfuse.cpp
@brief	joint chunked haplotypes into chromosome wide haplotypes
@author	Yi Wang
@date	04/01/2011
*/
#include	<iostream>
#include	<algorithm>
#include	<dirent.h>
#include	<stdint.h>
#include	<cstring>
#include	<fstream>
#include	<sstream>
#include	<cstdio>
#include	<string>
#include	<vector>
#include	<zlib.h>
#include	<list>
#include	<set>
using	namespace	std;
typedef	double	real;

struct	Site{
	vector<float>	hap;
	string	chr;
	uint32_t	pos;
	uint16_t	cov;
	char	all[2];
	
	static	bool	is_x;
	static	vector<bool>	is_male;

	void	write(FILE	*F);
};

bool	Site::is_x;
vector<bool>	Site::is_male;
	
void	Site::write(FILE	*F){
	fprintf(F,	"%s\t%u\t.\t%c\t%c\t100\tPASS\t.\tGT:AP",	chr.c_str(),	pos,	all[0],	all[1]);
	uint	in=hap.size()/2;
	real	k=1.0/cov;
	bool	is_par=(pos>=60001&&pos<=2699520)||(pos>=154931044&&pos<=155270560);
	for(uint	i=0;	i<in;	i++){
		uint	a,	b;
		real	p0=hap[i*2]*k,	p1=hap[i*2+1]*k;
		real	prr=(1-p0)*(1-p1),	pra=(1-p0)*p1+p0*(1-p1),	paa=p0*p1;
		if(is_x&&is_male[i]&&!is_par){
			if(prr>=paa)	a=b=0;
			else	a=b=1;	
		}
		else{
			if(prr>=pra&&prr>=paa)	a=b=0;
			else	if(pra>=prr&&pra>=paa){	
				if(p0>p1){	a=1;	b=0;	}
				else{	a=0;	b=1;	}
			}
			else	a=b=1;
		}		
		fprintf(F,	"\t%u|%u:%.3f,%.3f",	a,	b,	p0,	p1);
	}
	fprintf(F,	"\n");
}
	
class	hapfuse{
private:
	list<Site>	site;
	uint32_t	in;
	FILE	*vcf;
	vector<string>	file;
	vector<string>	name;
	vector<Site>	chunk;
	set<string>	male;

	void	vcf_head(const	char	*F);
	bool	load_chunk(const	char	*F);
public:
	static	void	document(void);
	bool	gender(char	*F);
	bool	load_dir(const	char	*D);
	void	work(const	char	*F);
};

bool	hapfuse::gender(char	*F){
	ifstream	fi(F);
	if(!fi){	cerr<<"fail to open "<<F<<endl;	return	false;	}
	string	s,	g;
	for(fi>>s>>g;	!fi.eof();	fi>>s>>g)	if(g=="male")	male.insert(s);
	fi.close();
	return	true;
}

bool	hapfuse::load_chunk(const	char	*F){
	char	buff[65536];
	string	temp;
	name.clear();	chunk.clear();
	gzFile	f=gzopen(F,	"rt");
	if(f==Z_NULL)	return	false;
	do	gzgets(f,	buff,	65536);	while(strstr(buff,	"#CHROM")==NULL);
	{
		istringstream	si(buff);
		for(uint	i=0;	i<9;	i++)	si>>temp;
		temp="";
		while(!si.eof()){
			si>>temp;	if(temp!="")	name.push_back(temp);
			temp="";
		}
	}
	
	in=name.size();
	Site::is_male.resize(in);
	for(uint	i=0;	i<in;	i++)	Site::is_male[i]=(male.find(name[i])!=male.end());
	
	Site	s;	s.hap.resize(in*2);	s.cov=1;
	string	a,	b;
	while(!gzeof(f)){
		gzgets(f,	buff,	65536);
		istringstream	si(buff);
		si>>s.chr>>s.pos>>a>>a>>b;
		s.all[0]=a[0];	s.all[1]=b[0];
		si>>a>>a>>a>>a;
		for(uint	i=0;	i<in;	i++){
			si>>a;	a[9]=' ';
			sscanf(a.c_str()+4,	"%f%f",	&s.hap[i*2],	&s.hap[i*2+1]);
		}
		chunk.push_back(s);
	}
	gzclose(f);	
	return	true;
}
	
void	hapfuse::vcf_head(const	char	*F){
	vcf=fopen(F,	"wt");
	fprintf(vcf,	"##fileformat=VCFv4.0\n");
	fprintf(vcf,	"##source=BCM:SNPTools:hapfuse\n");
	fprintf(vcf,	"##reference=1000Genomes-NCBI37\n");
	fprintf(vcf,	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	fprintf(vcf,	"##FORMAT=<ID=AP,Number=2,Type=Float,Description=\"Allelic Probability, P(Allele=1|Haplotype)\">\n");
	fprintf(vcf,	"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
}

bool	hapfuse::load_dir(const	char	*D){
	string	d=D;	if(d.size()&&d[d.size()-1]!='/')	d+='/';
	DIR	*dir=opendir(D);
	if(dir==NULL){	cerr<<"fail to open "<<D<<endl;	return	false;	}
	struct	dirent  *ptr;
	while((ptr=readdir(dir))!=NULL){
		string	s=d+ptr->d_name;
		if(s.find(".vcf.gz")!=string::npos)	file.push_back(s);
	}
	closedir(dir);
	sort(file.begin(),	file.end());
	return	true;
}
	
void	hapfuse::work(const	char	*F){
	vcf_head(F);
	vector<real>	sum;
	for(uint	i=0;	i<file.size();	i++){
		if(!load_chunk(file[i].c_str()))	return;
		if(!i){
			for(uint	i=0;	i<in;	i++)	fprintf(vcf,	"\t%s",	name[i].c_str());
			fprintf(vcf,	"\n");
		}
		while(!site.empty()&&site.front().pos<chunk[0].pos){
			site.front().write(vcf);
			site.pop_front();
		}
		//	find the correct phase
		sum.assign(in*2,	0);
		for(list<Site>::iterator	li=site.begin();	li!=site.end();	li++)	for(uint	m=0;	m<chunk.size();	m++)	if(li->pos==chunk[m].pos){
			float	*p=&(li->hap[0]),	*q=&(chunk[m].hap[0]);
			for(uint	j=0;	j<in;	j++){
				sum[j*2]+=p[j*2]*q[j*2]+p[j*2+1]*q[j*2+1];	sum[j*2+1]+=p[j*2]*q[j*2+1]+p[j*2+1]*q[j*2];
			}			
		}
		//	swap phase if needed
		for(uint	j=0;	j<in;	j++)	if(sum[j*2]<sum[j*2+1]){
			for(uint	m=0;	m<chunk.size();	m++){
				float	t=chunk[m].hap[j*2];	chunk[m].hap[j*2]=chunk[m].hap[j*2+1];	chunk[m].hap[j*2+1]=t;
			}	
		}
		//	add chunk to buffer		
		for(uint	m=0;	m<chunk.size();	m++){
			bool	found=false;
			for(list<Site>::iterator	li=site.begin();	li!=site.end();	li++)	if(li->pos==chunk[m].pos){
				found=true;	li->cov++;
				float	*p=&(li->hap[0]),	*q=&(chunk[m].hap[0]);
				for(uint	j=0;	j<in*2;	j++)	p[j]+=q[j];
			}
			if(!found)	site.push_back(chunk[m]);
		}
		cerr<<file[i]<<endl;
	}
	for(list<Site>::iterator	li=site.begin();	li!=site.end();	li++)	li->write(vcf);
	fclose(vcf);
}

void	hapfuse::document(void){
	cerr<<"\nhapfuse";
	cerr<<"\njoint chunked haplotypes into chromosome wide haplotypes";
	cerr<<"\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage	hapfuse <out.vcf dir> [gender]";
	cerr<<"\n\n";
	exit(0);		
}

int	main(int	ac,	char	**av){
	if(ac<3)	hapfuse::document();
	hapfuse	hf;
	if(ac<4)	Site::is_x=false;
	else{	Site::is_x=true;	hf.gender(av[3]);	}
	hf.load_dir(av[2]);
	hf.work(av[1]);
	return	0;
}

