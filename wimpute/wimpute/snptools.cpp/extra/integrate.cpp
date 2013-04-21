/**
@file	integrate.cpp
@brief	Integration of multiple data source
@author	Yi Wang
@date	03/30/2011
*/
#include	<iostream>
#include	<stdint.h>
#include	<cstdlib>
#include	<fstream>
#include	<sstream>
#include	<tabix.h>
#include	<string>
#include	<vector>
#include	<zlib.h>
#include	<cmath>
#include	<map>
using	namespace	std;
typedef	uint32_t	uint;
typedef	__uint128_t	quad;
const	float	miss=1.0f/3.0f;

struct	Site{
	char	all[4];
	uint	pos;
	char	chr[8];
	bool	operator()(Site	X,	Site	Y){	return	*((quad*)&X)<*((quad*)&Y);	}
}__attribute__	((aligned(16)));

struct	Data{
	string	name;	//	the name of the dataset
	string	type;	//	geno for genotype data in vcf format; prob for probability data
	double	conf;	//	confidence of the data
	string	file;	//	file or URL
	
	uint	in,	mn;
	vector<string>	samp;
	vector<Site>	site;
	vector<float>	prob;
	bool	load_geno(const	char	*R);
	bool	load_prob(const	char	*R);
};

bool	Data::load_geno(const	char	*R){
	cerr<<"reading\t"<<file<<endl;
	double	rest=0.5*(1-conf);
	tabix_t *t;	ti_iter_t	iter;	const	char	*s;	int	len;
	if((t=ti_open(file.c_str(), 0))==0)	return	false;
	{
		iter=ti_query(t, 0, 0, 0);
		while((s=ti_read(t,	iter,	&len))!=0){
			if(*s!='#') break;
			else	if(strstr(s,	"#CHROM")==NULL)	continue;
			istringstream	si(s);
			string	a;
			si>>a>>a>>a>>a>>a>>a>>a>>a>>a;
			while(!si.eof()){	si>>a;	samp.push_back(a);	}
		}
		ti_iter_destroy(iter);
		in=samp.size();
	}
	iter=ti_querys(t,	R);
	while((s=ti_read(t,	iter,	&len))!=0){
		istringstream	si(s);	string	a;
		Site	ts;	memset(&ts,	0,	sizeof(Site));
		si>>a;	strncpy(ts.chr,	a.c_str(),	7);
		si>>ts.pos>>a>>a;	ts.all[0]=a.size()==1?a[0]:'R';
		si>>a;	ts.all[1]=a.size()==1?a[0]:'V';
		si>>a>>a>>a>>a;
		site.push_back(ts);
		for(uint	i=0;	i<in;	i++){
			si>>a;
			if(a[0]=='.'||a[2]=='.'){	prob.push_back(miss);	prob.push_back(miss);	prob.push_back(miss);	}
			else	if(a[0]=='0'&&a[2]=='0'){	prob.push_back(conf);	prob.push_back(rest);	prob.push_back(rest);	}
			else	if(a[0]=='0'&&a[0]!=a[2]){	prob.push_back(rest);	prob.push_back(conf);	prob.push_back(rest);	}
			else{	prob.push_back(rest);	prob.push_back(rest);	prob.push_back(conf);	}
		}
	}
	mn=site.size();
	ti_iter_destroy(iter);
	ti_close(t);
	cerr<<"sites\t"<<mn<<"\nsample\t"<<in<<endl<<endl;
	return	true;
}

bool	Data::load_prob(const	char	*R){
	cerr<<"reading\t"<<file<<endl;
	tabix_t *t;	ti_iter_t	iter;	const	char	*s;	int	len;
	if((t=ti_open(file.c_str(), 0))==0)	return	false;
	{
		iter=ti_query(t, 0, 0, 0);
		s=ti_read(t,	iter,	&len);
		istringstream	si(s);
		string	a;
		si>>a>>a>>a;
		while(!si.eof()){	si>>a;	samp.push_back(a);	}
		ti_iter_destroy(iter);
		in=samp.size();
	}
	iter=ti_querys(t,	R);
	while((s=ti_read(t,	iter,	&len))!=0){
		istringstream	si(s);	string	a;
		Site	ts;	memset(&ts,	0,	sizeof(Site));
		si>>a;	strncpy(ts.chr,	a.c_str(),	7);
		si>>ts.pos>>a;	ts.all[0]=a[0];	ts.all[1]=a[1];
		site.push_back(ts);
		for(uint	i=0;	i<in;	i++){
			float	a,	h,	b;
			si>>h>>b;	a=fmaxf(1-h-b,0);
			a=powf(a,	conf);	h=powf(h,	conf);	b=powf(b,	conf);
			float	sum=1.0f/(a+h+b);
			prob.push_back(a*sum);	prob.push_back(h*sum);	prob.push_back(b*sum);
		}
	}
	mn=site.size();
	ti_iter_destroy(iter);
	ti_close(t);
	cerr<<"sites\t"<<mn<<"\nsample\t"<<in<<endl<<endl;
	return	true;
}
	
class	Integrate{
private:
	uint	dn,	in,	mn;
	vector<Data>	data;
public:
	static	void	document(void);	

	string	chr;
	uint	sta,	end;
	bool	load_list(const	char	*F);
	bool	fetch_data(void);
	void	make_union(void);
};

void	Integrate::document(void){
	cerr<<"\nIntegration of multiple data source";
	cerr<<"\nauthor\tYi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage\tintegrate integrate <datasets.list chr begin end>";
	cerr<<"\n\n";
	exit(0);
}

bool	Integrate::load_list(const	char	*F){
	ifstream	fi(F);
	if(!fi)	return	false;
	Data	temp;
	fi>>temp.file>>temp.file>>temp.file>>temp.file;
	for(fi>>temp.name>>temp.type>>temp.conf>>temp.file;	!fi.eof();	fi>>temp.name>>temp.type>>temp.conf>>temp.file){
		if(temp.name.size()&&temp.name[0]!='#')	data.push_back(temp);
	}	
	fi.close();
	return	true;
}

bool	Integrate::fetch_data(void){
	char	region[256];
	sprintf(region,	"%s:%u-%u",	chr.c_str(),	sta,	end);
	dn=data.size();
	for(uint	i=0;	i<dn;	i++){
		if(data[i].type=="geno"){
			if(!data[i].load_geno(region))	return	false;
		}	
		else	if(data[i].type=="prob"){
			if(!data[i].load_prob(region))	return	false;
		}
		else{	cerr<<"invalid data type:\t"<<data[i].name<<':'<<data[i].type<<endl;	return	false;	}
	}
	return	true;
}

void	Integrate::make_union(void){
	cerr<<"integrating...\n";	
	map<string,	uint>	imap;
	map<Site,	uint,	Site>	smap;
	for(uint	d=0;	d<dn;	d++){
		for(uint	i=0;	i<data[d].in;	i++)	imap[data[d].samp[i]]++;
		for(uint	m=0;	m<data[d].mn;	m++)	smap[data[d].site[m]]++;
	}
	in=mn=0;
	for(map<string,	uint>::iterator	mi=imap.begin();	mi!=imap.end();	++mi,	in++)	mi->second=in;
	for(map<Site,	uint,	Site>::iterator	mi=smap.begin();	mi!=smap.end();	++mi,	mn++)	mi->second=mn;
	vector<float>	prob;	prob.assign(in*mn*3,	1);
	cerr<<"sites\t"<<mn<<"\nsample\t"<<in<<endl<<endl;
	
	for(uint	n=0;	n<dn;	n++){
		Data	&d=data[n];
		vector<uint>	id(d.in);
		for(uint	i=0;	i<d.in;	i++)	id[i]=imap[d.samp[i]];
		for(uint	m=0;	m<d.mn;	m++){
			float	*p=&d.prob[m*d.in*3],	*q=&prob[smap[d.site[m]]*in*3];
			for(uint	i=0;	i<d.in;	i++,	p+=3){
				float	*s=q+id[i]*3;	s[0]*=p[0];	s[1]*=p[1];	s[2]*=p[2];
			}
		}
	}

	char	temp[256];
	sprintf(temp,	"%s_%09u_%09u.bin",	chr.c_str(),	sta,	end);
	gzFile	f=gzopen(temp,	"wt");
	gzprintf(f,	"#%u\tpos\tallele",	in);
	for(map<string,	uint>::iterator	mi=imap.begin();	mi!=imap.end();	++mi)	gzprintf(f,	"\t%s",	mi->first.c_str());
	mn=0;
	for(map<Site,	uint,	Site>::iterator	mi=smap.begin();	mi!=smap.end();	++mi,	mn++){
		gzprintf(f,	"\n%s\t%u\t%c%c",	mi->first.chr,	mi->first.pos,	mi->first.all[0],	mi->first.all[1]);
		float	*p=&prob[mn*in*3];
		for(uint	i=0;	i<in;	i++,	p+=3){
			double	sum=p[0]+p[1]+p[2];
			if(sum>0)	gzprintf(f,	"\t%g %g",	p[1]/sum,	p[2]/sum);
			else	gzprintf(f,	"\t%g %g",	miss,	miss);
		}
	}
	gzprintf(f,	"\n");
	gzclose(f);
}

int	main(int	ac,	char	**av){
	if(ac!=5)	Integrate::document();
	Integrate	ing;
	ing.chr=av[2];	ing.sta=atoi(av[3]);	ing.end=atoi(av[4]);
	if(!ing.load_list(av[1])){	cerr<<"fail to load "<<av[1]<<endl;	return	0;	}
	if(!ing.fetch_data()){	cerr<<"fail to fetch all data\n";	return	0;	}
	ing.make_union();
	return	0;
}
