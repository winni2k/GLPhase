/**
@file	probin.cpp
@brief	generate chunked .bin files for each chromosome
@author	Yi Wang
@date	04/01/2011
*/
#include	<iostream>
#include	<stdint.h>
#include	<cstring>
#include	<cstdlib>
#include	<cstdio>
#include	<string>
#include	<vector>
#include	<zlib.h>
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

class	probin{
private:
	uint64_t	offset;
	uint32_t	mn,	in;
	vector<Site>	site;
	vector<string>	name;

	void	load_head(FILE	*F);
	void	write_bin(FILE	*F,	uint64_t	S,	uint	N);
public:
	static	void	document(void);	
	
	string	chr,	fld;
	uint	bin,	str;   // size of bin and stride of binning
	
	void	work(const	char	*F);
};

void	probin::load_head(FILE	*F){
	if(!fread(&mn,	4,	1,	F))	return;  // comes from the prob files 
	if(!fread(&in,	4,	1,	F))	return;
	site.resize(mn);
	if(!fread(&site[0],	mn*sizeof(Site),	1,	F))	return;
	name.resize(in);	char	temp[64];
	for(uint	i=0;	i<in;	i++){
		if(!fread(temp,	64,	1,	F))	return;
		name[i]=temp;
	}
	offset=ftell(F);
}

void	probin::write_bin(FILE	*F,	uint64_t	S,	uint	N){
	vector<uint16_t>	buff(N*in*2);
	fseek(F,	offset+S*in*4,	SEEK_SET);
	if(!fread(&buff[0],	N*in*4,	1,	F))	return;

	char	fname[1024];
	sprintf(fname,	"%s%s_%09u_%09u.bin",	fld.c_str(),	chr.c_str(),	site[S].pos,	site[S+N-1].pos);
	gzFile	file=gzopen(fname,	"wt");
	fprintf(stderr,	"%s\n",	fname);
	gzprintf(file,	"chr\tpos\tallele");
	for(uint	i=0;	i<in;	i++)	gzprintf(file,	"\t%s",	name[i].c_str());
	for(uint	i=0;	i<N;	i++){
		gzprintf(file,	"\n%s\t%u\t%c%c",	site[S+i].chr,	site[S+i].pos,	site[S+i].all[0],	site[S+i].all[1]);  // outputs the probabilities
		uint16_t	*p=&buff[i*in*2];
		for(uint	j=0;	j<in;	j++,	p+=2)	gzprintf(file,	"\t%g %g",	2e-5*p[0],	2e-5*p[1]);
	}
	gzprintf(file,	"\n");
	gzclose(file);
}

void	probin::work(const	char	*F){  // open the prob file
	FILE	*f=fopen(F,	"rb");
	if(f==NULL){	fprintf(stderr,	"fail to open %s\n",	F);	exit(0);	}
	load_head(f);
	
	uint	sta=0xffffffff,	end=0;
	for(uint	i=0;	i<site.size();	i++)	if(chr==site[i].chr){
		if(i>end)	end=i;
		if(i<sta)	sta=i;
	}
	for(uint	i=sta;	i<=end;	i+=str){
		if(i+bin+str-1<=end)	write_bin(f,	i,	bin);
		else{	write_bin(f,	i,	end+1-i);	break;	}
	}
	fclose(f);
}

void	probin::document(void){
	cerr<<"\nprobin";
	cerr<<"\ngenerate chunked .bin files for each chromosome";
	cerr<<"\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage	probin [options] <in.prob chr>";
	cerr<<"\n	-b <INT>	number of SNPs in a bin (1024)";
	cerr<<"\n	-s <INT>	stride of binning (512)";
	cerr<<"\n	-f <STR>	folder (./)";
	cerr<<"\n\n";
	exit(0);	
}

int	main(int	ac,	char	**av){
	probin	pb;	pb.bin=1024;	pb.str=512;	pb.fld="./";
	int	opt;
	while((opt=getopt(ac,	av,	"b:s:f:"))>=0){
		switch(opt){
		case	'b':	pb.bin=atoi(optarg);	break;
		case	's':	pb.str=atoi(optarg);	break;
		case	'f':	pb.fld=optarg;	break;
		default:	probin::document();
		}
	}
	if(ac<optind+2)	probin::document();
	pb.chr=av[optind+1];	if(pb.fld.size()&&pb.fld[pb.fld.size()-1]!='/')	pb.fld+='/';
	pb.work(av[optind]);
	return	0;
}

