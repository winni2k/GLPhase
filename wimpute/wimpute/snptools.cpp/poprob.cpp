/**
@file	poprob.cpp
@brief	merge and transpose individual likelihood files to a population likelihood file
@author	Yi Wang
@date	04/01/2011
*/
// pool probability

#include	<algorithm>
#include	<iostream>
#include	<stdint.h>
#include	<cstdlib>
#include	<cstring>
#include	<fstream>
#include	<sstream>
#include	<cstdio>
#include	<string>
#include	<vector>
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

class	poprob{
private:
	uint	sn,	fn;
	vector<Site>	site;
	vector<string>	file;	
public:
	uint64_t	bn;
	static	void	document(void);
	bool	load_site(const	char	*F);
	bool	load_file(const	char	*F);
	bool	write_out(const	char	*F);
};

bool	poprob::load_site(const	char	*F){
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
	cerr<<"unique\t"<<sn<<endl;
	return	true;
}

bool	poprob::load_file(const	char	*F){
	ifstream	fi(F);
	if(!fi)	return	false;
	string	temp;
	for(fi>>temp;	!fi.eof();	fi>>temp)	file.push_back(temp);
	fi.close();
	cerr<<"listing "<<(fn=file.size())<<"\tfiles\n";	
	return	true;
}

bool	poprob::write_out(const	char	*F){
	bn=bn*1024*256/fn;	if(bn>sn)	bn=sn;
	vector<vector<uint16_t>	>	buff(fn);
	vector<uint16_t>	outb(fn*2);
	
	FILE	*outf=fopen(F,	"wb");
	if(outf==NULL)	return	false;
	
	fwrite(&sn,	4,	1,	outf);	fwrite(&fn,	4,	1,	outf);  // write out information 
	fwrite(&site[0],	sizeof(Site)*sn,	1,	outf);
	
	char	temp[64];
	for(uint	i=0;	i<fn;	i++){
		memset(temp,	0,	64);
		size_t	s=file[i].find_last_of('/');
		if(s==string::npos)	s=0;	else	s++;
		for(uint	j=s;	j<file[i].size()&&file[i][j]!='.';	j++)	temp[j-s]=file[i][j];
		fwrite(temp,	64,	1,	outf);
	}
	
	for(uint	m=0;	m<sn;	m+=bn){
		uint	len=m+bn<sn?bn:sn-m;
		cerr<<"reading\t"<<len<<'\t'<<site[m].chr<<':'<<site[m].pos<<"<=>"<<site[m+len-1].chr<<':'<<site[m+len-1].pos<<endl;
		for(uint	i=0;	i<fn;	i++){
			buff[i].resize(len*2);
			FILE	*f=fopen(file[i].c_str(),	"rb");
			if(f==NULL){	cerr<<"fail to open "<<file[i]<<endl;	return	false;	}
			fseek(f,	m*4,	SEEK_SET);
			if(!fread(&buff[i][0],	len*4,	1,	f)){	cerr<<"fail to read "<<file[i]<<endl;	return	false;	}
			fclose(f);
		}
		for(uint	j=0;	j<len;	j++){
			uint	q=j*2;
			uint16_t	*p=&outb[0];
			for(uint	i=0;	i<fn;	i++,	p+=2){	p[0]=buff[i][q];	p[1]=buff[i][q+1];	}
			fwrite(&outb[0],	4*fn,	1,	outf);s
		}
	}
	fclose(outf);
	return	true;
}
	
void	poprob::document(void){
	cerr<<"\npoprob";
	cerr<<"\nmerge and transpose individual likelihood files to a population likelihood file";
	cerr<<"\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage	poprob [options] <site.vcf raws.list out.prob>";
	cerr<<"\n	-b <INT>	buffer size in MB (1024)";
	cerr<<"\n\n";
	exit(0);
}

int	main(int	ac,	char	**av){
	poprob	pp;	pp.bn=1024;
	
	int	opt;
	while((opt=getopt(ac,	av,	"b:"))>=0){
		switch(opt){
		case	'b':	pp.bn=atoi(optarg);	break;
		default:	poprob::document();
		}
	}
	if(ac<optind+3)	poprob::document();
	if(!pp.load_site(av[optind])){	cerr<<"fail to open "<<av[optind]<<endl;	return	0;	}
	if(!pp.load_file(av[optind+1])){cerr<<"fail to open "<<av[optind+1]<<endl;	return	0;	}
	if(!pp.write_out(av[optind+2])){cerr<<"fail to open "<<av[optind+2]<<endl;	return	0;	}
	return	0;
}
