#include	<iostream>
#include	<tabix.h>
#include	<bgzf.h>
#include	<cstdio>
#include	<vector>
#include	<cmath>
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

class	prob2vcf{
private:
	uint32_t	mn,	in;
	uint64_t	offset;
	vector<Site>	site;
	vector<string>	name;	
		
	void	load_head(FILE	*F);
	void	vcf_head(BGZF	*F,	char	*B);
public:
	static	void	document(void);
	void	work(const	char	*I,	const	char	*O,	const	char	*C);
};

void	prob2vcf::document(void){
	cerr<<"\nprob2vcf";
	cerr<<"\nconvert binary .prob file to tabixed .vcf.gz file";
	cerr<<"\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage	prob2vcf <in.prob out.vcf.gz chr>";
	cerr<<"\n\n";
	exit(0);
}

void	prob2vcf::load_head(FILE	*F){
	if(!fread(&mn,	4,	1,	F))	return;
	if(!fread(&in,	4,	1,	F))	return;
	site.resize(mn);
	if(!fread(&site[0],	mn*sizeof(Site),	1,	F))	return;
	name.resize(in);	char	temp[64];
	for(uint	i=0;	i<in;	i++){
		if(!fread(temp,	64,	1,	F))	return;
		name[i]=temp;
	}
	offset=ftell(F);
	cerr<<"sites\t"<<mn<<endl;	
	cerr<<"sample\t"<<in<<endl;
}

void	prob2vcf::vcf_head(BGZF	*F,	char	*B){
	char	*p;
	p=B;	p+=sprintf(p,	"##fileformat=VCFv4.0\n");	bgzf_write(F,	B,	p-B);
	p=B;	p+=sprintf(p,	"##source=BCM:SNPTools:prob2vcf\n");	bgzf_write(F,	B,	p-B);
	p=B;	p+=sprintf(p,	"##reference=GRCh37\n");	bgzf_write(F,	B,	p-B);
	p=B;	p+=sprintf(p,	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"genotype\">\n");	bgzf_write(F,	B,	p-B);
	p=B;	p+=sprintf(p,	"##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"three log10-scaled likelihoods for RR,RA,AA genotypes\">\n");	bgzf_write(F,	B,	p-B);
	p=B;	p+=sprintf(p,	"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
	for(uint	i=0;	i<in;	i++)	p+=sprintf(p,	"\t%s",	name[i].c_str());	
	p+=sprintf(p,	"\n");	bgzf_write(F,	B,	p-B);
}

void	prob2vcf::work(const	char	*I,	const	char	*O,	const	char	*C){
	FILE	*fi=fopen(I,	"rb");
	if(fi==NULL){	fprintf(stderr,	"fail to open %s\n",	I);	exit(0);	}
	load_head(fi);
	
	BGZF	*bgzf=bgzf_open(O,	"w");
	char	*line=new	char[1024*1024];
	vcf_head(bgzf,	line);
	vector<uint16_t>	buff(in*2);
	
	uint64_t	beg=mn,	end=0;
	for(uint	m=0;	m<mn;	m++)	if(!strcmp(site[m].chr,	C)){
		if(m<beg)	beg=m;	if(m>end)	end=m;
	}
	fseek(fi,	offset+beg*in*4,	SEEK_SET);
	for(uint	m=beg;	m<=end;	m++){
		if(!fread(&buff[0],	4*in,	1,	fi)){	fprintf(stderr,	"fail to read %s\n",	I);	exit(0);	}
		char	*p=line;	
		p+=sprintf(p,	"%s\t%u\t.\t%c\t%c\t100\tPASS\t.\tGT:GL",	site[m].chr,	site[m].pos,	site[m].all[0],	site[m].all[1]);
		for(uint	i=0;	i<in;	i++){
			float	ph=2e-5f*buff[i*2],	pb=2e-5f*buff[i*2+1],	pa=1.0f-ph-pb;
			p+=sprintf(p,	"\t./.:%g,%g,%g",	log10f(fmaxf(pa,	1e-5f)),	log10f(fmaxf(ph,	1e-5f)),	log10f(fmaxf(pb,	1e-5f)));
		}
		p+=sprintf(p,	"\n");
		bgzf_write(bgzf,	line,	p-line);
		if(!(m%100000))	cerr<<1e-6*m<<"M\r";
	}
	delete	[]	line;
	bgzf_close(bgzf);
	fclose(fi);
}

int	main(int	ac,	char	**av){
	if(ac!=4)	prob2vcf::document();
	prob2vcf	pv;
	pv.work(av[1],	av[2],	av[3]);
	return	0;
}
