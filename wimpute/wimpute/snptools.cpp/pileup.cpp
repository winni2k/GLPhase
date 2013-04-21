#include	<iostream>
#include	<stdint.h>
#include	<bzlib.h>
#include	<cstdlib>
#include	<vector>
#include	<cmath>
#include	<ctime>
#include	<sam.h>
using	namespace	std;
typedef	unsigned	uint;
#define	BlockSize	65536

// Uses library from samtools http://samtools.sourceforge.net/samtools/sam/Functions/Functions.html#//apple_ref/c/func/sampileup

// define pileup class
class	pileup{
private:
	vector<uint16_t>	ebd;
	uint32_t	*length,	tid;
	char	**name;	 // pointer to pointer
	FILE	*file;
	gzFile	index;
	time_t	start;	
	bool	save_chr(void);
public:
	static	float	phred[256];
	static	uint8_t	table[16];
	static	void	make_code(void);
	static	int	pileup_func(uint tid,	uint pos,	int n,	const	bam_pileup1_t	*pl,	void	*data);
	static	void	document(void);
	
	bool	scan_bam(const	char	*F);
};

float	pileup::phred[256];
uint8_t	pileup::table[16];

//not sure what the purpose of the table is... has values 01424443444444444
void	pileup::make_code(void){
	for(uint	i=0;	i<256;	i++)	phred[i]=1-exp(-0.1*i);
	for(uint	i=0;	i<16;	i++)	table[i]=4;
	table[1]=0;	table[2]=1;	table[4]=2;	table[8]=3;
}

// startup document

void	pileup::document(void){
	cerr<<"\npileup";
	cerr<<"\nextract read depth from BAM";
	cerr<<"\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
	cerr<<"\nusage	pileup [1.bam 2.bam ...]";
	cerr<<"\n\n";
	exit(0);	
}

int	pileup::pileup_func(uint tid,	uint pos,	int n,	const	bam_pileup1_t	*pl,	void	*data){  // this pileup function is what computes the ebd
	pileup	*pu=(pileup*)data;
	if(tid!=pu->tid){
		pu->save_chr();	pu->tid=tid;	pu->ebd.assign(pu->length[tid]*6,	0);
		cerr<<'\t'<<pu->name[tid]<<"\tpileuping";
		pu->start=time(NULL);
	}
	float	cnt[6]={0,0,0,0,0,0};
	const	bam_pileup1_t	*p=pl;
	for(int	i=0;	i<n;	i++,	p++)	if(p->b->core.qual){
		float	wmap=phred[p->b->core.qual];  // what is phred score
		if(p->indel>0)	cnt[4]+=wmap;
		else	if(p->indel<0)	cnt[5]+=wmap;
		else	if(!p->is_del){
			uint8_t	c=table[bam1_seqi(bam1_seq(p->b), p->qpos)];  //pointer to alignment... get the base
			if(c<4)	cnt[c]+=wmap*phred[bam1_qual(p->b)[p->qpos]];  // get the quality string
		}
	}
	uint16_t	*pebd=&pu->ebd[pos*6];
	for(uint	i=0;	i<6;	i++)	pebd[i]=(uint16_t)fminf(cnt[i]*10+0.5,	65535.0f);     // find the minimum... 65535.0f is some sort of scaling
	return	0;
}

bool	pileup::save_chr(void){
	if(tid==0xffffffff)	return	true;
	cerr<<"\tcompressing";
	uint32_t	total=length[tid];	
	vector<char>	buff((uint)(BlockSize*12*1.01+1024));
	for(uint	i=0;	i<total;	i+=BlockSize){
		uint	len=i+BlockSize<total?BlockSize:total-i,	bufn=buff.size();
		BZ2_bzBuffToBuffCompress(&buff[0],	&bufn,	(char*)&ebd[i*6],	len*12,	9,	0,	0);
		gzprintf(index,	"%s\t%lu\t%u\n",	name[tid],	ftell(file),	bufn);
		if(!fwrite(&buff[0],	bufn,	1,	file))	return	false;
	}
	cerr<<"\t"<<time(NULL)-start<<"s\n";
	return	true;
}

bool	pileup::scan_bam(const	char	*F){
	tid=0xffffffff;  // unsigned integer 
	samfile_t *bam;
	bam=samopen(F, "rb", 0); // read bam
	if(bam==0)	return	false; // in case a filename is entered with directory structure
	const	char	*p=strrchr(F,	'/');	if(p==NULL)	p=F;	else	p++;
	char	temp[256];

// create the ebd and ebi files
	sprintf(temp,	"%s.ebd",	p);	if((file=fopen(temp,	"wb"))==NULL)	return	false;
	sprintf(temp,	"%s.ebi",	p);	if((index=gzopen(temp,	"wt"))==NULL)	return	false;
	cerr<<endl<<F<<endl;
	name=bam->header->target_name;	length=bam->header->target_len;
	sampileup(bam,	-1,	pileup_func,	this);  // run pileup function
	save_chr();
	gzclose(index);	fclose(file);	samclose(bam);
	return	true;
}


int	main(int	ac,	char	**av){
	if(ac<2)	pileup::document();
	pileup::make_code();
		
	pileup	pu;
	for(int	i=1;	i<ac;	i++)	pu.scan_bam(av[i]);
	return	0;
}

