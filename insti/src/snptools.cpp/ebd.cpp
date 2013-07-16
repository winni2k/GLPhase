#include	<stdint.h>
#include	<bzlib.h>
#include	<cstdlib>
#include	<cstring>
#include	<string>
#include	<vector>
#include	<zlib.h>
using	namespace	std;

struct	EBI{
	char	chr[20];
	uint64_t	off;
	uint32_t	len;
};

struct	EBD{
	vector<uint16_t>	cnt;
	vector<EBI>	ebi;
	string	file;
	
	bool	bind(const	char	*C);
	bool	load(uint32_t	B);
};

bool	EBD::bind(const	char	*C){
	string	s=file;	*s.rbegin()='i';
	gzFile	f=gzopen(s.c_str(),	"rt");
	if(f==Z_NULL)	return	false;
	char	temp[1024];
	EBI	idx;
	while(gzgets(f,	temp,	1024)!=Z_NULL){
		sscanf(temp,	"%s%lu%u",	idx.chr,	&idx.off,	&idx.len);
		if(C==NULL||!strcmp(idx.chr,	C))	ebi.push_back(idx);
	}
	gzclose(f);
	return	true;
}

bool	EBD::load(uint32_t	B){
	if(B>=ebi.size())	return	false;
	cnt.assign(65536*6,	0);
	vector<char>	temp(ebi[B].len);
	FILE	*f=fopen(file.c_str(),	"rb");
	if(f==NULL)	return	false;
	fseek(f,	ebi[B].off,	SEEK_SET);
	if(!fread(&temp[0],	ebi[B].len,	1,	f))	return	false;
	fclose(f);
	uint32_t	len=65536*12;
	BZ2_bzBuffToBuffDecompress((char*)&cnt[0],	&len,	&temp[0],	ebi[B].len,	0,	0);
	return	true;
}

