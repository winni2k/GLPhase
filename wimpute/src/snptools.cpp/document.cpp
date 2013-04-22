#include	<algorithm>
#include	<iostream>
#include	<dirent.h>
#include	<cstdlib>
#include	<cstring>
#include	<fstream>
#include	<vector>
#include	<string>
using	namespace	std;

void	print_one(string&	F){
	char	buffer[1024];
	ifstream	fi(F.c_str());
	if(!fi)	return;
	for(fi.getline(buffer,	1024);	!fi.eof();	fi.getline(buffer,	1024)){
		char	*p=strstr(buffer,	"cerr<<\"\\n");
		if(p==NULL)	continue;
		for(p+=9;	p<buffer+1024;	p++){
			if(*p=='"'){	printf("\n");	break;	}
			if(*p=='\\'){	printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n");	break;	}
			putchar(*p);
		}
	}
	fi.close();
}	

int	main(void){
	vector<string>	name;
	
	DIR	*dir=opendir(".");
	if(dir==NULL){	fprintf(stderr,	"can't open current directory\n");	return	0;	}
	struct	dirent  *ptr;
	while((ptr=readdir(dir))!=NULL){
		string	s=ptr->d_name;
		if(s.find(".cpp")!=string::npos)	name.push_back(s);
	}
	closedir(dir);
	sort(name.begin(),	name.end());
	
	cout<<"SNPTools	1.0\n";
	cout<<"pileline for NGS SNP analysis\n";
	cout<<"author	Yi Wang @ Fuli Yu' Group\n";
	cout<<"Baylor College of Medicine Human Genome Sequencing Center.\n";
	cout<<"All rights reserved.\n\n";
	cout<<"Workflow summary:\n";
	cout<<"pileup->varisite->bamodel->poprob->probin->impute->hapfuse\n";
	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
		
	for(uint	i=0;	i<name.size();	i++)	print_one(name[i]);
		
	cout<<"Installation\n";
	cout<<"The software is based on samtools's API library: libbam.a\n";
	cout<<"The software is based on tabix's API library: libtabix.a\n";
	cout<<"The software is based on GNU Scientific library\n";		
	cout<<"User should make sure that they are installed before make SNPTools\n";
	cout<<"Edit makefile if they are not installed in default directories\n";
	cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	cout<<"License\n";
	cout<<"Permission is hereby granted, free of charge, to any person obtaining a copy\n";
	cout<<"of this software and associated documentation files (the 'Software'), to deal\n";
	cout<<"in the Software without restriction, including without limitation the rights\n";
	cout<<"to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n";
	cout<<"copies of the Software, and to permit persons to whom the Software is\n";
	cout<<"furnished to do so, subject to the following conditions:\n\n";
	cout<<"The above copyright notice and this permission notice shall be included in\n";
	cout<<"all copies or substantial portions of the Software.\n\n";
	cout<<"THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n";
	cout<<"IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n";
	cout<<"FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n";
	cout<<"AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n";
	cout<<"LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n";
	cout<<"OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n";
	cout<<"THE SOFTWARE.\n";
	return	0;
}
