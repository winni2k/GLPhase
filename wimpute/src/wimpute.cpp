
#include "wimpute.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;
void Wimpute::SetLog( const string & sLogFile )
{
    m_sLogFile = sLogFile;
    
    // open logFileStream for append if it's not yet open
    if( !m_ofsLogFileStream ){
        m_ofsLogFileStream.open(m_sLogFile);
    }
    //otherwise close and open again
    else
    {
        m_ofsLogFileStream.close();
        m_ofsLogFileStream.open(m_sLogFile);
    }
        
    // exit if the file cannot be opened
    if( !m_ofsLogFileStream.is_open() ){
        cerr << "could not open log file "<< m_sLogFile <<" for writing" << endl;
        exit(1);
    }
};

void Wimpute::WriteToLog( const string & sInput )
{
    // open logFileStream for append if it's not yet open
    if( !m_ofsLogFileStream ){
        m_ofsLogFileStream.open(m_sLogFile);
    }

    // exit if the file cannot be opened
    if( !m_ofsLogFileStream.is_open() ){
        cerr << "could not open log file "<< m_sLogFile <<" for writing" << endl;
        exit(1);
    }

    // write input to log file
    m_ofsLogFileStream << sInput;
    m_ofsLogFileStream.flush();
};

void    Wimpute::document(void) {
    cerr << "\nimpute";
    cerr << "\nhaplotype imputation by cFDSL distribution";
    cerr << "\nauthor	Yi Wang @ Fuli Yu' Group @ BCM-HGSC";
    cerr << "\nusage	impute [options] 1.bin 2.bin ...";
    cerr << "\n	-d <density>	relative SNP density to Sanger sequencing (1)";
    cerr << "\n	-b <burn>	burn-in generations (56)";
    cerr << "\n	-l <file>	list of input files";
    cerr << "\n	-m <mcmc>	sampling generations (200)";
    cerr << "\n	-n <fold>	sample size*fold of nested MH sampler iteration (2)";
    cerr << "\n	-t <thread>	number of threads (0=MAX)";
    cerr << "\n	-v <vcf>	integrate known genotype in VCF format";
    cerr << "\n	-c <conf>	confidence of known genotype (0.9998)";
    cerr << "\n	-x <gender>	impute x chromosome data";
    cerr << "\n	-e <file>	write log to file";
    cerr << "\n\n";
    exit(0);
}
