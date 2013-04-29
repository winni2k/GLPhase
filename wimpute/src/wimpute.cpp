
#include "wimpute.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

// return the probability of the model given the input haplotypes P and
// emission and transition matrices of individual I
// call Impute::hmm_like and print out result
fast    Wimpute::hmm_like(uint I, uint *P) {

//    cerr << "calling Wimpute::hmm_like()\n";

    // call Impute::hmm_like
    fast curr = Impute::hmm_like( I, P );

    // print result: iteration individual likelihood
    stringstream message;
    message << m_nIteration << "\t" << I << "\t" <<  curr << "\n";
    WriteToLog( message.str() );
    
    return curr;
}

void Wimpute::SetLog( const string & sLogFile )
{
    m_sLogFile = sLogFile;
    m_bLogIsGz = 0;

    size_t found = sLogFile.find(".gz", sLogFile.length() - 3);
    if( found!=std::string::npos ){
        m_bLogIsGz = 1;
        m_gzLogFileStream = gzopen(m_sLogFile.c_str(), "wt");
    }
    else{

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
    }
};

void    Wimpute::estimate(void) {
    cerr.setf(ios::fixed);
    cerr.precision(3);
    cerr << "iter\tpress\tlike\tfold\n";
    
    // n is number of cycles = burnin + sampling cycles
    // increase penalty from 2/bn to 1 as we go through burnin
    // iterations.    
    for (uint n = 0; n < bn + sn; n++) {
        m_nIteration = n;
        fast sum = 0, pen = fminf(2 * (n + 1.0f) / bn, 1), iter = 0;
        pen *= pen;  // pen = 1 after bn/2 iterations
        for (uint i = 0; i < in; i++) {
            uint len = nn * in;  // nn is number of folds, in = num individuals
            sum += solve(i, len, pen, n >= bn);  // call solve=> inputs the sample number,
            iter += len;
        }
        swap(hnew, haps);
        if (n >= bn) for (uint i = 0; i < in; i++) replace(i);  // call replace
        cerr << n << '\t' << pen << '\t' << sum / in / mn << '\t' << iter / in / in << '\r';
    }
    cerr << endl;
    result();    // call result
}


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
