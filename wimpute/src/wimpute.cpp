
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

void Wimpute::WriteToLog( const stringstream & tInput )
{

        
    if(m_bLogIsGz){
        gzprintf(m_gzLogFileStream, tInput.str().c_str());
    }
    else{
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
    m_ofsLogFileStream << tInput.str();
    m_ofsLogFileStream.flush();
    }
};


// this part of the code seems to be responsible for:
// A - finding a set of four haps that are close to the current individual
// B - running the HMM and udating the individual I's haplotypes
// A takes much longer than B

/* CHANGES from impute.cpp:
   moved logging to solve from hmm_like()

*/
fast Wimpute::solve(uint I, uint    &N, fast S, bool P) {  // solve(i,	len,	pen,	n>=bn)

    // pick 4 haplotype indices at random not from individual
    uint p[4];
    for (uint j = 0; j < 4; j++) {
        do p[j] = gsl_rng_get(rng) % hn; while (p[j] / 2 == I);
    }

    // get a probability of the model for individual I given p
    fast curr = hmm_like(I, p);

    // pick a random haplotype to replace with another one from all
    // haplotypes.  calculate the new probability of the model given
    // those haplotypes.
    // accept new set if probability has increased.
    // otherwise, accept with penalized probability
    for (uint n = 0; n < N; n++) {  // fixed number of iterations

        uint rp = gsl_rng_get(rng) & 3, oh = p[rp];    
        do p[rp] = gsl_rng_get(rng) % hn; while (p[rp] / 2 == I);
        fast prop = hmm_like(I, p);
        bool bAccepted = false;
        if (prop > curr || gsl_rng_uniform(rng) < expf((prop - curr) * S)) {
            curr = prop;
            bAccepted = true;
        }
        else p[rp] = oh;

        // log results
        stringstream message;
        message << m_nIteration << "\t" << I << "\t" <<  curr << "\t" << bAccepted << "\n";
        WriteToLog( message );

    }

    // if we have passed the burnin cycles (n >= bn)
    // start sampling the haplotypes for output
    if (P) {
        uint16_t *pa = &pare[I * in];
        for (uint i = 0; i < 4; i++) pa[p[i] / 2]++;
    }
    hmm_work(I, p, S);
    return curr;
}

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
