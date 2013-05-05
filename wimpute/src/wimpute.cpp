
#include "wimpute.h"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;


//initializyng static member variables
int Wimpute::s_iEstimator;
uint Wimpute::s_uParallelChains;
uint Wimpute::s_uCycles;

// return the probability of the model given the input haplotypes P and
// emission and transition matrices of individual I
// call Impute::hmm_like and print out result
fast    Wimpute::hmm_like(uint I, uint *P) {

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


bool    Wimpute::load_bin(const char *F) {

    bool retval = Impute::load_bin(F);
    if(retval == false){ return false; }
    
    // setting number of cycles to use
    // here is best place to do it because in is defined in load_bin()

    if( s_uCycles > 0){
        m_uCycles = s_uCycles;
    }
    else{
        m_uCycles = nn * in;  // this was how snptools does it
    }

    return true;
}

// this part of the code seems to be responsible for:
// A - finding a set of four haps that are close to the current individual
// B - running the HMM and udating the individual I's haplotypes
// A takes much longer than B

/* CHANGES from impute.cpp:
   moved logging to solve from hmm_like()
   cycles is now stored in a private member variable and defined after load_bin() time
*/

// solve(individual, number of cycles, penalty, burnin?)
fast Wimpute::solve(uint I, uint    &N, fast S, bool P) {

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
        message << m_nIteration << "\t" << I << "\t" <<  prop << "\t" << bAccepted << "\n";
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

/* CHANGES from impute.cpp
   added a member variable for n so
*/

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
            sum += solve(i, m_uCycles, pen, n >= bn);  // call solve=> inputs the sample number,
            iter += m_uCycles;
        }
        swap(hnew, haps);
        if (n >= bn) for (uint i = 0; i < in; i++) replace(i);  // call replace
        cerr << n << '\t' << pen << '\t' << sum / in / mn << '\t' << iter / in / in << '\r';
    }
    cerr << endl;
    result();    // call result
}

/* estimate_EMC -- Evolutionary Monte Carlo
   Here we try to increase the speed of convergence by
   running a parallel chain evolutionary monte carlo scheme.

   The reference for this implementation is
   "Advanced Markov Choin Monte Carlo Methods" by Liang, Liu and Carroll
   first edition?, 2010, pp. 128-132
 */

// solve(individual, number of cycles, penalty, burnin?)
fast Wimpute::solve_EMC(const uint I, const uint    &N, const fast S, const bool P) {

    cerr << "Entering solve_EMC...\n";
    // for lack of a better place, define free parameters here
    fast fMutationRate = 0.3f; // see p.134 of Liang et al.
    fast fSelectTemp = 1 / S;
    
    // initialize emc chains with increasing temperatures
    vector <EMCChain> vcChains;
    uint uNumChains = 5;
    vector <uint> vuChainTempHierarchy; // index of Chains sorted by temperature, ascending
    for (uint i = 1; i<=uNumChains; i++){
        EMCChain cChain( i / S, fSelectTemp, I, in);

        // initialize current likelihood
        cChain.setLike( hmm_like(cChain.m_uI, cChain.m_auParents) );
        vcChains.push_back(cChain);
        vuChainTempHierarchy.push_back(i);
    }


    // pick a random haplotype to replace with another one from all
    // haplotypes.  calculate the new probability of the model given
    // those haplotypes.
    // accept new set if probability has increased.
    // otherwise, accept with penalized probability
    for (uint n = 0; n < N; n++) {  // fixed number of iterations

        cerr << "\tIteration " << n << "\n";
        //  now choose whether to mutate or crossover
        bool bMutate =  gsl_rng_uniform(rng) > fMutationRate;
        bool bMutCrossAccepted = false;
        fast prop;
        if(bMutate){
            // mutate

            cerr << "\tMutating...\n";
            // choose chain randomly (uniform)
            EMCChain &rcChain = vcChains[gsl_rng_get(rng) % uNumChains];
            
            fast curr = rcChain.getLike();

            // choose parent hap (rp) to mutate
            // replaced hap is stored in oh (Original Hap)
            uint rp = gsl_rng_get(rng) & 3, oh = rcChain.m_auParents[rp];

            // mutate parent hap
            do rcChain.m_auParents[rp] = gsl_rng_get(rng) % hn; while (rcChain.m_auParents[rp] / 2 == I);

            // calculate acceptance probability
            prop = hmm_like(rcChain.m_uI, rcChain.m_auParents);
            if (prop > curr || gsl_rng_uniform(rng) < expf((prop - curr) / rcChain.m_fTemp)) {
                rcChain.setLike(prop);
                bMutCrossAccepted = true;
            }
            else rcChain.m_auParents[rp] = oh;
            
        }
        else{
            //crossover

            cerr << "\tCrossing Over...\n";
            // 1. choose random chain to work on by roulette wheel selection
            fast fTotalProb = 0; // always positive
            for( auto icChain: vcChains){
                fTotalProb += icChain.getCrossProb();
            }
            fast fStopPoint = gsl_rng_uniform(rng) * fTotalProb;
            auto icFirstChain = vcChains.begin();
            while(icFirstChain != vcChains.end()){
                fStopPoint += icFirstChain->getCrossProb();
                if(fStopPoint < 0) break;
                else ++icFirstChain;           
            }

            // 2. select second chain at random (uniform) from remaining chains
            uint uSecondChain;
            do uSecondChain = gsl_rng_get(rng) % uNumChains; while (vcChains[uSecondChain].m_uChainID != icFirstChain->m_uChainID);
            EMCChain &rcSecondChain = vcChains[uSecondChain];

            // only crossover with certain probability depending
            // on crossover probabilities of parents           
            if( gsl_rng_uniform(rng) * fTotalProb * (uNumChains -1)
                > icFirstChain->getCrossProb() + rcSecondChain.getCrossProb() )
                break;

            bMutCrossAccepted = true;
            
            // uniform crossover: find haps to keep
            word cSelection = static_cast<word>( gsl_rng_get(rng) & 15);
            for(uint i = 0; i < 4; i++){
                
                // if bit at location i is 1, exchange haps
                if( test( &cSelection, i) ) {
                    uint oh = icFirstChain->m_auParents[i];
                    icFirstChain->m_auParents[i] = rcSecondChain.m_auParents[i];
                    rcSecondChain.m_auParents[i] = oh;
                }
            }
        }

        // now try uNumChains exchanges
        cerr << "\tExchanging...\n";
        uint uNumExchanges = 0;
        for( uint i = 0; i < uNumChains; i++){
            uint uFirstChainIndex = gsl_rng_get(rng) % uNumChains;
            auto rcFirstChain = vcChains[ vuChainTempHierarchy[ uFirstChainIndex ] ];

            // selecting second chain
            uint uSecondChainIndex;
            if (uFirstChainIndex == 0)
                uSecondChainIndex = uFirstChainIndex + 1;
            else if ( uFirstChainIndex == uNumChains - 1)
                uSecondChainIndex = uFirstChainIndex - 1;
            else if( gsl_rng_get(rng) & 1 )
                uSecondChainIndex = i - 1;
            else
                uSecondChainIndex = i + 1;
            
            auto rcSecondChain = vcChains[ vuChainTempHierarchy[ uSecondChainIndex ] ];

            // MH step for exchange
            fast fAcceptProb = fminf( expf( (rcFirstChain.getLike() - rcSecondChain.getLike())
                                            * ( (1/rcFirstChain.m_fTemp) - (1/rcSecondChain.m_fTemp)))
                                      , 1);

            // exchange with acceptance probability
            if( gsl_rng_uniform(rng) < fAcceptProb){

                // exchange temperatures
                fast fTemp = rcFirstChain.m_fTemp;
                rcFirstChain.setTemp(rcSecondChain.m_fTemp);
                rcSecondChain.setTemp(fTemp);

                // exchange location in vcChains
                std::swap(vuChainTempHierarchy[uFirstChainIndex], vuChainTempHierarchy[uSecondChainIndex]);
                ++ uNumExchanges;
            }
                
        }
        
        // log results
        stringstream message;
        message << m_nIteration << "\t" << I << "\t" <<  prop << "\t"
                << bMutCrossAccepted << "\t" << bMutate << "\t" << uNumExchanges << "\n";
        WriteToLog( message );

    }

    // now select a chain for sampling according to roulette wheel selection
    fast fTotalProb = 0; // always positive
    for( auto icChain: vcChains){
        fTotalProb += icChain.getCrossProb();
    }
    fast fStopPoint = gsl_rng_uniform(rng) * fTotalProb;
    auto icFirstChain = vcChains.begin();
    while(icFirstChain != vcChains.end()){
        fStopPoint += icFirstChain->getCrossProb();
        if(fStopPoint < 0) break;
        else ++icFirstChain;           
    }
    
    // if we have passed the burnin cycles (n >= bn)
    // start sampling the haplotypes for output
    if (P) {
        uint16_t *pa = &pare[I * in];
        for (uint i = 0; i < 4; i++) pa[icFirstChain->m_auParents[i] / 2]++;
    }

    cerr << "Updating individual " << I << "\n";
    // update haplotypes of I
    hmm_work(I, icFirstChain->m_auParents, S);
    return icFirstChain->getLike();
}

/* estimate_EMC -- Evolutionary Monte Carlo
   Here we try to increase the speed of convergence by
   running a parallel chain evolutionary monte carlo scheme.

   The reference for this implementation is
   "Advanced Markov Choin Monte Carlo Methods" by Liang, Liu and Carroll
   first edition?, 2010, pp. 128-132
 */

void    Wimpute::estimate_EMC(void) {
    cerr.setf(ios::fixed);
    cerr.precision(3);
    cerr << "Running Evolutionary Monte Carlo\n";
    cerr << "iter\tpress\tlike\tfold\n";
   
    // n is number of cycles = burnin + sampling cycles
    // increase penalty from 2/bn to 1 as we go through burnin
    // iterations.    
    for (uint n = 0; n < bn + sn; n++) {
        m_nIteration = n;
        fast sum = 0, pen = fminf(2 * (n + 1.0f) / bn, 1), iter = 0;
        pen *= pen;  // pen = 1 after bn/2 iterations
        for (uint i = 0; i < in; i++) {           
            sum += solve_EMC(i, m_uCycles, pen, n >= bn);  // call solve=> inputs the sample number,
            iter += m_uCycles;
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
    cerr << "\n	-E <integer>	choice of estimation algorithm (0)";
    cerr << "\n			0 - Metropolis Hastings with simulated annealing";
    cerr << "\n			1 - Evolutionary Monte Carlo with -p parallel chains";
    cerr << "\n\t-p <integer>	number of parallel chains to use in parallel estimation algorithms";
    cerr << "\n\t		(at least 2, default 5)";
    cerr << "\n\t-C <integer>	number of cycles to estimate an individual's parents before updating";
    cerr << "\n\n";
    exit(1);
}
