#include "insti.hpp"

//require c++11
static_assert(__cplusplus > 199711L, "Program requires C++11 capable compiler");

using namespace std;

/*
  #define DEBUG
  #define DEBUG2
  #define DEBUG3
*/

#ifdef DEBUG
#define DEBUG_MSG(str) do { std::cerr << str; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#ifdef DEBUG2
#define DEBUG_MSG2(str) do { std::cerr << str; } while( false )
#else
#define DEBUG_MSG2(str) do { } while ( false )
#endif

#ifdef DEBUG3
#define DEBUG_MSG3(str) do { std::cerr << str; } while( false )
#else
#define DEBUG_MSG3(str) do { } while ( false )
#endif

//initializyng static member variables

// Metropolis Hastings with Annealing is default
int Insti::s_iEstimator = 0;

// number of parallel chains to use for parallel estimators
unsigned Insti::s_uParallelChains = 5;
unsigned Insti::s_uCycles = 0;
bool Insti::s_bIsLogging = false;
bool Insti::s_bKickStartFromRef = false;
string Insti::s_sRefLegendFile = "";
string Insti::s_sRefHapFile = "";
string Insti::s_scaffoldHapsFile = "";
string Insti::s_scaffoldSampleFile = "";
double Insti::s_scaffoldFreqCutoff = 0.05;
unsigned Insti::s_uNumClusters = 0;
unsigned Insti::s_uClusterType = 0;
unsigned Insti::s_uSABurninGen = 28;
unsigned Insti::s_uNonSABurninGen = 28;

// start clustering after simulated annealing burnin
unsigned Insti::s_uStartClusterGen = Insti::s_uNonSABurninGen;

// return the probability of the model given the input haplotypes P and
// emission and transition matrices of individual I
// call Impute::hmm_like and print out result
fast    Insti::hmm_like(unsigned I, uint* P) {
    return Impute::hmm_like(I, P);
}

void Insti::SetLog(const string & sLogFile) {

    s_bIsLogging = true;
    m_sLogFile = sLogFile;
    m_bLogIsGz = false;
    size_t found = sLogFile.find(".gz", sLogFile.length() - 3);

    if (found != std::string::npos) {
        m_bLogIsGz = 1;
        m_gzLogFileStream = gzopen(m_sLogFile.c_str(), "wt");
    } else {

        // open logFileStream for append if it's not yet open
        if (!m_ofsLogFileStream)
            m_ofsLogFileStream.open(m_sLogFile);

        //otherwise close and open again
        else {
            m_ofsLogFileStream.close();
            m_ofsLogFileStream.open(m_sLogFile);
        }

        // exit if the file cannot be opened
        if (!m_ofsLogFileStream.is_open()) {
            cerr << "could not open log file " << m_sLogFile << " for writing" << endl;
            exit(1);
        }
    }

    cerr << "Logging to:\t" << m_sLogFile << endl;
};

void Insti::WriteToLog(const string & tInput) {
    if (! s_bIsLogging)
        return;

    if (m_bLogIsGz)
        gzprintf(m_gzLogFileStream, tInput.c_str());
    else {

        // open logFileStream for append if it's not yet open
        if (!m_ofsLogFileStream)
            m_ofsLogFileStream.open(m_sLogFile);

        // exit if the file cannot be opened
        if (!m_ofsLogFileStream.is_open()) {
            cerr << "could not open log file " << m_sLogFile << " for writing" << endl;
            exit(1);
        }

        // write input to log file
        m_ofsLogFileStream << tInput;
        m_ofsLogFileStream.flush();
    }

    DEBUG_MSG("wrote something" << endl);
};

// what to log when given an EMCChain
void Insti::WriteToLog(const EMCChain& rcChain, const bool bMutate) {
    stringstream message;
    message << m_nIteration << "\t" << rcChain.m_uI << "\t" <<  rcChain.getLike() <<
            "\t"
            << rcChain.m_uChainID << "\t" << rcChain.m_fTemp << "\t"
            << bMutate << endl;
    WriteToLog(message.str());
}

// Roulette Wheel Selection, returns index of chain selected
int Insti::RWSelection(const vector <EMCChain> &rvcChains) {
    double dTotalProb = 0; // always positive, could  be larger than 1

    for (const auto& icChain : rvcChains) {
        dTotalProb += icChain.getSelection();
        DEBUG_MSG2("\ttotal prob:\t" << dTotalProb << endl);
    }

    assert(dTotalProb > 0);
    assert(dTotalProb < std::numeric_limits<double>::max());
    double dStopPoint = gsl_rng_uniform(rng) * dTotalProb;
    assert(dStopPoint > 0);
    DEBUG_MSG2("\t..." << endl);
    DEBUG_MSG2("\ttotalProb:\t" << dTotalProb << endl);
    int iChainIndex = 0;

    do
        dStopPoint -= rvcChains[iChainIndex].getSelection();

    while (dStopPoint > 0 && ++iChainIndex);

    assert(iChainIndex >= 0);
    assert(iChainIndex < static_cast<int>(rvcChains.size()));
    return iChainIndex;
}


bool    Insti::load_bin(const char *F) {
    bool bRetVal = Impute::load_bin(F);

    if (bRetVal  == false) return false;

    // setting number of cycles to use
    // here is best place to do it because in is defined in load_bin()
    if (s_uCycles > 0)
        m_uCycles = s_uCycles;
    else {
        m_uCycles = nn * in;  // this was how snptools does it
    }

    return true;
}

// HAPS/SAMPLE sample file
void Insti::OpenSample(string sampleFile, vector<string>& IDs) {

    // read in sample file
    ifile sampleFD(sampleFile);
    string buffer;

    // discard header
    unsigned uLineNum = 0;

    while (getline(sampleFD, buffer, '\n')) {
        uLineNum++;
        vector<string> tokens;
        sutils::tokenize(buffer, tokens);

        // make sure header start is correct
        if (uLineNum == 1) {
            vector<string> headerTokenized;
            string header = "ID_1 ID_2 missing";
            sutils::tokenize(header, headerTokenized);

            for (unsigned i = 0; i != headerTokenized.size(); i++) {
                try {
                    if (tokens[i] != headerTokenized[i])
                        throw myException("Error in sample file (" + sampleFile +
                                          "): header start does not match:\n\t" + header +
                                          "\n.  Instead the first line of the header is:\n\t" + buffer + "\n");
                } catch (exception& e) {
                    cout << e.what() << endl;
                    exit(2);
                }
            }

            continue;
        }

        // now add sample id to vector
        IDs.push_back(tokens[0]);
    }

    try {
        if (name.size() != IDs.size())
            throw myException("bin and scaffold IDs list lengths (" + sutils::uint2str(
                                  name.size()) + " and " + sutils::uint2str(IDs.size()) +
                              " respectively) do not agree.");
    } catch (exception& e) {
        cout << e.what() << endl;
        exit(2);
    }

    // make sure ids match names pulled from bin files
    for (unsigned i = 0; i < name.size(); i++) {
        try {
            if (name[i] != IDs[i])
                throw myException("ID number " + sutils::uint2str(i) +
                                  " (0-based): bin and IDs from sample file do not match (" + name[i] + " and " +
                                  IDs[i] + " respectively) do not agree.");
        } catch (exception& e) {
            cout << e.what() << endl;
            exit(2);
        }
    }
}

//read in the haps file
unsigned Insti::OpenHaps(string hapFile, vector<uint64_t> & vHaps) {

    ifile hapsFD(hapFile);
    string buffer;
    unsigned uLineNum = 0;
    unsigned uNumHaps = 0;

    while (getline(hapsFD, buffer, '\n')) {
        uLineNum++;
        vector<string> tokens;
        sutils::tokenize(buffer, tokens);

        // make sure header start is correct
        if (uLineNum == 1) {

            // count number of haps
            uNumHaps = tokens.size();

            // resize the vector for holding haps
            vHaps.resize(uNumHaps * wn);
        }

        try {
            if (tokens.size() != uNumHaps)
                throw myException("Every row of haplotypes file must have the same number of columns");
        } catch (exception& e) {
            cout << e.what() << endl;
            exit(1);
        }

        // store haplotypes
        for (unsigned i = 0; i != tokens.size(); i++) {
            int val = atoi(tokens[i].c_str());

            if (val == 0)
                set0(&vHaps[i * wn], uLineNum - 1);
            else if (val == 1)
                set1(&vHaps[i * wn], uLineNum - 1);
            else {
                cerr << "All alleles are not 0 or 1 In haplotypes file " << hapFile <<
                     " at line number " << uLineNum << endl;
                throw 1;
            }
        }
    }

    if (uNumHaps == 0) {
        cerr << "num haps is 0.  Haps file empty?";
        throw 2;
    }

    return uNumHaps;

}


bool Insti::LoadHapsSamp(string sampleFile, string hapsFile,
                         PanelType panelType) {

    // make sure both files are defined
    if (sampleFile.size() == 0) {
        cerr << "Need to define a sample file\n";
        document();
    }

    if (hapsFile.size() == 0) {
        cerr << "Need to define a haps file\n";
        document();
    }

    // load haps file
    CheckPanelPrereqs(panelType);

    // get sample information if we are using a scaffold
    if (panelType == PanelType::SCAFFOLD)
        OpenSample(sampleFile, m_vsScaffoldIDs);

    vector<uint64_t> vHaps;
    unsigned numHaps = OpenHaps(hapsFile, vHaps);
    LoadHaps(vHaps, numHaps, panelType);

    return true;
}

void Insti::LoadHaps(vector<uint64_t> & vHaps, unsigned numHaps,
                     PanelType panelType) {

    // store the haplotypes in the correct place based on what type of panel we are loading
    switch (panelType) {
    case PanelType::REFERENCE:
        vHaps.swap(m_vRefHaps);
        m_uNumRefHaps = numHaps;
        cerr << "Reference panel haplotypes\t" << m_uNumRefHaps << endl;
        break;

    case PanelType::SCAFFOLD:
        vHaps.swap(m_vScaffoldHaps);
        m_uNumScaffoldHaps = numHaps;
        cerr << "Scaffold haplotypes\t" << m_uNumScaffoldHaps << endl;

        try {
            if (m_uNumScaffoldHaps != hn)
                throw myException("Error while reading scaffold: Scaffold needs to have two haplotypes for every input sample");
        } catch (exception& e) {
            cout << e.what() << endl;
            exit(1);
        };

        break;

    default:
        assert(false); // this should not happen
    }
}

void Insti::CheckPanelPrereqs(PanelType panelType) {
    switch (panelType) {
    case PanelType::REFERENCE:
        m_bUsingRefHaps = true;
        assert(m_vRefHaps.size() == 0);
        break;

    case PanelType::SCAFFOLD:
        m_bUsingScaffold = true;
        assert(m_vuScaffoldPositions.size() == 0);
        assert(m_vScaffoldHaps.size() == 0);
        assert(m_vsScaffoldIDs.size() == 0);
        break;

    default:
        assert(false);  // this should not happen
    }
}

unsigned Insti::OpenLegend(string legendFile) {

    // read in legend file
    ifile legendFD(legendFile);
    string buffer;

    // discard header
    unsigned uLineNum = 0;

    while (getline(legendFD, buffer, '\n')) {
        uLineNum++;
        vector<string> tokens;
        sutils::tokenize(buffer, tokens);

        // make sure header start is correct
        if (uLineNum == 1) {
            vector<string> headerTokenized;
            string header = "id position a0 a1";
            sutils::tokenize(header, headerTokenized);

            for (unsigned i = 0; i != headerTokenized.size(); i++) {
                try {
                    if (tokens[i] != headerTokenized[i])
                        throw myException("Error in legend file (" + legendFile +
                                          "): header start does not match:\n\t" + header +
                                          "\n.  Instead the first line of the header is:\n\t" + buffer + "\n");
                } catch (exception& e) {
                    cout << e.what() << endl;
                    exit(2);
                }
            }

            continue;
        }

        // now make sure legend file matches sites
        try {
            if (tokens[1] != sutils::uint2str(site[uLineNum - 2].pos))
                throw myException("Error at Line " + sutils::uint2str(
                                      uLineNum) + " in legend file:\tPosition " + tokens[1] +
                                  " in legend file needs to match position in probin file: " + sutils::uint2str(
                                      site[uLineNum - 2].pos));
        } catch (exception& e) {
            cout << e.what() << endl;
            exit(1);
        }

    }

    // legend should have as many rows (sites) as in the probin file
    try {
        if (site.size() != uLineNum -1)
            throw myException("Alleles in legend file (" + sutils::uint2str(
                                  uLineNum -1 ) + ") need to match current data (" + sutils::uint2str(
                                  site.size()) + ")");
    } catch (exception& e) {
        cout << e.what() << endl;
        exit(1);
    }

    return uLineNum; // return unmber of sites
}

unsigned Insti::OpenHap(string hapFile, vector<uint64_t> & tempHaps) {

    // read in the hap file
    ifile hapsFD(hapFile);
    string buffer;
    unsigned uLineNum = 0;
    unsigned uNumHaps = 0;

    while (getline(hapsFD, buffer, '\n')) {
        vector<string> tokens;
        sutils::tokenize(buffer, tokens);

        // make sure header start is correct
        if (uLineNum == 0) {

            // count number of haps
            uNumHaps = tokens.size();

            // resize the vector for holding haps
            tempHaps.resize(uNumHaps * wn);
        }

        try {
            if (tokens.size() != uNumHaps)
                throw myException("Every row of hap file must have the same number of columns");
        } catch (exception& e) {
            cout << e.what() << endl;
            exit(1);
        }

        // store haplotypes
        for (unsigned i = 0; i != tokens.size(); i++) {
            int val = atoi(tokens[i].c_str());

            if (val == 0)
                set0(&tempHaps[i * wn], uLineNum);
            else if (val == 1)
                set1(&tempHaps[i * wn], uLineNum);
            else {
                cerr << "All alleles are not 0 or 1 In haplotypes file " << hapFile <<
                     " at line number " << uLineNum + 1 << endl;
                throw 1;
            }
        }

        uLineNum++;
    }

    if (uNumHaps == 0) {
        cerr << "num haps is 0.  Haps file empty?";
        throw 2;
    }

    return uNumHaps;
}

bool Insti::LoadHapLegSamp(string legendFile, string hapFile, string sampleFile,
                           PanelType panelType) {

    // make sure required files are defined
    if (legendFile.size() == 0) {
        cerr << "Need to define a legend file if defining a hap file\n";
        document();
    }

    if (hapFile.size() == 0) {
        cerr << "Need to define a hap file if defining a legend file\n";
        document();
    }

    if (sampleFile.size() == 0 && panelType == PanelType::SCAFFOLD) {
        cerr << "Need to define a sample file if using scaffold\n";
        document();
    }

    // for now sample loading is not implemented
    if (sampleFile.size() > 0) {
        cerr << "for now sample loading is not implemented in the sampleghap paradigm";
        document();
    }

    CheckPanelPrereqs(panelType);

    // just make sure the legend matches what we have for now
    OpenLegend(legendFile);

    vector<uint64_t> tempHaps;
    unsigned numHaps = OpenHap(hapFile, tempHaps);
    LoadHaps(tempHaps, numHaps, panelType);

    return true;
}

/* CHANGES from impute.cpp
   support for reference haplotypes
*/

void Insti::initialize() {

    // copied from SNPTools Impute
    // all haplotypes are saved in 64 bit unsigned ints (a word), where each bit represents a position
    // first, figure out how many words we'll need to store a hap and save in wn
    // if total sites overlaps with 00111111,  then wn = mn number of sites shifter to right....
    // we define a minimum block size of 64.
    wn = (mn & WordMod) ? (mn >> WordShift) + 1 : (mn >> WordShift);
    hn = in * 2; //number of haps
    haps.resize(hn * wn);  // space to store all haplotypes
    hnew.resize(hn *
                wn); // number of haplotypes = 2 * number of samples  ... haps mn is # of sites,
    hsum.assign(hn * mn,
                0); // one unsigned for every hap's site - what for?  To estimate allele probs

    //    pare.assign(in * in, 0);  // in x in matrix, one uint16 for every pair of individuals
    pn = 3 * mn;  // set the number of transitions.  three transitions for every site
    tran.resize(pn);  // tran looks like the transition matrix,

    // i.e. recombination rate
    vector<fast> temp(in *
                      pn);   // transitions between 3 types of genotypes P(RR), P(RA) P(AA)

    // initialize emission matrix
    // 4 emissions, for each site 4 emissions * number of samples.  (0|0, 1|1 0|1 1|0)
    en = 4 * mn;
    emit.resize(in * en);

    // is_par defines whether a site is in the paralogous region
    // if the site is on chromosome X
    // these magic numbers should probably get their own #define statement...
    // move away from vector of bools to vector of chars
    vector<bool> is_par(mn);

    if (posi.size() == mn) {
        for (unsigned m = 0; m < mn; m++)
            is_par[m] = (posi[m] >= 60001 && posi[m] <= 2699520) || (posi[m] >= 154931044
                        && posi[m] <= 155270560);
    }

    if (posi.size() != mn) {
        posi.resize(mn);

        for (unsigned m = 0; m < mn; m++) posi[m] = m;
    }    // if sites not stored

    // initialize the mutation rate mu:
    // If S is a harmonic series of length hn (number of haplotypes),
    // then mu = 1/S ( hn + 1/S)
    // initialize recombination rate rho based on SNP density
    fast mu = 0, rho;

    for (unsigned i = 1; i < hn; i++) mu += 1.0 / i;

    mu = 1 / mu;
    rho = 0.5 * mu * (mn - 1) / (posi[mn - 1] - posi[0]) / density;
    mu = mu / (hn + mu);  // rho is recombination rate?  mu is mutation rate

    // initialzie the site transition matrix tran
    // posi is recombination between its position and previous
    // r < 1; the larger the number of haplotypes, the smaller r gets
    // tran is a site's recombination probability matrix
    // r therefore must be a recombination rate estimate
    for (unsigned m = mn - 1; m; m--) {
        posi[m] = (posi[m] - posi[m - 1]) * rho;
        fast r = posi[m] / (posi[m] + hn);
        tran[m * 3] = (1 - r) * (1 - r);
        tran[m * 3 + 1] = r * (1 - r);
        tran[m * 3 + 2] = r *
                          r;  // for each position, transition.  r= alternative, 1-r= refrence? 4 state HMM with three transitions at each position
    }

    // initialize site mutation probability matrix
    // diagonal is chance of no mutation
    // the diagonal "rotated by 90 degrees" is the chance of both positions mutating
    // all other entries are chance of just one mutation
    pc[0][0] = pc[1][1] = pc[2][2] = pc[3][3] = (1 - mu) * (1 -
                                     mu); //	  probability of mutating no positions for each parents haplotype
    pc[0][1] = pc[0][2] = pc[1][0] = pc[1][3] = pc[2][0] = pc[2][3] = pc[3][1] =
                                         pc[3][2] = mu * (1 -
                                                 mu);  //	  probability of mutating one position for each parental haplotype
    pc[0][3] = pc[1][2] = pc[2][1] = pc[3][0] = mu *
                                     mu;  //	  probability of mutating both positions for each parental haplotype

    // initialize individual haplotypes
    for (unsigned i = 0; i < in; i++) {

        // define pointers to an individual's two haplotypes: ha and hb
        uint64_t *ha = &haps[i * 2 * wn], *hb = ha + wn;

        // here an individual's transition and probability matrices are
        // pulled out of the set of all individuals' matrices
        // t = genotype transition matrix? or genotype probability matrix?
        // e = phase emission matrix
        // p = a site's genotype probability (initially GL)
        // now iterate through each site
        fast *t = &temp[i * pn], *e = &emit[i * en], *p = &prob[i * 2];

        for (unsigned m = 0; m < mn; m++, t += 3, e += 4, p += hn) {

            // initialize genotype probabilities as genotype likelihoods
            if (is_x && male.find(name[i]) != male.end()
                    && !is_par[m]) { /// treat it differently for genders
                t[0] = max(1 - p[0] - p[1], 0.0f);
                t[1] = 0;
                t[2] = max(p[1], 0.0f);

                if (t[0] + t[2]) {
                    t[0] /= t[0] + t[2];
                    t[2] = 1 - t[0];
                } else t[0] = t[2] = 0.5;
            } else {

                // initial prob is the GL.
                t[0] = max(1 - p[0] - p[1], 0.0f);
                t[1] = max(p[0], 0.0f);
                t[2] = max(p[1], 0.0f);
            }

            // set each hap's bit randomly to 0 or 1 according to GLs, assuming flat prior for GLs
            // leave both 0 if we have a ref hom (most common)
            if (gsl_rng_uniform(rng) < t[0]);

            // choose one of the haps to be 1 and the other 0 if het
            // random phase
            else if (gsl_rng_uniform(rng) < t[1] / (t[1] + t[2])) {
                if (gsl_rng_uniform(rng) < 0.5) set1(ha, m);
                else set1(hb, m);
            }

            // hom
            else {
                set1(ha, m);
                set1(hb, m);
            }

            // initial emit assumes all states are by random mutation only.  basic state is ref/ref
            for (unsigned j = 0; j < 4; j++)
                e[j] = pc[j][0] * t[0] + pc[j][1] * t[1] + pc[j][2] * t[1] + pc[j][3] * t[2];
        }
    }

    swap(temp, prob);  // swap the assignments to each vector

    // end of copy from SNPTools Impute
    // load ref haps
    if (s_sRefLegendFile.size() > 0 || s_sRefHapFile.size() > 0)
        LoadHapLegSamp(s_sRefLegendFile, s_sRefHapFile, "", PanelType::REFERENCE);

    if (m_bUsingRefHaps) {

        // add ref haplotypes to sample haps
        haps.insert(haps.end(), m_vRefHaps.begin(), m_vRefHaps.end());

        // enlarge hnew so haps and hnew can be swapped
        // the ref haps are never updated, so they'll stick around forever
        hnew.insert(hnew.end(), m_vRefHaps.begin(), m_vRefHaps.end());
    }

    // load the scaffold
    if (s_scaffoldSampleFile.size() > 0 || s_scaffoldHapsFile.size() > 0)
        LoadHapsSamp(s_scaffoldSampleFile, s_scaffoldHapsFile, PanelType::SCAFFOLD);
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
fast Insti::solve(unsigned I, unsigned N, fast pen, Relationship &oRel) {

    // write log header
    stringstream message;
    message << "##iteration\tindividual\tproposal" << endl;
    WriteToLog(message.str());

    // pick 4 haplotype indices at random not from individual
    unsigned p[4];

    for (unsigned j = 0; j < 4; j++)
        p[j] = oRel.SampleHap(I, rng);

    // get a probability of the model for individual I given p
    fast curr = hmm_like(I, p);

    // pick a random haplotype to replace with another one from all
    // haplotypes.  calculate the new probability of the model given
    // those haplotypes.
    // accept new set if probability has increased.
    // otherwise, accept with penalized probability
    for (unsigned n = 0; n < N; n++) {  // fixed number of iterations
        unsigned rp = gsl_rng_get(rng) & 3, oh = p[rp];

        // kickstart phasing and imputation by only sampling haps
        // from ref panel in first round
        if (s_bKickStartFromRef && n == 0)
            p[rp] = oRel.SampleHap(I, rng, true);
        else
            p[rp] = oRel.SampleHap(I, rng);

        fast prop = hmm_like(I, p);
        bool bAccepted = false;

        if (prop > curr || gsl_rng_uniform(rng) < exp((prop - curr) * pen)) {
            curr = prop;
            bAccepted = true;
        } else p[rp] = oh;

        // Update relationship graph with proportion pen
        oRel.UpdateGraph(p, bAccepted, I, pen);

        // log accepted proposals
        if (bAccepted) {
            stringstream message;
            message << m_nIteration << "\t" << I << "\t" <<  prop << endl;
            WriteToLog(message.str());
        }
    }

    // if we have passed the burnin cycles (n >= bn)
    // start sampling the haplotypes for output
    /*
    if (P) {
        uint16_t *pa = &pare[I * in];
        for (unsigned i = 0; i < 4; i++) pa[p[i] / 2]++;
    }
    */
    hmm_work(I, p, pen);
    return curr;
}

/* CHANGES from impute.cpp
   added a member variable for n so
*/

void    Insti::estimate() {

    // the uniform relationship "graph" will be used until
    // -M option says not to.
    Relationship oUniformRel(2, in, hn + m_uNumRefHaps);

    switch (s_iEstimator) {

    // simulated annealing MH
    case 0:

        // cluster!
        if (Insti::s_uNumClusters > 0) {
            switch (Insti::s_uClusterType) {

            // use kmedoids
            case 0:
            case 1:
                m_oRelationship.init(3, haps, wn, mn, rng);
                break;

            // use kmeans
            case 2:
                if (UsingScaffold())
                    m_oRelationship.init(4, m_vScaffoldHaps, wn, m_uNumScaffoldHaps,
                                         s_scaffoldFreqCutoff);
                else
                    m_oRelationship.init(4, haps, wn, mn, rng);

                break;

            default:
                cout << "unexpected cluster type: " << s_uClusterType << endl;
                document();
            }
        }

        // just sample uniformly
        else
            m_oRelationship.init(2, in, hn + m_uNumRefHaps);

        break;

    case 1:
        estimate_EMC();
        return;

    case 2:
        estimate_AMH(0);
        return;

    case 3:
        estimate_AMH(1);
        return;

    default:
        document();
    }

    cerr.setf(ios::fixed);
    cerr.precision(3);
    cerr << "iter\tpress\tlike\tfold\trunTime\texpectedRunTime" << endl;

    // n is number of cycles = burnin + sampling cycles
    // increase penalty from 2/bn to 1 as we go through burnin
    // iterations.
    for (unsigned n = 0; n < bn + sn; n++) {
        m_nIteration = n;
        fast sum = 0, iter = 0;
        fast pen = min<fast> ((n + 1.0f) / Insti::s_uSABurninGen, 1);
        pen *= pen;  // pen = 1 after bn/2 iterations

        // Phase each individual based on the rest of the individuals
        for (unsigned i = 0; i < in; i++) {

            // re-update graph based on current haplotypes
            if (n < Insti::s_uStartClusterGen)
                sum += solve(i, m_uCycles, pen, oUniformRel);
            else
                sum += solve(i, m_uCycles, pen, m_oRelationship);

            iter += m_uCycles;
        }

        swap(hnew, haps);

        if (n >= bn) for (unsigned i = 0; i < in; i++) replace(i);  // call replace

        m_oRelationship.UpdateGraph(&haps);
        cerr << n << '\t' << pen << '\t' << sum / in / mn << '\t' << iter / in / in <<
             endl;
    }

    cerr << endl;
    result();    // call result
}

/* estimate_EMC -- Evolutionary Monte Carlo
   Here we try to increase the speed of convergence by
   running a parallel chain evolutionary monte carlo scheme.

   The reference for this implementation is
   "Advanced Markov Chain Monte Carlo Methods" by Liang, Liu and Carroll
   first edition?, 2010, pp. 128-132
*/

// solve(individual, number of cycles, penalty, burnin?)
fast Insti::solve_EMC(unsigned I, unsigned  N, fast S) {
    DEBUG_MSG("Entering solve_EMC..." << endl);

    // for lack of a better place, define free parameters here
    fast fMutationRate = 0.3f; // see p.134 of Liang et al.
    fast fSelectTemp = 10000;
    fast fMaxTemp = Insti::s_uParallelChains;

    // write log header
    stringstream message;
    message << "##iteration\tindividual\tproposal\tchainID\tchainTemp"
            << "\tmutation" << endl;
    WriteToLog(message.str());

    // initialize emc chains with increasing temperatures
    vector <EMCChain> vcChains;
    vector <uint>
    vuChainTempHierarchy; // index of Chains sorted by temperature, ascending

    for (unsigned i = 0; i < Insti::s_uParallelChains; i++) {
        vcChains.push_back(EMCChain((i + 1) * fMaxTemp / Insti::s_uParallelChains,
                                    fSelectTemp, I, in, i));
        DEBUG_MSG2("\tlength of vcChains\t" << vcChains.size() << endl);

        // initialize current likelihood
        // randomize parent haps
        for (unsigned j = 0; j < 4; j++) {
            do {
                unsigned uChosenHap = gsl_rng_get(rng) % vcChains[i].m_uHapNum;
                vcChains[i].setParent(j, uChosenHap);
                DEBUG_MSG3("\t\tchosen hap:\t" << uChosenHap << endl << "\t\tchosen parent:\t"
                           << vcChains[i].getParent(j) << endl);
            } while (vcChains[i].getParent(j) / 2 == vcChains[i].m_uI);
        }

        DEBUG_MSG2("\tsetting likelihood" << endl;);

        // set likelihood
        vcChains[i].setLike(hmm_like(vcChains[i].m_uI, vcChains[i].getParents()));
        vuChainTempHierarchy.push_back(i);
    }

    // pick a random haplotype to replace with another one from all
    // haplotypes.  calculate the new probability of the model given
    // those haplotypes.
    // accept new set if probability has increased.
    // otherwise, accept with penalized probability
    for (unsigned n = 0; n < N; n++) {  // fixed number of iterations
        DEBUG_MSG2("\tIteration " << n << endl);

        //  now choose whether to mutate or crossover
        bool bMutate =  gsl_rng_uniform(rng) > fMutationRate;
        fast prop;

        if (bMutate) {

            // mutate
            DEBUG_MSG2("\tMutating...");

            // choose chain randomly (uniform)
            unsigned j = gsl_rng_get(rng) % Insti::s_uParallelChains;
            DEBUG_MSG2("\t" << "Temp: " << vcChains[j].m_fTemp);
            fast curr = vcChains[j].getLike();

            // choose parent hap (rp) to mutate
            // replaced hap is stored in oh (Original Hap)
            unsigned rp = gsl_rng_get(rng) & 3, oh = vcChains[j].getParent(rp);

            // mutate parent hap
            do vcChains[j].setParent(rp, gsl_rng_get(rng) % hn);

            while (vcChains[j].getParent(rp) / 2 == I);

            // calculate acceptance probability
            prop = hmm_like(vcChains[j].m_uI, vcChains[j].getParents());

            if (prop > curr
                    || gsl_rng_uniform(rng) < exp((curr - prop) / vcChains[j].m_fTemp)) {
                vcChains[j].setLike(prop);
                WriteToLog(vcChains[j], bMutate);
            } else vcChains[j].setParent(rp, oh);
        } else {

            //crossover
            DEBUG_MSG2("\tCrossing Over...");

            // 1. choose random chain to work on by roulette wheel selection
            int iFirstChainIndex = RWSelection(vcChains);
            DEBUG_MSG2("\t\tFirst Chain:\t" << vcChains[iFirstChainIndex].m_uChainID <<
                       endl);
            DEBUG_MSG2("\t\tSelecting second chain:");

            // 2. select second chain at random (uniform) from remaining chains
            int iSecondChain;

            do {
                iSecondChain = gsl_rng_get(rng) % Insti::s_uParallelChains;
            } while (iSecondChain == iFirstChainIndex);

            //---- only crossover with certain probability depending
            //---- on crossover probabilities of parents
            // A. save parents and likelihoods of original chains
            EMCChain cFirstOrigChain = vcChains[iFirstChainIndex];
            EMCChain cSecondOrigChain = vcChains[iSecondChain];
            bool bFirstOrigChainHigherLike = cFirstOrigChain.getLike() >
                                             cSecondOrigChain.getLike();

            // B. cross over
            // uniform crossover: find haps to keep
            uint64_t cSelection = static_cast<uint64_t>(gsl_rng_get(rng) & 15);

            for (unsigned i = 0; i < 4; i++) {

                // if bit at location i is 1, exchange haps
                if (test(&cSelection, i)) {
                    unsigned oh = vcChains[iFirstChainIndex].getParent(i);
                    vcChains[iFirstChainIndex].setParent(i, vcChains[iSecondChain].getParent(i));
                    vcChains[iSecondChain].setParent(i, oh);
                }
            }

            // update likelihoods of crossed over chains
            auto& rcFirstChain = vcChains[iFirstChainIndex];
            rcFirstChain.setLike(hmm_like(rcFirstChain.m_uI, rcFirstChain.getParents()));
            auto& rcSecondChain = vcChains[iSecondChain];
            rcSecondChain.setLike(hmm_like(rcSecondChain.m_uI,
                                           rcSecondChain.getParents()));
            bool const bFirstChainHigherLike = cFirstOrigChain.getLike() >
                                               cSecondOrigChain.getLike();

            // C. deterimen if cross-over was successful
            bool bCrossAccepted = false;

            // the order of likelihoods is not the same
            if (bFirstOrigChainHigherLike ^ bFirstChainHigherLike) {
                bCrossAccepted = gsl_rng_uniform(rng)
                                 <= exp(((cSecondOrigChain.getLike() - rcFirstChain.getLike()) /
                                         cSecondOrigChain.m_fTemp)
                                        + (cFirstOrigChain.getLike() - rcSecondChain.getLike()) /
                                        cFirstOrigChain.m_fTemp);
            }

            // the order of the likelihoods matches
            else {
                bCrossAccepted = gsl_rng_uniform(rng)
                                 <= exp(((cSecondOrigChain.getLike() - rcSecondChain.getLike()) /
                                         cSecondOrigChain.m_fTemp)
                                        + (cFirstOrigChain.getLike() - rcFirstChain.getLike()) /
                                        cFirstOrigChain.m_fTemp);
            }

            // replace old chains if cross-over was not accepted
            if (!bCrossAccepted) {
                vcChains[iFirstChainIndex] = cFirstOrigChain;
                vcChains[iSecondChain] = cSecondOrigChain;
                stringstream message;
                message << "# Unsuccessful Crossover\tChainIDs:\t" << rcFirstChain.m_uChainID <<
                        "\t" << rcSecondChain.m_uChainID << endl;
                WriteToLog(message.str());
            } else {

                // otherwise log changes to likelihood
                WriteToLog(vcChains[iFirstChainIndex], bMutate);
                WriteToLog(vcChains[iSecondChain], bMutate);
            }
        }

        // now try Insti::s_uParallelChains exchanges
        DEBUG_MSG2("\tExchanging..." << endl);
        unsigned uNumExchanges = 0;

        for (unsigned i = 0; i < Insti::s_uParallelChains; i++) {
            unsigned uFirstChainIndex = gsl_rng_get(rng) % Insti::s_uParallelChains;
            DEBUG_MSG3("\t\tfirstChainIndex " << uFirstChainIndex);
            unsigned uFirstChainHierarchyIndex = vuChainTempHierarchy[ uFirstChainIndex ];
            DEBUG_MSG3("\tfirst chain: " << vuChainTempHierarchy[ uFirstChainIndex ]);

            // selecting second chain
            unsigned uSecondChainIndex;

            if (uFirstChainIndex == 0)
                uSecondChainIndex = uFirstChainIndex + 1;
            else if (uFirstChainIndex == Insti::s_uParallelChains - 1)
                uSecondChainIndex = Insti::s_uParallelChains - 2;
            else if (gsl_rng_get(rng) & 1)
                uSecondChainIndex = uFirstChainIndex - 1;
            else
                uSecondChainIndex = uFirstChainIndex + 1;

            unsigned uSecondCHI = vuChainTempHierarchy[ uSecondChainIndex ];
            DEBUG_MSG3("\tsecond chain: " << vuChainTempHierarchy[ uSecondChainIndex ]);

            // MH step for exchange
            fast fAcceptProb = min<fast> (exp((
                                                  vcChains[uFirstChainHierarchyIndex].getLike() - vcChains[uSecondCHI].getLike())
                                              * ((1 / vcChains[uFirstChainHierarchyIndex].m_fTemp) -
                                                 (1 / vcChains[uSecondCHI].m_fTemp)))
                                          , 1);
            DEBUG_MSG3("\taccept prob: " << fAcceptProb);

            // exchange with acceptance probability
            if (gsl_rng_uniform(rng) < fAcceptProb) {

                // exchange temperatures
                fast fTemp = vcChains[uFirstChainHierarchyIndex].m_fTemp;
                vcChains[uFirstChainHierarchyIndex].setTemp(vcChains[uSecondCHI].m_fTemp);
                vcChains[uSecondCHI].setTemp(fTemp);

                // exchange location in vcChains
                std::swap(vuChainTempHierarchy[uFirstChainIndex],
                          vuChainTempHierarchy[uSecondChainIndex]);
                ++ uNumExchanges;
            }

            DEBUG_MSG3("\tnumExchanges: " << uNumExchanges << endl);
        }

        // keep track of number of exchanges
        stringstream message;
        message << "# Number of Exchanges out of total:\t" << uNumExchanges <<
                "\t" << Insti::s_uParallelChains << endl;
        WriteToLog(message.str());
    }

    // now select a chain for sampling according to roulette wheel selection
    int iFirstChainIndex = RWSelection(vcChains);
    auto& rcFirstChain = vcChains[iFirstChainIndex];

    // if we have passed the burnin cycles (n >= bn)
    // start sampling the haplotypes for output
    /*
    if (P) {
        uint16_t *pa = &pare[I * in];
        for (unsigned i = 0; i < 4; i++) pa[rcFirstChain.getParent(i) / 2]++;
    }
    */
    DEBUG_MSG("Updating individual " << I << "\n");

    // update haplotypes of I
    //    cerr << "parents (2 lines):" << endl;
    //    cerr << rcFirstChain.m_auParents[1] << endl;
    //    cerr << vcChains[ iFirstChainIndex ].m_auParents[1] << endl;
    hmm_work(I, rcFirstChain.getParents(), S);
    return rcFirstChain.getLike();
}

/* estimate_EMC -- Evolutionary Monte Carlo
   Here we try to increase the speed of convergence by
   running a parallel chain evolutionary monte carlo scheme.

   The reference for this implementation is
   "Advanced Markov Choin Monte Carlo Methods" by Liang, Liu and Carroll
   first edition?, 2010, pp. 128-132
*/

void    Insti::estimate_EMC() {
    cerr.setf(ios::fixed);
    cerr.precision(3);
    cerr << "Running Evolutionary Monte Carlo\n";
    cerr << "iter\tpress\tlike\tfold\n";

    // n is number of cycles = burnin + sampling cycles
    // increase penalty from 2/bn to 1 as we go through burnin
    // iterations.
    for (unsigned n = 0; n < bn + sn; n++) {
        m_nIteration = n;
        fast sum = 0, pen = min<fast> (2 * (n + 1.0f) / bn, 1), iter = 0;
        pen *= pen;  // pen = 1 after bn/2 iterations

        for (unsigned i = 0; i < in; i++) {
            sum += solve_EMC(i, m_uCycles, pen);  // call solve=> inputs the sample number,
            iter += m_uCycles;
        }

        swap(hnew, haps);

        if (n >= bn) for (unsigned i = 0; i < in; i++) replace(i);  // call replace

        cerr << n << '\t' << pen << '\t' << sum / in / mn << '\t' << iter / in / in <<
             '\r';
    }

    cerr << endl;
    result();    // call result
}


/* estimate_AMH -- Adaptive Metropolis Hastings
   Here we try to increase the speed of convergence by
   keeping track of which individuals tend to copy from each other

   The reference for parts of this implementation is
   "Advanced Markov Choin Monte Carlo Methods" by Liang, Liu and Carroll
   first edition?, 2010, pp. 309
*/
void    Insti::estimate_AMH(unsigned uRelMatType) {
    cerr.setf(ios::fixed);
    cerr.precision(3);
    cerr << "Running Adaptive Metropolis Hastings\n";
    cerr << "iter\tpress\tlike\tfold" << endl;

    // create a relationshipGraph object
    // initialize relationship matrix
    // create an in x uSamplingInds matrix
    m_oRelationship.init(uRelMatType, in, hn + m_uNumRefHaps);

    // n is number of cycles = burnin + sampling cycles
    // increase penalty from 2/bn to 1 as we go through burnin
    // iterations.
    for (unsigned n = 0; n < bn + sn; n++) {

        //        cerr << "iter\t" << n << endl;
        m_nIteration = n;
        fast sum = 0, pen = min<fast> (2 * (n + 1.0f) / bn, 1), iter = 0;
        pen *= pen;  // pen = 1 after bn/2 iterations

        // update all individuals once
        for (unsigned i = 0; i < in; i++) {

            //            if( i % 1000 == 0)
            //                cerr << "cycle\t" << i << endl;
            sum += solve(i, m_uCycles, pen, m_oRelationship);
            iter += m_uCycles;
        }

        swap(hnew, haps);

        if (n >= bn) for (unsigned i = 0; i < in; i++) replace(i);  // call replace

        cerr << n << '\t' << pen << '\t' << sum / in / mn << '\t' << iter / in / in <<
             endl;
    }

    cerr << endl;
    result();    // call result
}

void Insti::save_relationship_graph(string sOutputFile) {
    vector<string> vsSampNames;
    vsSampNames.insert(vsSampNames.end(), name.begin(), name.end());

    for (unsigned i = 0; i < ceil(m_uNumRefHaps / 2);
            i++)  vsSampNames.push_back(string("refSamp") + sutils::uint2str(i));

    m_oRelationship.Save(sOutputFile, vsSampNames);
}

void    Insti::save_vcf(const char *F) {
    string temp = F;
    temp += ".vcf.gz";
    gzFile f = gzopen(temp.c_str(), "wt");
    gzprintf(f, "##fileformat=VCFv4.0\n");
    stringstream source;
    source <<  "##source=WTCHG:INSTIv" <<  VERSION_MAJOR << "." << VERSION_MINOR <<
           "." << VERSION_REVISION << "\n";
    gzprintf(f, source.str().c_str());

    //    gzprintf(f, "##reference=1000Genomes-NCBI37\n");
    gzprintf(f, "##iteration=%u\n", sn);
    gzprintf(f,
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");

    //    gzprintf(f, "##FORMAT=<ID=AP,Number=2,Type=Float,Description=\"Allelic Probability, P(Allele=1|Haplotype)\">\n");
    gzprintf(f,
             "##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Phred-scaled genotype posterior probability\">\n");
    gzprintf(f, "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");

    for (unsigned i = 0; i < in; i++) gzprintf(f, "\t%s", name[i].c_str());

    for (unsigned m = 0; m < mn; m++) {
        gzprintf(f, "\n%s\t%u\t.\t%c\t%c\t100\tPASS\t.\tGT:GP",
                 site[m].chr.c_str(), site[m].pos, site[m].all[0], site[m].all[1]);
        fast *p = &prob[m * hn];

        for (unsigned i = 0; i < in; i++, p += 2) {
            fast prr = (1 - p[0]) * (1 - p[1]), pra = (1 - p[0]) * p[1] + p[0] * (1 - p[1]),
                 paa = p[0] * p[1];

            if (prr >= pra && prr >= paa) gzprintf(f, "\t0|0");
            else if (pra >= prr && pra >= paa) {
                if (p[0] > p[1]) gzprintf(f, "\t1|0");
                else gzprintf(f, "\t0|1");
            } else gzprintf(f, "\t1|1");

            // test for any p being zero
            vector < fast > vfProb = { prr, pra, paa };

            for (unsigned j = 0; j < vfProb.size(); j++) {
                if (vfProb[j] == 0) vfProb[j] = 999;
                else if (vfProb[j] == 1) vfProb[j] = 0;
                else if (vfProb[j] > 1 || vfProb[j] < 0) {
                    cerr << "error: probability "  << vfProb[j] << " is not between 0 and 1";
                    exit(1);
                } else vfProb[j] = -10 * log10(vfProb[j]);
            }

            //-10 * log10(prr), -10 * log10(pra), -10 * log10(paa));
            gzprintf(f, ":%.3f,%.3f,%.3f", vfProb[0], vfProb[1], vfProb[2]);
        }
    }

    gzprintf(f, "\n");
    gzclose(f);
}


void    Insti::document(void) {
    cerr << "Author\tWarren W. Kretzschmar @ Marchini Group @ Universiy of Oxford - Statistics";
    cerr << "\n\nThis code is based on SNPTools impute:";
    cerr << "\nhaplotype imputation by cFDSL distribution";
    cerr << "\nAuthor\tYi Wang @ Fuli Yu' Group @ BCM-HGSC";
    cerr << "\n\nusage\timpute [options] 1.bin 2.bin ...";
    cerr << "\n\t-d <density>    relative SNP density to Sanger sequencing (1)";

    //    cerr << "\n\t-b <burn>       burn-in generations (56)";
    cerr << "\n\t-l <file>       list of input files";
    cerr << "\n\t-n <fold>       sample size*fold of nested MH sampler iteration (2)";

    //    cerr << "\n\t-t <thread>     number of threads (0=MAX)";
    cerr << "\n\t-v <vcf>        integrate known genotype in VCF format";
    cerr << "\n\t-c <conf>       confidence of known genotype (0.9998)";
    cerr << "\n\t-x <gender>     impute x chromosome data";
    cerr << "\n\t-e <file>       write log to file";
    cerr << "\n\n    GENERATION OPTIONS";
    cerr << "\n\t-m <mcmc>       sampling generations (200)";
    cerr << "\n\t-C <integer>    number of cycles to estimate an individual's parents before updating";
    cerr << "\n\t-B <integer>    number of simulated annealing generations (28)";
    cerr << "\n\t-i <integer>    number of non-simulated annealing burnin generations (28)";
    cerr << "\n\t-M <integer>    generation number at which to start clustering, 0-based (28)";
    cerr << "\n\n    HAPLOTYPE ESTIMATION OPTIONS";
    cerr << "\n\t-E <integer>    choice of estimation algorithm (0)";
    cerr << "\n\t                0 - Metropolis Hastings with simulated annealing";
    cerr << "\n\t                1 - Evolutionary Monte Carlo with -p parallel chains";
    cerr << "\n\t                2 - Adaptive Metropolis Hastings - sample/sample matrix";
    cerr << "\n\t                3 - Adaptive Metropolis Hastings - sample/haplotype matrix";
    cerr << "\n\t-p <integer>    number of parallel chains to use in parallel estimation algorithms";
    cerr << "\n\t                (at least 2, default 5)";
    cerr << "\n\t-K <integer>    number of clusters to use for haplotypes clustering (0 = option is off).";
    cerr << "\n\t                Does not currently work with -k option";
    cerr << "\n\t-t <integer>    Cluster type (0)";
    cerr << "\n\t                0 - k-Medoids -- PAM";
    cerr << "\n\t                1 - k-Medoids -- Park and Jun 2008";
    cerr << "\n\t                2 - k-Means -- IMPUTE2 (-K is the number of haplotypes to keep)";
    cerr << "\n\n    REFERENCE PANEL OPTIONS";
    cerr << "\n\t-H <file>       IMPUTE2 style HAP file";
    cerr << "\n\t-L <file>       IMPUTE2 style LEGEND file";
    cerr << "\n\t-k              Kickstart phasing by using only ref panel in first iteration";
    cerr << "\n\n    SCAFFOLD OPTIONS";
    cerr << "\n\t-h <file>       WTCCC style HAPS file";
    cerr << "\n\t-s <file>       WTCCC style SAMPLE file";
    cerr << "\n\t-q <float>      Minor allele frequency ([0-1], default 0.05) above which sites are fixed (using -c) according to scaffold.";
    cerr << "\n\n";
    exit(1);
}














































