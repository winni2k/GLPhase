//
//  main.cpp
//  insti
//
//  Created by warren kretzschmar on 19/04/2013.
//  Copyright (c) 2013 warren kretzschmar. All rights reserved.
//

// 

#include "version.hpp"
#include "impute.hpp"
#include "insti.hpp"
#include <chrono>

int main(int ac, char **av) {

    cerr<< "INSTI -- v" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_REVISION << endl;
    Impute::sn = 200;
    Impute::nn = 2;
    Impute::density = 1.0;
    Impute::conf = 0.9998;
    Impute::is_x = false;
    Impute::is_y = false;
    
//    uint threads = 0;
    vector <string> file;
    
    string sLogFile;
    int opt;
    while ((opt = getopt(ac, av, "d:l:m:n:v:c:x:e:E:p:C:L:H:kK:t:B:i:M:")) >= 0) {
        switch (opt) {
        case    'd':
            Impute::density = atof(optarg);
            break;
/*        case 'b':
            Impute::bn = atoi(optarg);
            break; */
        case    'm':
            Impute::sn = atoi(optarg);
            break;
        case    'n':
            Impute::nn = atoi(optarg);
            break;
        case    'v':
            Impute::vcf_file.push_back(optarg);
            break;
        case    'c':
            Impute::conf = atof(optarg);
            break;
        case    'x':
            Impute::is_x = true;
            Impute::gender(optarg);
            break;
        case    'l': {
            char temp[256];
            FILE *f = fopen(optarg, "rt");
            while (fscanf(f, "%s", temp) != EOF) file.push_back(temp);
            fclose(f);
        }
            break;
        case 'e':
            Insti::s_bIsLogging = true;
            sLogFile = optarg;
            break;
        case 'E':
            Insti::s_iEstimator = atoi(optarg);
            if(Insti::s_iEstimator > 3){
                cerr << "-E needs to be between 0 and 3" << endl;
                Insti::document();
            }
            break;
        case 'p':{           
            uint uP = atoi(optarg);
            if(uP < 2)
                Insti::document();
            Insti::s_uParallelChains = uP;
            break;
        }
        case 'C':
            Insti::s_uCycles = atoi(optarg);
            break;
        case 'L':
            Insti::s_sRefLegendFile = optarg;
            break;
        case 'H':
            Insti::s_sRefHapFile = optarg;
            break;
        case 'k':
            Insti::s_bKickStartFromRef = true;
            break;
/*        case    't':
            threads = atoi(optarg);
            break;
*/
        case 'K':
            Insti::s_uNumClusters = atoi(optarg);
            break;
        case 't':
            Insti::s_uClusterType = atoi(optarg);
            break;
        case 'B':
            Insti::s_uSABurninGen = atoi(optarg);
            break;
        case 'i':
            Insti::s_uNonSABurninGen = atoi(optarg);
            break;
        case 'M':
            Insti::s_uStartClusterGen = atoi(optarg);
            break;
        case 'h':
            Insti::s_sScaffoldHapFile = optarg;
            break;
        case 's':
            Insti::s_sScaffoldLegendFile = optarg;
            break;
        case 'S':
            Insti::s_sScaffoldSampleFile = optarg;
            break;
        case 'q':
            Insti::s_dScaffoldFreqCutoff = std::stod(optarg);
            break;
        default:
            Insti::document();
        }
    }

    // need to specify burnin generations as sum of SA and non-SA gens
    Impute::bn = Insti::s_uSABurninGen + Insti::s_uNonSABurninGen;
    
    // need to specify ref panel if kickstarting
    if(Insti::s_bKickStartFromRef){
        if( Insti::s_sRefLegendFile.size() == 0){
            cerr << endl << "error: Need to specify ref panel if kickstarting." << endl;
            Insti::document();
        }
    }
    
//    if (threads) omp_set_num_threads(threads);
    // read in files
    for (int i = optind; i < ac; i++) file.push_back(av[i]);
    sort(file.begin(), file.end());
    uint fn = unique(file.begin(), file.end()) - file.begin();
    if (!fn) Insti::document();

    // Die if more than one file was specified on command line
    if(fn != 1){
        cerr << endl << "INSTI only accepts one input .bin file" << endl << endl;
        Insti::document();
    }
    
//#pragma omp parallel for
    for (uint i = 0; i < fn; i++) {

        // keep track of time - these things are important!
        timeval sta, end;
        gettimeofday(&sta, NULL);

        // create an Insti instance!
        Insti lp;

        if( Insti::s_bIsLogging )
            lp.SetLog(sLogFile);

        // print date to start of log
        auto tt = std::chrono::system_clock::to_time_t ( std::chrono::system_clock::now());
        stringstream log;
        log << ctime(&tt) << endl;
        lp.WriteToLog(log.str());

        // load gls
        // add a reserve of space
        if (!lp.load_bin(file[i].c_str())) {
            cerr << "fail to load " << file[i] << endl;
            continue;
        }

        /* debugging
        cout << "name\tsite\tprob\tposi"<<endl;
        for ( auto x: lp.name) cout << x << " ";
        cout << endl;
        for ( auto x: lp.site) cout << x.chr << "\t" << x.all << "\t" << x.pos << endl;
        for ( auto x: lp.prob) cout << x << " ";
        cout << endl;
        for ( auto x: lp.posi) cout << x << " ";
        cout << endl;
        */
        for (uint j = 0; j < Impute::vcf_file.size(); j++)
            cerr << Impute::vcf_file[j] << '\t' << lp.load_vcf(Impute::vcf_file[j].c_str()) << endl;
        cerr << "initializing..\n";

        lp.initialize();
        cerr << "estimating..\n";

        // choose which estimation method to use
        lp.estimate();

        // save results of estimation
        lp.save_vcf(file[i].c_str());
        lp.save_relationship_graph(file[i]);
//        char temp[256];
//        sprintf(temp, "mv %s %s.ok", file[i].c_str(), file[i].c_str());
//        cerr << "rename\t" << system(temp) << endl;
        gettimeofday(&end, NULL);
        cerr << "time\t" << end.tv_sec - sta.tv_sec + 1e-6 * (end.tv_usec - sta.tv_usec) << endl << endl;
    }
    return 0;
}

