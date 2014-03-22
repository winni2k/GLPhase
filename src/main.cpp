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

using namespace std;

int main(int ac, char **av) {

  cerr << "INSTI -- v" << VERSION_MAJOR << "." << VERSION_MINOR << "."
       << VERSION_REVISION << endl;
  try {
    stringstream commandLine;
    commandLine << av[0];
    if (ac) {
      for (int i = 1; i < ac; ++i) {
        commandLine << " " << av[i];
      }
    }

    Impute::sn = 200;
    Impute::nn = 2;
    Impute::density = 1.0;
    Impute::conf = 0.9998;
    Impute::is_x = false;
    Impute::is_y = false;

    //    uint threads = 0;
    vector<string> file;
    string outBase;

    string sLogFile;
    int opt;
    while ((opt = getopt(
                ac, av, "Vd:l:m:n:v:c:x:e:E:p:C:L:H:kK:t:B:i:M:h:s:q:fo:DT")) >=
           0) {
      switch (opt) {
      case 'd':
        Impute::density = atof(optarg);
        break;
      /*        case 'b':
                  Impute::bn = atoi(optarg);
                  break; */
      case 'm':
        Impute::sn = atoi(optarg);
        break;
      case 'n':
        Impute::nn = atoi(optarg);
        break;
      case 'v':
        Impute::vcf_file.push_back(optarg);
        break;
      case 'c':
        Impute::conf = atof(optarg);
        break;
      case 'x':
        Impute::is_x = true;
        Impute::gender(optarg);
        break;
      case 'l': {
        char temp[256];
        FILE *f = fopen(optarg, "rt");
        while (fscanf(f, "%s", temp) != EOF)
          file.push_back(temp);
        fclose(f);
      } break;
      case 'e':
        Insti::s_bIsLogging = true;
        sLogFile = optarg;
        break;
      case 'E':
        Insti::s_iEstimator = atoi(optarg);
        if (Insti::s_iEstimator > 3) {
          cerr << "-E needs to be between 0 and 3" << endl;
          Insti::document();
        }
        break;
      case 'D':
        Insti::s_MHSamplerType = MHType::DRMH;
        break;
      case 'p': {
        uint uP = atoi(optarg);
        if (uP < 2)
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
        Insti::s_scaffoldHapsFile = optarg;
        break;
      case 's':
        Insti::s_scaffoldSampleFile = optarg;
        break;
      case 'q':
        Insti::s_scaffoldFreqCutoff = std::stod(optarg);
        break;
      case 'f':
        Insti::s_initPhaseFromScaffold = true;
        break;
      case 'o':
        outBase = optarg;
        break;
      case 'V':
        exit(0);
        break;
      case 'T':
        Insti::s_clusterDistanceMetric = kNNDistT::tracLen;
        break;
      default:
        Insti::document();
      }
    }

    cerr << "Call: " << commandLine.str() << endl;

    // need to specify burnin generations as sum of SA and non-SA gens
    Impute::bn = Insti::s_uSABurninGen + Insti::s_uNonSABurninGen;

    // need to specify ref panel if kickstarting
    if (Insti::s_bKickStartFromRef) {
      if (Insti::s_sRefLegendFile.size() == 0) {
        cerr << endl << "error: Need to specify ref panel if kickstarting."
             << endl;
        Insti::document();
      }
    }

    //    if (threads) omp_set_num_threads(threads);
    // read in files
    for (int i = optind; i < ac; i++)
      file.push_back(av[i]);
    sort(file.begin(), file.end());
    uint fn = unique(file.begin(), file.end()) - file.begin();
    if (!fn)
      cerr << "input files are not unique";

    // Die if more than one file was specified on command line
    if (fn != 1) {
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

      if (Insti::s_bIsLogging)
        lp.SetLog(sLogFile);

      // print date to start of log
      auto tt = std::chrono::system_clock::to_time_t(
          std::chrono::system_clock::now());
      stringstream log;
      log << ctime(&tt) << endl;
      lp.WriteToLog(log.str());

      // load gls
      // add a reserve of space
      if (!lp.load_bin(file[i].c_str())) {
        cerr << "fail to load " << file[i] << endl;
        continue;
      }

      for (uint j = 0; j < Impute::vcf_file.size(); j++)
        cerr << Impute::vcf_file[j] << '\t'
             << lp.load_vcf(Impute::vcf_file[j].c_str()) << endl;
      cerr << lp.m_tag << ": initializing..\n";

      lp.initialize();
      cerr << lp.m_tag << ": estimating..\n";

      // choose which estimation method to use
      lp.estimate();

      // save results of estimation
      if (outBase.empty())
        outBase = file[i];

      lp.save_vcf(outBase.c_str(), commandLine.str());
      lp.save_relationship_graph(outBase);

      // printing out run time
      gettimeofday(&end, NULL);
      cerr << lp.m_tag << ": time\t"
           << end.tv_sec - sta.tv_sec + 1e-6 * (end.tv_usec - sta.tv_usec)
           << endl << endl;
    }
  }
  catch (exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
