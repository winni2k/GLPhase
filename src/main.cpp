//
//  main.cpp
//  insti
//
//  Created by warren kretzschmar on 19/04/2013.
//  Copyright (c) 2013 warren kretzschmar. All rights reserved.
//

//

#include <chrono>
#include "version.hpp"
#include "impute.hpp"
#include "insti.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int ac, char **av) {

  cerr << "INSTI -- v" << VERSION_MAJOR << "." << VERSION_MINOR << "."
       << VERSION_XSTR(VERSION_REVISION) << endl;
  try {
    stringstream commandLine;
    commandLine << av[0];
    if (ac) {
      for (int i = 1; i < ac; ++i) {
        commandLine << " " << av[i];
      }
    }

    // init is used to pass initialization options to Insti
    InstiHelper::Init init;

    Impute::sn = 200;
    Impute::nn = 2;
    Impute::density = 1.0; // this is not used anymore
    Impute::conf = 0.9998;
    Impute::is_x = false;
    Impute::is_y = false;

    string outBase;

    string sLogFile;
    int opt;
    bool optMSet = false;
    while ((opt = getopt(ac, av, "V:m:n:v:c:x:e:E:p:C:L:H:kK:t:B:i:M:h:s:q:Q:S:"
                                 "fo:DTr:P:ag:F:R:O:u:w:W:AI:")) >= 0) {
      switch (opt) {

      /*      case 'd':
              Impute::density = atof(optarg);
              break;
                    case 'b':
                        Impute::bn = stoul(optarg);
                        break; */
      case 'm':
        Impute::sn = stoul(optarg);
        break;
      case 'n':
        Impute::nn = stoul(optarg);
        break;
      /*
            case 'v':
                Impute::vcf_file.push_back(optarg);
                break;
            case 'c':
              Impute::conf = atof(optarg);
              break;
      */

      case 'x':
        Impute::is_x = true;
        Impute::gender(optarg);
        break;
      case 'e':
        Insti::s_bIsLogging = true;
        sLogFile = optarg;
        break;
      case 'E':
        init.estimator = stoul(optarg);
        if (init.estimator > 3) {
          cerr << "-E needs to be between 0 and 3" << endl;
          Insti::document();
        }
        break;

      case 'D':
        Insti::s_MHSamplerType = MHType::DRMH;
        break;
      case 'p': {
        uint uP = stoul(optarg);
        if (uP < 2)
          Insti::document();
        Insti::s_uParallelChains = uP;
        break;
      }
      case 'C':
        Insti::s_uCycles = stoul(optarg);
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
      case 'P':
        init.numThreads = stoul(optarg);
        break;
      case 'K':
        Insti::s_uNumClusters = stoul(optarg);
        break;
      case 't':
        Insti::s_uClusterType = stoul(optarg);
        break;
      case 'B':
        Insti::s_uSABurninGen = stoul(optarg);
        break;
      case 'i':
        Insti::s_uNonSABurninGen = stoul(optarg);
        break;
      case 'M':
        Insti::s_uStartClusterGen = stoul(optarg);
        optMSet = true;
        break;
      case 'h':
        init.scaffoldFiles["h"] = optarg;
        break;
      case 's':
        init.scaffoldFiles["s"] = optarg;
        break;
      case 'l':
        init.scaffoldFiles["l"] = optarg;
        break;
      case 'q':
        init.scaffoldFreqLB = std::stod(optarg);
        break;
      case 'Q':
        init.scaffoldFreqUB = std::stod(optarg);
        break;
      case 'a':
        init.scaffoldUsingMAF = true;
        break;
      case 'f':
        init.initPhaseFromScaffold = true;
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
      case 'r':
        init.reclusterEveryNGen = stoul(optarg);
        break;
      case 'g':
        init.geneticMap = optarg;
        break;
      case 'F':
        init.inputGLFileType = optarg;
        break;
      case 'R':
        init.inputGLRegion = optarg;
        break;
      case 'O':
        init.scaffoldExtraRegionSize = stoul(optarg);
        break;
      case 'S':
        init.inputSamplesKeep = optarg;
        break;
      case 'u':
        init.reclusterStage = optarg;
        if (init.reclusterStage.size() != 1) {
          cerr << "-u needs to be a single char" << endl;
          exit(1);
        }
        break;
      case 'w':
        init.fixPhaseFromScaffold = true;
        init.fixPhaseFrequencyLowerBound = stod(optarg);
        break;
      case 'W':
        init.fixPhaseFromScaffold = true;
        init.fixPhaseFrequencyUpperBound = stod(optarg);
        break;
      case 'A':
        init.fixPhaseFromScaffold = true;
        init.fixPhaseReferenceAlleleFrequency = false;
        break;
      case 'I':
        init.fixPhaseFromScaffold = true;
        init.fixPhaseAlwaysKeepVariantsFile = optarg;
        break;
      default:
        Insti::document();
      }
    }

    cerr << "Call: " << commandLine.str() << endl;

    // need to specify burnin generations as sum of SA and non-SA gens
    Impute::bn = Insti::s_uSABurninGen + Insti::s_uNonSABurninGen;
    if (optMSet == false)
      Insti::s_uStartClusterGen = Insti::s_uSABurninGen;
    clog << "Clustering will start at generation " << Insti::s_uStartClusterGen
         << endl;

    // need to specify ref panel if kickstarting
    if (Insti::s_bKickStartFromRef) {
      if (Insti::s_sRefLegendFile.size() == 0) {
        cerr << endl << "error: Need to specify ref panel if kickstarting."
             << endl;
        Insti::document();
      }
    }

    // Die if more than one file was specified on command line
    if (ac - optind != 1) {
      cerr << endl << "INSTI only accepts one input GL file" << endl << endl;
      Insti::document();
    }

    init.inputGLFile = av[optind];

    // keep track of time - these things are important!
    timeval sta, end;
    gettimeofday(&sta, NULL);

    // create an Insti instance!
    Insti lp(init);

    if (Insti::s_bIsLogging)
      lp.SetLog(sLogFile);

    // print date to start of log
    auto tt =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    stringstream log;
    log << "##" << ctime(&tt) << endl;
    lp.WriteToLog(log.str());

    /*
    impute::load_vcf needs to be reimplemented
  for (uint j = 0; j < Impute::vcf_file.size(); j++)
    cerr << Impute::vcf_file[j] << '\t'
         << lp.load_vcf(Impute::vcf_file[j].c_str()) << endl;
  */
    cerr << lp.m_tag << ": initializing..\n";

    lp.initialize();
    cerr << lp.m_tag << ": estimating..\n";

    // choose which estimation method to use
    lp.estimate();

    // save results of estimation
    if (outBase.empty())
      outBase = init.inputGLFile;

    outBase += ".phased";
    lp.save_vcf(outBase.c_str(), commandLine.str());
    try {
      lp.save_relationship_graph(outBase);
    }
    catch (exception &e) {
      cerr << e.what() << endl;
    }

    // printing out run time
    gettimeofday(&end, NULL);
    cerr << lp.m_tag << ": time\t"
         << end.tv_sec - sta.tv_sec + 1e-6 * (end.tv_usec - sta.tv_usec) << endl
         << endl;
  }
  catch (const exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  catch (...) {
    cerr << "Caught exception of unknown type" << endl;
    exit(2);
  }

#ifndef NCUDA
  // this is here so I can use the profiler...
  cudaDeviceReset();
#endif

  return 0;
}
