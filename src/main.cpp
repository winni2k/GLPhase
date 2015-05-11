//
//  main.cpp
//  insti
//
//  Created by warren kretzschmar on 19/04/2013.
//  Copyright (c) 2013 warren kretzschmar. All rights reserved.
//

//

#include <chrono>
#include "boost/program_options.hpp"
#include "version.hpp"
#include "impute.hpp"
#include "insti.hpp"

using namespace std;
namespace po = boost::program_options;

int main(int ac, char **av) {

  cout << "INSTI -- v" << VERSION_MAJOR << "." << VERSION_MINOR << "."
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

    string inputFileType = "bin";
    string inputFile;
    string outBase;
    string sexFile;

    string sLogFile;
    bool optMSet = false;

    po::options_description general("General options");
    general.add_options()(

        "help", "Print help messages")(

        "gls", po::value<string>(&inputFile)->required(),
        "Genotype likelihoods file ('bin','VCF','BCF' format)")(

        "fold,n", po::value<unsigned>(&Impute::nn)->default_value(2),
        "Fold: Number of iterations of nested MH sampler = sample size*fold")(

        "prefix,o", po::value<string>(&outBase),
        "Prefix to use for output files")(

        "threads,P", po::value<size_t>(&init.numThreads)->default_value(1),
        "Number of threads (0=MAX)")(

        "confidence,c", po::value<double>(&Impute::conf)->default_value(0.9998),
        "Confidence of known genotype")(

        "sex,x", po::value<string>(&sexFile),
        "Sex file. Impute x chromosome data")(

        "input-file-type,I",
        po::value<string>(&inputFileType)->default_value("bin"),
        "Input file type")(

        "log-file,e", po::value<string>(&sLogFile), "Write to log file")(

        "gmap,g", po::value<string>(&init.geneticMap), "Genetic map")(

        "gmap-inflation",
        po::value<double>(&init.geneticMapInflationFactor)->default_value(1),
        "Genetic map distance inflation factor")(

        "region,R", po::value<string>(&init.glSubsetRegion),
        string("Region to subset GLs to. Format: <chrom:start-end>. Only the "
               "first " +
               to_string(NUMSITES) +
               " sites are used.  A fatal error is thrown "
               "if more than 5% of sites in a region are "
               "thrown out in this way.").c_str());

    po::options_description generation("Generation options");
    generation.add_options()(

        "sampling-gens,m", po::value<unsigned>(&Impute::sn)->default_value(200),
        "Number of sampling generations")(

        "simulated-annealing-burnin-gens,B",
        po::value<unsigned>(&Insti::s_uSABurninGen)->default_value(28),
        "Number of simulated annealing burnin generations")(

        "burnin-gens,i",
        po::value<unsigned>(&Insti::s_uNonSABurninGen)->default_value(28),
        "Number of non-simulated annealing burnin generations")(

        "nested-gens,C",
        po::value<unsigned>(&Insti::s_uNonSABurninGen)->default_value(0),
        "Number of generations of nested MH sampler. Overrides --fold if "
        "greater than 0.");

    po::options_description estimates("Haplotype estimation options");
    estimates.add_options()(
        "hap-estimation-algorithm,E",
        po::value<unsigned>(&init.estimator)->default_value(0),
        "Haplotype estimation algorithm to use:\n"
        " 0: \tMetropolis Hastings with simulated "
        "annealing\n"
        " 1: \tEvolutionary Monte Carlo with "
        "--num-parallel-chains parallel chains\n"
        " 2: \tAdaptive Metropolis Hastings - "
        "sample/sample matrix\n"
        " 3: \tAdaptive Metropolis Hastings - "
        "sample/haplotype matrix")(

        "num-parallel-chains", po::value<unsigned>(&Insti::s_uParallelChains),
        "Number of parallel chains to use in parallel estimation algorithms")(

        "num-clusters,K",
        po::value<unsigned>(&Insti::s_uNumClusters)->default_value(0),
        "Number of clusters/nearest neighbors to use for haplotype "
        "clustering/kNN search. Does not currently work with --kickstart "
        "option.")(

        "cluster-type,t",
        po::value<unsigned>(&Insti::s_uClusterType)->default_value(0),
        "Cluster type:\n"
        " 0: \tk-Medoids -- PAM\n"
        " 1: \tk-Medoids -- Park and Jun 2008\n"
        " 2: \tk-Nearest Neighbors -- IMPUTE2 (-K is the "
        "number of haplotypes to keep)")(

        "use-tract-length,T", "Use shared tract length as distance metric for "
                              "clustering")(

        "recluster-every-n-gen,r", po::value<size_t>(&init.reclusterEveryNGen),
        "Recluster every n generations. Only works if cluster-type = 2")(

        "cluster-start-gen,M",
        po::value<unsigned>(&Insti::s_uStartClusterGen)->default_value(28),
        "Generation number (0-based) at which to start clustering")(

        "delayed-rejection-MH", "Use delayed reject Metropolis Hastings");

    po::options_description refPanel("Reference panel options");
    refPanel.add_options()("i2-haps,H",
                           po::value<string>(&Insti::s_sRefHapFile),
                           "IMPUTE2 style HAP file")(

        "i2-legend,L", po::value<string>(&Insti::s_sRefLegendFile),
        "IMPUTE2 style LEGEND file")(

        "kickstart,k", "Kickstart phasing by using only ref panel in "
                       "first iteration");

    po::options_description preHaps("Pre-existing haplotype options");
    preHaps.add_options()(

        "wtccc-haps,h", po::value<string>(&init.scaffoldHapsFile),
        "WTCCC style HAPS file")(

        "wtccc-samples,s", po::value<string>(&init.scaffoldSampleFile),
        "WTCCC style SAMPLES file")(

        "varaf-LB,q",
        po::value<double>(&init.scaffoldFreqLB)->default_value(0.0f),
        "Lower bound of variant allele frequency [0-1] above which sites are "
        "used for clustering from scaffold.")(

        "varaf-UB,Q",
        po::value<double>(&init.scaffoldFreqUB)->default_value(1.0f),
        "Upper bound of variant allele frequency [0-1] above which sites are "
        "used for clustering from scaffold.")(

        "use-minor-varaf,a",
        "Use minor allele frequency instead of variant "
        "allele frequency for clustering and applying --varaf-LB and "
        "--varaf-UB below which sites are used for clustering from scaffold.")(

        "fix-phase", "Fix phase according to pre-existing haplotypes");

    po::options_description all("allowed options");
    all.add(general).add(generation).add(estimates).add(refPanel).add(preHaps);

    po::positional_options_description positionalOptions;
    positionalOptions.add("gls", 1);

    po::variables_map vm;
    try {
      po::store(po::command_line_parser(ac, av)
                    .options(all)
                    .positional(positionalOptions)
                    .run(),
                vm); // can throw

      /** --help option
       */
      if (vm.count("help")) {
        cerr << general << endl;
        cerr << generation << endl;
        cerr << estimates << endl;
        cerr << refPanel << endl;
        cerr << preHaps << endl;
        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
                      // there are any problems
    } catch (po::error &e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;

      cerr << general << endl;
      cerr << generation << endl;
      cerr << estimates << endl;
      cerr << refPanel << endl;
      cerr << preHaps << endl;

      return 1;
    }

    // perform extra validations
    if (vm.count("sex"))
      Impute::is_x = true;

    if (vm.count("log-file"))
      Insti::s_bIsLogging = true;

    if (init.estimator > 3) {
      cerr << "ERROR: --hap-estimation-algorithm needs to be between 0 and 3"
           << endl;
      exit(1);
    }

    if (Insti::s_uParallelChains < 2) {
      cerr << "ERROR: --num-parallel-chains needs to be > 1" << endl;
      exit(1);
    }

    Insti::s_bKickStartFromRef = vm.count("kickstart") > 0;
    optMSet = vm.count("cluster-start-gen") > 0;
    init.scaffoldUsingMAF = vm.count("use-minor-varaf") > 0;
    init.initPhaseFromScaffold = vm.count("fix-phase");
    if (vm.count("use-tract-length"))
      Insti::s_clusterDistanceMetric = kNNDistT::tracLen;

    /*
    while (
        (opt = getopt(
             ac, av,
             "Vm:n:v:c:x:e:E:p:C:L:H:kK:t:B:i:M:h:s:q:Q:fo:DTr:P:ag:I:R:")) >=
        0) {
      switch (opt) {
    */
    /*      case 'd':
            Impute::density = atof(optarg);
            break;
                  case 'b':
                      Impute::bn = stoul(optarg);
                      break; */
    /*
  case 'm':
    Impute::sn = stoul(optarg);
    break;
  case 'n':
    Impute::nn = stoul(optarg);
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
    break;*/
    /*
      This option makes no sense if we are only allowing one input file
  case 'l': {
    char temp[256];
    FILE *f = fopen(optarg, "rt");
    while (fscanf(f, "%s", temp) != EOF)
      file.push_back(temp);
    fclose(f);
  } break;
    */
    /*  case 'I':
        inputFileType = optarg;
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
        init.scaffoldHapsFile = optarg;
        break;
      case 's':
        init.scaffoldSampleFile = optarg;
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
      case 'R':
        init.geneticMapInflationFactor = std::stoul(optarg);
        break;
      default:
        Insti::document();
      }
    }

    */

    cout << "Call: " << commandLine.str() << endl;

    // need to specify burnin generations as sum of SA and non-SA gens
    Impute::bn = Insti::s_uSABurninGen + Insti::s_uNonSABurninGen;
    if (optMSet == false)
      Insti::s_uStartClusterGen = Insti::s_uSABurninGen;

    // need to specify ref panel if kickstarting
    if (Insti::s_bKickStartFromRef) {
      if (Insti::s_sRefLegendFile.size() == 0) {
        cerr << endl << "error: Need to specify ref panel if kickstarting."
             << endl;
        Insti::document();
      }
    }

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

    // load gls
    if (inputFileType == "bin") {
      try {
        lp.load_bin(inputFile);
      } catch (std::exception &e) {
        cerr << "[main] While loading .bin file: " << inputFile << endl
             << e.what() << endl;
        exit(1);
      }
    } else if (inputFileType == "b") {
      try {
        lp.load_bcf(inputFile);
      } catch (exception &e) {
        cerr << "[main] While loading bcf/vcf file: " << inputFile << endl
             << e.what() << endl;
        exit(1);
      }
    } else
      throw runtime_error("[main] Unexpected input file type encountered: [" +
                          inputFileType + "]");

    /*
    load_vcf is broken
    for (uint j = 0; j < Impute::vcf_file.size(); j++)
      cerr << Impute::vcf_file[j] << '\t'
           << lp.load_vcf(Impute::vcf_file[j].c_str()) << endl;
    */
    cout << lp.m_tag << ": initializing.." << endl;

    lp.initialize();
    cout << lp.m_tag << ": estimating.." << endl;

    // choose which estimation method to use
    lp.estimate();

    // save results of estimation
    if (outBase.empty())
      outBase = inputFile;

    lp.save_vcf(outBase, commandLine.str());

    // output relationship graph if applicable
    if (init.estimator == 2 || init.estimator == 3) {
      try {
        lp.save_relationship_graph(outBase);
      } catch (exception &e) {
        cerr << e.what() << endl;
      }
    }

    // printing out run time
    gettimeofday(&end, NULL);
    cout << lp.m_tag << ": time\t"
         << end.tv_sec - sta.tv_sec + 1e-6 * (end.tv_usec - sta.tv_usec)
         << endl;
  } catch (exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

#ifndef NCUDA
  // this is here so I can use the profiler...
  cudaDeviceReset();
#endif

  return 0;
}
