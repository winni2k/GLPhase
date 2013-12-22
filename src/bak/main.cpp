//
//  main.cpp
//  wimpute
//
//  Created by warren kretzschmar on 19/04/2013.
//  Copyright (c) 2013 warren kretzschmar. All rights reserved.
//

// 

#include "impute.h"

int main(int ac, char **av) {
    Impute::bn = 56;
    Impute::sn = 200;
    Impute::nn = 2;
    Impute::density = 1.0;
    Impute::conf = 0.9998;
    Impute::is_x = false;
    Impute::is_y = false;
    uint threads = 0;
    vector <string> file;

    int opt;
    while ((opt = getopt(ac, av, "d:b:l:m:n:t:v:c:x:")) >= 0) {
        switch (opt) {
            case    'd':
                Impute::density = atof(optarg);
                break;
            case    'b':
                Impute::bn = atoi(optarg);
                break;
            case    'm':
                Impute::sn = atoi(optarg);
                break;
            case    'n':
                Impute::nn = atoi(optarg);
                break;
            case    't':
                threads = atoi(optarg);
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
            default:
                Impute::document();
        }
    }
    if (threads) omp_set_num_threads(threads);
    for (int i = optind; i < ac; i++) file.push_back(av[i]);
    sort(file.begin(), file.end());
    uint fn = unique(file.begin(), file.end()) - file.begin();
    if (!fn) Impute::document();

#pragma omp parallel for
    for (uint i = 0; i < fn; i++) {
        timeval sta, end;
        gettimeofday(&sta, NULL);
        Impute lp;
        if (!lp.load_bin(file[i].c_str())) {
            cerr << "fail to load " << file[i] << endl;
            continue;
        }
        for (uint j = 0; j < Impute::vcf_file.size(); j++)
            cerr << Impute::vcf_file[j] << '\t' << lp.load_vcf(Impute::vcf_file[j].c_str()) << endl;
        lp.initialize();
        lp.estimate();
        lp.save_vcf(file[i].c_str());
        lp.save_pare(file[i].c_str());
        char temp[256];
        sprintf(temp, "mv %s %s.ok", file[i].c_str(), file[i].c_str());
        gettimeofday(&end, NULL);
        cerr << "rename\t" << system(temp) << endl;
        cerr << "time\t" << end.tv_sec - sta.tv_sec + 1e-6 * (end.tv_usec - sta.tv_usec) << endl << endl;
    }
    return 0;
}

