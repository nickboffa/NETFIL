#include <iostream>
#include <fstream>
#include <cstring>
#include "mda.h"
using namespace std;

int count_mda_scenarios(string filename){
    //counting the number of different MDA scenerios that we are testing!

    ifstream in;
    string line;
    
    in.open(filename.c_str());
    if(!in){
        cout << "open " << filename << " failed" << endl;
        exit(1);
    }
    int MDAScenarioNum {-1}; //Count the number of lines (not including the first, which is just the column titles)

    while (getline(in, line)){
        ++MDAScenarioNum;
    }
    in.close();
    return MDAScenarioNum;
}

mda_strat get_mda_strat(string filename, int N)
{
    ifstream in;
    
    in.open(filename.c_str());
    if(!in){
        cout << "open " << filename << " failed" << endl;
        exit(1);
    }
    
    string line;

    //skip N lines
    for(int i = 0; i < N; ++i){
        getline(in, line);
    }
    
    getline(in,line);

    // The below could maybe be modernised using
    char *str = new char[line.size()+1];
    strcpy(str, line.c_str());
    char *p = NULL;

    p = strtok(str, ",");       double MDACoverage = atof(p);
    p = strtok(NULL, ",");      double MDAKillProb = atof(p);
    p = strtok(NULL, ",");      double MDAFullSterProb = atof(p);
    p = strtok(NULL, ",");      double MDAPartSterProb = atof(p);
    p = strtok(NULL, ",");      double MDASterDur = atof(p);
    p = strtok(NULL, ",");      double MDAPartSterMagnitude = atof(p);
    p = strtok(NULL, ",");      int MinAge = atoi(p);
    p = strtok(NULL, ",");      int MDAStartYear = atoi(p);
    p = strtok(NULL, ",");      int MDANumRound = atoi(p);
    p = strtok(NULL, ",");      int MDAYearsBetweenRound = atoi(p);
    p = strtok(NULL, ",");      int NumSims = atoi(p);
   
    delete []str;
    
    //creating drug profile
    drugs drug {MDAKillProb, MDAFullSterProb, MDAPartSterProb, MDASterDur, MDAPartSterMagnitude};
   
    //creating mda strate
    mda_strat strat {MDACoverage,drug, MinAge, MDAStartYear, MDANumRound, MDAYearsBetweenRound,NumSims};

    //printing MDA and drug stats
    //drug.print_drugs();
    //strat.print_mda_strat();

    return strat;
}