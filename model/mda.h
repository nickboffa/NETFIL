#ifndef mda_h
#include "params.h"
#define mda_h
#include <vector>
#include <algorithm>

using namespace std;

class MDAStrat;
class Drugs;


int count_mda_scenarios(string filename);
MDAStrat get_mda_strat(string filename, int N);


class Drugs{ //Drug profile for MDA
public:
    double kill_prob; // prob will kill worms
    double full_ster_prob; // prob will fully sterilise
    double part_ster_prob; // prob will partially sterilise
    double ster_dur; // duration of sterilisation
    double part_ster_magnitude; //magnitude of partial sterilisation

    
    Drugs(double K=0.0, double FSP=0.0, double PSP=0.0, double SD=0.0, double PSM=0.0){
        kill_prob = K;
        full_ster_prob = FSP;
        part_ster_prob = PSP;
        ster_dur = SD;
        part_ster_magnitude = PSM;
    }
    
    void print_drugs(ostream& out = cout){
        out << "Drug's kill_prob: " << kill_prob << endl;
        out << "Drug's full_ster_prob: " << full_ster_prob << endl;
        out << "Drug's part_ster_prob: " << part_ster_prob << endl;
        out << "Drug's ster_dur: " << ster_dur << endl;
        out << "Drug's part_ster_magnitude: " << part_ster_magnitude << endl;
    }
    
};

class MDAStrat{
public:
    double coverage; // MDA coverage (whole population, not just target)
    Drugs drug; // Drug profile
    int min_age; // Min age to take MDA
    int mda_start_year; // MDA start year
    int n_mda_rounds; // Number of rounds
    int years_between_rounds; // Time between rounds
    vector<int> mda_years; // Vector storing all mda years
    int n_sims; // Number of simulations

    MDAStrat(double C, Drugs D, int MN, int S, int N, int Y, int NS){
        coverage = C;
        drug = D;
        min_age = MN;
        mda_start_year = S;
        n_mda_rounds = N;
        years_between_rounds = Y;
        n_sims = NS;

        mda_years.resize(n_mda_rounds);
        for(int i = 0; i<n_mda_rounds; ++i){
            mda_years[i] = mda_start_year + i * years_between_rounds;
            cout << mda_years[i] << " is a MDA year" << endl;
        }   
    }

    void print_mda_strat(ostream& out = cout){
        out << "Coverage: " << coverage * 100 << "%" << endl;
        out << "Minimum age: " << min_age << endl;
        out << "First year of MDA: " << mda_start_year << endl;
        out << "Number of rounds: " << n_mda_rounds << endl;
        out << "Years between rounds: " << years_between_rounds << endl;
        out << "Number of simulations: " << n_sims << endl;
        drug.print_drugs(out);
    }

    bool is_mda_year(int Year){ // Returns true if given year is an MDA year for the strategy
        return find(mda_years.begin(), mda_years.end(), Year) != mda_years.end();
    }

};

#endif /* mda_h */