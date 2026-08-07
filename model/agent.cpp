#include "agent.h"
#include "mda.h"
#include "network.h"
#include <limits>

//Constructer of agent
Agent::Agent(int aid,  double bite_shape, int age){
    this->aid = aid;
    this->age = age;
   
    changed_epi_today = false;
    status = 'S';

    worm_strength = 0;

    last_mworm_time = - std::numeric_limits<double>::infinity();

    bite_scale = bite_gamma(bite_shape, 1/bite_shape);
}

Agent::~Agent(){
    dgp = NULL;
    ngp = NULL;
    for(int i = 0; i < wvec.size(); ++i){
        delete wvec[i];
    }
    wvec.clear();

}

void Agent::sim_bites(double c, double worktonot, bool single,
                      double imm_mean, double imm_std,
                      double mat_mean, double mat_std,
                      double prop_male_worm){

    int total_bites;

    if (single) {
        total_bites = poisson(c *  ngp->night_strength * bite_scale);
    } else {
        int day_bites;
        int night_bites;

        day_bites   = poisson(c * dgp->day_strength * bite_scale * worktonot);
        night_bites = poisson(c * ngp->night_strength * bite_scale * (1.0 - worktonot));

        total_bites = day_bites + night_bites;
    }

    for (int i = 0; i < total_bites; ++i) { //looping through infective bites and assigning worms
        int immature_period = normal(imm_mean, imm_std); //immature period of worm
        int mature_period = normal(mat_mean, mat_std); //mature period of worm

        if (random_real() < prop_male_worm) { // worm is male
            wvec.push_back(new Worm('P', immature_period, mature_period, 'M'));
        } else { // worm is female
            wvec.push_back(new Worm('P', immature_period, mature_period, 'F'));
        }
    }

    if(total_bites > 0 && status == 'S') status = 'E';
}

void Agent::mda(Drugs drug){
    if(wvec.size() > 0){ //if person has worms
        double rr = random_real(); //same thing will occur to all worms!
        for(int i = 0; i < wvec.size(); i++){ // looping through worms
            wvec[i]->age_mda = drug.ster_dur*365;
            if (rr <= drug.kill_prob){
                wvec[i]->status = 'D';
            }
            else if(rr <= drug.kill_prob + drug.full_ster_prob){ // sterilise worms with probability full_ster_prob
                wvec[i]->mda_sterile = 0.0; //worm is sterile!
            }
            else if (rr <= drug.kill_prob + drug.full_ster_prob + drug.part_ster_prob){ //Partially sterilise with probability part_ster_prob
                wvec[i]->mda_sterile = min(wvec[i]->mda_sterile, 1 - drug.part_ster_magnitude); // worm is partially sterile, we also ensure that mda does NOT increase infectivity of an already sterilised worm
            }
        }
    }
}

void Agent::update(int year, int day, int dt){
    // First update status of all worms in the body
    if(wvec.size() > 0){
        for(int i = 0; i < wvec.size();){ 
            wvec[i]->update(dt);
            
            // Remove dead worms
            if(wvec[i]->status == 'D'){ 
                delete wvec[i];
                wvec.erase(wvec.begin() + i);
            }
            else ++i; // iterator increment placed here because
            // if a worm is removed from the list the next worm 
            // will now have the index of deleted worm
        }
    }

    char prevstatus = status; // Agent's previous infection status

    bool mature_worm = false; // Has mature worms of either sex (fertile does not matter)

    double worm_strength_female = 0.0; // Counter to track strength of females
    double worm_strength_male = 0.0; // Counter to track strength of males

    // Update agent's epi status
    if (wvec.size() == 0) { // no worms => susceptible/uninfected
        status = 'S';
        worm_strength = 0.0;
    } else {
        for (int i = 0; i < wvec.size(); ++i) { // iterating over worms
            if (wvec[i]->status == 'M') { // Mature
                mature_worm = true;

                if (wvec[i]->sex == 'M') { // Male worm
                    if (wvec[i]->mda_sterile > 0) {
                        worm_strength_male += wvec[i]->mda_sterile;
                    }
                } else if (wvec[i]->sex == 'F') { // Female worm
                    if (wvec[i]->mda_sterile > 0) {
                        worm_strength_female += wvec[i]->mda_sterile;
                    }
                } 
            }
        }
        if((worm_strength_female > 0) && (worm_strength_male > 0)){ //agent is infectious!
            status = 'I'; // person is infectious
            
            if (worm_strength_male > 1.0) worm_strength_male = 1.0; //We assume polygamous worms that relies on females

            worm_strength = worm_strength_male * worm_strength_female; // total worm strength

        } else { 
            worm_strength = 0.0;
            status = (mature_worm) ? 'U' : 'E';
            // U: Agent has worm(s) but not a mature fertile set
            // E: Agent has only immature worm(s)
        }
    }

    // Record if worm has died
    if ((prevstatus == 'U' || prevstatus == 'I') && (status == 'S' || status == 'E')) {
        last_mworm_time = year * 365 + day*dt;
    }
}

// Update worms
void Worm::update(int dt){
    // Three statuses: P: immature, M: mature, D: dead 

    if (status == 'P') {
        if (age_immature <= 0) { 
            status = 'M'; // Mature
        } else {
            age_immature -= dt; // Countdown
        }
    }

    if (status == 'M') {
        if (age_mature <= 0) {
            status = 'D'; // Dead
        } else {
            age_mature -= dt;
        }
    }

    if (age_mda > 0) {
        age_mda -= dt; // Now looking at sterilisation status of worm
    } else {
        mda_sterile = 1.0; // Sterile period now over
    }
}