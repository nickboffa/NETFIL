#include "network.h"
#include <stdio.h>
#include <cstdlib>
#include <algorithm>

void Region::implement_mda(int year, MDAStrat strat){
    int n_pop = 0;
    int n_treated = 0;
    int n_under_min = 0;

    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //for every group
        
        Group *grp = j->second;
        
        n_pop += grp->group_pop.size();
        
        for(map<int, Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //for every person 
            
            double age = k->second->age/365.0; // agent's age

            if (age<strat.min_age) ++n_under_min; 

        }
    }

    double target_prop = 1 - n_under_min /(double)n_pop;

    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //for every group
        
        Group *grp = j->second;
        
        for(map<int, Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //for all people 
            
            double age = k->second->age/365.0; // agent's age
            if(age >= strat.min_age){
                if(random_real() <= strat.coverage/(double)target_prop){
                    ++n_treated;
                    k->second->mda(strat.drug);
                }
            }
        }
    }

    number_treated[year] = n_treated;
    achieved_coverage[year] = n_treated/(double)n_pop;
}

void Region::handle_commute(int year){
    
    //firstly need to clear previous storage
   
    if (year % RECALC_YEARS == 0){    
        
        if (groups.size() > 1){
            radt_model(DISTANCE_TYPE); //generating commuting network
            for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //now using the network
                Group *grp = j->second;
                int no_commute_id = grp->gid;
                double commuter_prop = grp->total_commute / (double) grp->group_pop.size();
                //now iterating over all group members
                for(map<int, Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
                    Agent *agt = k->second; //our agent
                    
                    double cum_sum_floor = 0;
                    if(random_real() > commuter_prop){ //Will not commute!
                        grp->day_population.insert(pair<int, Agent*>(agt->aid, agt)); //storing them in current group for day population
                        agt->dgp = groups[no_commute_id];
                    } 
                    else{ //person will commute
                        int commute_id;
                        double commute_dest = random_real();
                        //but commute where?
                        for(map<int, double>::iterator i = grp->commuting_cumsum.begin(); i != grp->commuting_pop.end(); ++i){
                            
                            if ((cum_sum_floor < commute_dest) && (commute_dest <= i->second)){
                                commute_id = i->first;
                                agt->dgp = groups[commute_id]; //assigning agent to day group
                                groups[commute_id]->day_population.insert(pair<int, Agent*>(agt->aid, agt)); //storing them in commuting group for day population
                                goto found_commute;
                            } 
                            else {
                                cum_sum_floor = i->second;
                            }   
                        }
                    }
                    found_commute:;
                }
            }
        }
        //we will also update all agents biting probs 

    }
    
    rpop = 0;
    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        Group *grp = j->second;   
        rpop += grp->group_pop.size();
    }
}

void Region::calc_risk(){
    
    char form = 'l'; //l for limitation, f for facilation, or anything else for linear 
    bool single = false;

    if(groups.size() == 1){
        single = true;
    }

    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //setting the force of transmission in each group to 0
        Group *grp = j->second;
        
        grp->day_strength = 0;
        grp->night_strength = 0;
        grp->night_bites = 0;
        grp->day_bites = 0;
        double db = 0;
        double nb = 0;

        for(map<int, Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            Agent *agt =k->second;
            int age = int(agt->age / 365);
            double c = 1.0;
            if(age <= 15) c = exposure_by_age[age];
            nb += agt->bite_scale*c;    
        }
        if(!single){
            for(map<int, Agent*>::iterator k = grp->day_population.begin(); k != grp->day_population.end(); ++k){
                Agent *agt =k->second;
                int age = int(agt->age / 365);
                double c = 1.0;
                if(age <= 15) c = exposure_by_age[age];
                db += agt->bite_scale*c;    
            }
            grp->day_bites = db;
        }else{
            grp->day_bites = 0;
        }
        grp->night_bites = nb;
       
    }
    //Finding strength of infection in each group
    //now looping over all infected agents
    for(map<int, Agent*>::iterator j = inf_indiv.begin(); j != inf_indiv.end(); ++j){
        
        Agent *agt =j->second;
        Group *ngrp = agt->ngp; //infected agents nightime group

        int age = int(agt->age / 365);
        double c = 1.0;

        if(age <= 15) c = exposure_by_age[age];
        
        ngrp->night_strength += (c*agt->bite_scale*mf_functional_form(form, j->second->worm_strength)) / ngrp->night_bites;
        
        if(!single){
            Group *dgrp = agt->dgp; //infected agents daytime group
            dgrp->day_strength += (c*agt->bite_scale*mf_functional_form(form, j->second->worm_strength)) / dgrp->day_bites;
        }
    }

    //Now finding infective bites
    
    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //looping over groups
        Group *grp = j->second;
        //looping over all people!
        
        for(map<int, Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //looping over all people will do both night and day bites in same loop
        
            Agent *agt = k->second; //our person
            char prev_status = agt ->status;
            double c = 1; //
            int age = int(agt->age / 365);
            if (age <= 15) c = exposure_by_age[age];
            
            agt->sim_bites(c, worktonot, single); // simulating the bites!

            if(agt->status == 'E' && prev_status == 'S'){
                pre_indiv.insert(pair<int, Agent*>(agt->aid, agt));
            }
        }
    }
}

double Region::mf_functional_form(char form, double worm_strength){
    if(form == 'l'){ // limitation
       
        return theta1*theta3*(1-exp(-theta2*worm_strength));
    }
    else if(form == 'f'){ //facilitation
        
        return theta1*(worm_strength - theta2*worm_strength / (1 + theta3*worm_strength));
    }
    else{ //asumme linear
        return theta1*worm_strength;
    }
}

void Region::update_epi_status(int year, int day, int dt){

    for(map<int, Agent*>::iterator j = pre_indiv.begin(); j != pre_indiv.end();){ //looking at all agents with immature worms but not a set!

        Agent *agt = j->second;
        agt->update(year,day,dt);

        if(agt->status != 'E'){ //agent has left this stage of infection

            agt->changed_epi_today = true;
            
            pre_indiv.erase(j++);

            if(agt->status == 'I'){ //infective now!
                inf_indiv.insert(pair<int, Agent*>(agt->aid, agt));
            }
            else if(agt->status == 'U'){ //do not have a set of mature worms!
                uninf_indiv.insert(pair<int, Agent*>(agt->aid, agt));
            }
        }
        else ++j;
    }
           

    for(map<int, Agent*>::iterator j = uninf_indiv.begin(); j != uninf_indiv.end();){ //looking at all agents with onlymature immature worms!
        
        Agent *agt = j->second;

        if(!agt->changed_epi_today) agt->update(year, day,dt);

        else agt->changed_epi_today = false;

        if(agt->status != 'U'){
            uninf_indiv.erase(j++);
            if(agt->status == 'I'){
                agt->changed_epi_today = true;
                inf_indiv.insert(pair<int, Agent*>(agt->aid, agt));
            }
            else if(agt->status == 'E'){
                pre_indiv.insert(pair<int, Agent*>(agt->aid, agt));
            }
        }
        else ++j;
    }

    for(map<int, Agent*>::iterator j = inf_indiv.begin(); j != inf_indiv.end();){
        Agent *agt = j->second;
        if(!agt->changed_epi_today) agt->update(year, day,dt);
        else agt->changed_epi_today = false;
        if(agt->status != 'I'){
            inf_indiv.erase(j++);
            if(agt->status == 'E'){
                pre_indiv.insert(pair<int, Agent*>(agt->aid, agt));
            }
            else if(agt->status == 'U'){
                uninf_indiv.insert(pair<int, Agent*>(agt->aid, agt));
            }
        }
        else ++j;
    }

}

void Region::renew_pop(int year, int day, int dt){
    //handleing deaths!
    vector<Agent*> deaths;

    for(map<int, Group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
        Group *grp = j->second;

        for(map<int, Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //going through group members
            Agent *agt = k->second;
            int index = int(int(agt->age/365)/5);
            if(index > 15) index = 15; //all 75+ the same

            double prob = 1 - exp(-mortality_rate[index]*dt);
            if(random_real() < prob) deaths.push_back(agt); //seeing if agent dies depending on age
            else agt->age += dt; //increase everyones age
        }
    }
    while(deaths.size() > 0){ //now removing agents that have died
        Agent *agt = deaths.back();
        remove_agent(agt);
        deaths.pop_back();
    }       
}

void Region::remove_agent(Agent *agt){

    //removing agent from lists of infected
   
    if(agt->status == 'E') pre_indiv.erase(agt->aid);
    else if(agt->status == 'I') inf_indiv.erase(agt->aid);
    else if(agt->status == 'U') uninf_indiv.erase(agt->aid);
    
    
    //nightime group    
    Group *ngrp = agt->ngp;
    ngrp->group_pop.erase(agt->aid);
    
    if (groups.size() >1 ){
        //daytime group
        Group *dgrp = agt->dgp;
        dgrp->day_population.erase(agt->aid);
    }
   
    
    delete agt;
}

void Region::handle_birth(int year, int day, int dt){ //deal with births
    
    int total_births  = 0;

    for(map<int,Group*>::iterator j = groups.begin(); j != groups.end(); j++){//looping over groups
        Group *grp = j->second;
        for(map<int, Agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){//over agents
           
            Agent *agt = k->second;
            if(agt->age >= 15*365 && agt->age < 50*365){
                int index = int((int(agt->age/365))/5);
                double prob = 1 - exp(-birth_rate[index]*dt);
                
                if(random_real() < prob) ++total_births; 
            }
            
        }
        //now assigning births
        while (total_births > 0) {
            Agent *bby = new Agent(next_aid++, agg_param, 0); //have birth!
            bby->ngp = grp;
            bby->dgp = grp;
            grp->add_member(bby); //assigning baby to group
            grp->day_population.insert(pair<int, Agent*>(bby->aid, bby));//assume baby stays within group during day
            
            --total_births;
        }
    }
}