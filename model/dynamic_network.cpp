#include "network.h"
#include "params.h"   // for USE_NONCOMPLIANCE, NONADH_* etc.
#include "agent.h"    // for agent::is_antigen_positive()
#include "rng.h"
#include <stdio.h>
#include <cstdlib>
#include <algorithm>
#include <unordered_set>

void region::implement_MDA(int year, mda_strat strat){
    int n_pop = 0;
    int n_treated = 0;
    int n_under_min = 0;

    // --- Pass 1: count population and under-min-age individuals ---
    for (auto j = groups.begin(); j != groups.end(); ++j) {
        group *grp = j->second;
        n_pop += grp->group_pop.size();

        for (auto k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k) {
            double age = k->second->age / 365.0;
            if (age < strat.min_age) ++n_under_min;
        }
    }

    double target_prop = 1 - n_under_min / (double)n_pop;

    // --- Pass 2: apply MDA with non-compliance ---
    for (auto j = groups.begin(); j != groups.end(); ++j) {
        group *grp = j->second;

        for (auto k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k) {
            agent *a = k->second;
            double age = a->age / 365.0;
            if (age <= strat.min_age) continue;

            // Coverage selection (programmatic target)
            if (random_real() > strat.Coverage / (double)target_prop) continue;

#if USE_NONCOMPLIANCE
            // --- Non-compliance module (persistent or per-round) ---
            bool agpos = a->is_antigen_positive(); // helper in agent.h
            double pr_nonadh = agpos ? NONADH_AGPOS : NONADH_AGNEG;
            bool refuse_this_round = false;

    #if PERSISTENT_NONCOMPLIANCE
            if (!a->mda_profile_init) {
                a->mda_never_taker = (random_real() < pr_nonadh);
                a->mda_profile_init = true;
            }
            refuse_this_round = a->mda_never_taker;
    #else
            refuse_this_round = (random_real() < pr_nonadh);
    #endif

            if (refuse_this_round) continue;
#endif

            // --- If reached, agent is treated ---
            ++n_treated;
            a->mda(strat.drug);
        }
    }

    number_treated[year] = n_treated;
    achieved_coverage[year] = n_treated / (double)n_pop;

    std::cout << "[MDA] Year " << (year + start_year)
              << " treated=" << n_treated
              << " achieved=" << std::fixed << std::setprecision(3)
              << achieved_coverage[year]
              << (PERSISTENT_NONCOMPLIANCE ? " (persistent)\n" : " (per-round)\n");
}

void region::implement_MDA_subset(int year, mda_strat strat, const std::vector<int>& target_gids) {
    // If no list given, fallback to full-pop MDA
    if (target_gids.empty()) return implement_MDA(year, strat);

    std::unordered_set<int> allow(target_gids.begin(), target_gids.end());
    int n_pop = 0, n_under_min = 0, n_treated = 0;

    // Count eligible in target groups only
    for (auto &jk : groups) {
        if (!allow.count(jk.first)) continue;
        group* grp = jk.second;
        n_pop += (int)grp->group_pop.size();
        for (auto &kv : grp->group_pop) {
            double age = kv.second->age / 365.0;
            if (age < strat.min_age) ++n_under_min;
        }
    }
    if (n_pop == 0) return;

    double target_prop = 1.0 - (n_under_min / (double)n_pop);

    // Treat only inside target groups
    for (auto &jk : groups) {
        if (!allow.count(jk.first)) continue;
        group* grp = jk.second;

        for (auto &kv : grp->group_pop) {
            agent *a = kv.second;
            double age = a->age / 365.0;
            if (age <= strat.min_age) continue;

            if (random_real() > strat.Coverage / (double)target_prop) continue;

        #if USE_NONCOMPLIANCE
            bool agpos = a->is_antigen_positive();
            double pr_nonadh = agpos ? NONADH_AGPOS : NONADH_AGNEG;
            bool refuse = false;
            #if PERSISTENT_NONCOMPLIANCE
                if (!a->mda_profile_init) { a->mda_never_taker = (random_real() < pr_nonadh); a->mda_profile_init = true; }
                refuse = a->mda_never_taker;
            #else
                refuse = (random_real() < pr_nonadh);
            #endif
            if (refuse) continue;
        #endif

            ++n_treated;
            a->mda(strat.drug);
        }
    }

    number_treated[year] += n_treated;        // accumulate within-year
    achieved_coverage[year] += n_treated / (double)rpop; // coarse region-level fraction
}

void region::handl_commute(int year){
    
    int recalc_years = 100; //how often we want to recalc commuters
    char distance_type = 'r'; // r for road distance, e for euclidean 

    //firstly need to clear previous storage
   
    if (year % recalc_years == 0){    
        
        if (groups.size() > 1){
            radt_model(distance_type); //generating commuting network
            for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //now using the network
                group *grp = j->second;
                int no_commute_id = grp->gid;
                double commuter_prop = grp->total_commute / (double) grp->group_pop.size();
                //now iterating over all group members
                for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
                    agent *cur = k->second; //our agent
                    
                    double cum_sum_floor = 0;
                    if(random_real() > commuter_prop){ //Will not commute!
                        grp->day_population.insert(pair<int, agent*>(cur->aid, cur)); //storing them in current group for day population
                        cur->dgp = groups[no_commute_id];
                    } 
                    else{ //person will commute
                        int commute_id;
                        double commute_dest = random_real();
                        //but commute where?
                        for(map<int, double>::iterator i = grp->commuting_cumsum.begin(); i != grp->commuting_pop.end(); ++i){
                            
                            if ((cum_sum_floor < commute_dest) && (commute_dest <= i->second)){
                                commute_id = i->first;
                                cur->dgp = groups[commute_id]; //assigning agent to day group
                                groups[commute_id]->day_population.insert(pair<int, agent*>(cur->aid, cur)); //storing them in commuting group for day population
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
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){
        group *grp = j->second;   
        rpop += grp->group_pop.size();
    }
}

void region::calc_risk(){
    
    char form = 'l'; //l for limitation, f for facilation, or anything else for linear 
    bool single = false;

    if(groups.size() == 1){
        single = true;
    }

    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //setting the force of transmission in each group to 0
        group *grp = j->second;
        
        grp->day_strength = 0;
        grp->night_strength = 0;
        grp->night_bites = 0;
        grp->day_bites = 0;
        double db = 0;
        double nb = 0;

        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){
            agent *cur =k->second;
            int age = int(cur->age / 365);
            double c = 1.0;
            if(age <= 15) c = exposure_by_age[age];
            nb += cur->bite_scale*c;    
        }
        if(!single){
            for(map<int, agent*>::iterator k = grp->day_population.begin(); k != grp->day_population.end(); ++k){
                agent *cur =k->second;
                int age = int(cur->age / 365);
                double c = 1.0;
                if(age <= 15) c = exposure_by_age[age];
                db += cur->bite_scale*c;    
            }
            grp->day_bites = db;
        }else{
            grp->day_bites = 0;
        }
        grp->night_bites = nb;
       
    }
    //Finding strength of infection in each group
    //now looping over all infected agents
    for(map<int, agent*>::iterator j = inf_indiv.begin(); j != inf_indiv.end(); ++j){
        
        agent *cur =j->second;
        group *ngrp = cur->ngp; //infected agents nightime group

        int age = int(cur->age / 365);
        double c = 1.0;

        if(age <= 15) c = exposure_by_age[age];
        
        ngrp->night_strength += (c*cur->bite_scale*mf_functional_form(form, j->second->worm_strength)) / ngrp->night_bites;
        
        if(!single){
            group *dgrp = cur->dgp; //infected agents daytime group
            dgrp->day_strength += (c*cur->bite_scale*mf_functional_form(form, j->second->worm_strength)) / dgrp->day_bites;
        }
    }

    //Now finding infective bites
    
    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //looping over groups
        group *grp = j->second;
        //looping over all people!
        
        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //looping over all people will do both night and day bites in same loop
        
            agent *cur = k->second; //our person
            char prev_status = cur ->status;
            double c = 1; //
            int age = int(cur->age / 365);
            if (age <= 15) c = exposure_by_age[age];
            
            cur->sim_bites(c, worktonot, single); // simulating the bites!

            if(cur->status == 'E' && prev_status == 'S'){
                pre_indiv.insert(pair<int, agent *>(cur->aid, cur));
            }
        }
    }
}

double region::mf_prev_region() const {
    // proportion infectious at time of call; convert to %
    if (rpop <= 0) return 0.0;
    return 100.0 * (double)inf_indiv.size() / (double)rpop;
}

std::map<int,double> region::mf_prev_by_group() const {
    std::map<int,double> out;
    for (auto &gkv : groups) {
        const group* g = gkv.second;
        double n = (double)g->group_pop.size();
        if (n <= 0) { out[gkv.first] = 0.0; continue; }
        // Count infectious in group
        int infc = 0;
        for (auto &kv : g->group_pop) {
            const agent* a = kv.second;
            if (a->status == 'I') ++infc;
        }
        out[gkv.first] = 100.0 * infc / n;
    }
    return out;
}

double region::tas_estimate_antigen_15plus(const std::vector<int>& target_gids,
                                           int sample_size,
                                           int current_year) const {
    std::unordered_set<int> allow(target_gids.begin(), target_gids.end());
    std::vector<const agent*> frame;
    frame.reserve(4000);

    // Build sampling frame: adults 15+ in chosen areas
    for (auto &gkv : groups) {
        if (!allow.empty() && !allow.count(gkv.first)) continue;
        const group* g = gkv.second;
        for (auto &kv : g->group_pop) {
            const agent* a = kv.second;
            int years = a->age / 365;
            if (years >= 15) frame.push_back(a);
        }
    }
    if (frame.empty()) return 0.0;

    // Sample without replacement
    std::vector<int> idx(frame.size());
    std::iota(idx.begin(), idx.end(), 0);
    shuffle(idx.begin(), idx.end(), gen); // 'gen' from rng.h / rng.cpp
    int n = std::min<int>(sample_size, (int)frame.size());

    // Antigen-positive proxy: I or U OR lingering antigen (DailyProbLoseAntigen)
    int pos = 0;
    for (int i = 0; i < n; ++i) {
        const agent* a = frame[idx[i]];
        bool agpos = (a->status == 'I' || a->status == 'U');
        if (!agpos) {
            double days_since_last_mature = current_year * 365.0 - a->lastwormtime;
            if (days_since_last_mature > 0) {
                agpos = (random_real() < pow(DailyProbLoseAntigen, days_since_last_mature));
            }
        }
        if (agpos) ++pos;
    }

    return 100.0 * (double)pos / (double)n;
}

double region::mf_functional_form(char form, double worm_strength){
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

void region::update_epi_status(int year, int day, int dt){

    for(map<int, agent*>::iterator j = pre_indiv.begin(); j != pre_indiv.end();){ //looking at all agents with immature worms but not a set!

        agent *p = j->second;
        p->update(year,day,dt);

        if(p->status != 'E'){ //agent has left this stage of infection

            p->ChangedEpiToday = true;
            
            pre_indiv.erase(j++);

            if(p->status == 'I'){ //infective now!
                inf_indiv.insert(pair<int, agent*>(p->aid, p));
            }
            else if(p->status == 'U'){ //do not have a set of mature worms!
                uninf_indiv.insert(pair<int, agent*>(p->aid, p));
            }
        }
        else ++j;
    }
           

    for(map<int, agent*>::iterator j = uninf_indiv.begin(); j != uninf_indiv.end();){ //looking at all agents with onlymature immature worms!
        
        agent *p = j->second;

        if(!p->ChangedEpiToday) p->update(year, day,dt);

        else p->ChangedEpiToday = false;

        if(p->status != 'U'){
            uninf_indiv.erase(j++);
            if(p->status == 'I'){
                p->ChangedEpiToday = true;
                inf_indiv.insert(pair<int, agent*>(p->aid, p));
            }
            else if(p->status == 'E'){
                pre_indiv.insert(pair<int, agent*>(p->aid, p));
            }
        }
        else ++j;
    }

    for(map<int, agent*>::iterator j = inf_indiv.begin(); j != inf_indiv.end();){
        agent *p = j->second;
        if(!p->ChangedEpiToday) p->update(year, day,dt);
        else p->ChangedEpiToday = false;
        if(p->status != 'I'){
            inf_indiv.erase(j++);
            if(p->status == 'E'){
                pre_indiv.insert(pair<int, agent*>(p->aid, p));
            }
            else if(p->status == 'U'){
                uninf_indiv.insert(pair<int, agent*>(p->aid, p));
            }
        }
        else ++j;
    }

}

void region::renew_pop(int year, int day, int dt){
    //handleing deaths!
    vector<agent*> deaths;

    for(map<int, group*>::iterator j = groups.begin(); j != groups.end(); ++j){ //going through groups
        group *grp = j->second;

        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){ //going through group members
            agent *cur = k->second;
            int index = int(int(cur->age/365)/5);
            if(index > 15) index = 15; //all 75+ the same

            double prob = 1 - exp(-mortality_rate[index]*dt);
            if(random_real() < prob) deaths.push_back(cur); //seeing if agent dies depending on age
            else cur->age += dt; //increase everyones age
        }
    }
    while(deaths.size() > 0){ //now removing agents that have died
        agent *cur = deaths.back();
        remove_agent(cur);
        deaths.pop_back();
    }       
}

void region::remove_agent(agent *p){

    //removing agent from lists of infected
   
    if(p->status == 'E') pre_indiv.erase(p->aid);
    else if(p->status == 'I') inf_indiv.erase(p->aid);
    else if(p->status == 'U') uninf_indiv.erase(p->aid);
    
    
    //nightime group    
    group *ngrp = p->ngp;
    ngrp->group_pop.erase(p->aid);
    
    if (groups.size() >1 ){
        //daytime group
        group *dgrp = p->dgp;
        dgrp->day_population.erase(p->aid);
    }
   
    
    delete p;
}

void region::hndl_birth(int year, int day, int dt){ //deal with births
    
    int total_births  = 0;

    for(map<int,group*>::iterator j = groups.begin(); j != groups.end(); j++){//looping over groups
        group *grp = j->second;
        for(map<int, agent*>::iterator k = grp->group_pop.begin(); k != grp->group_pop.end(); ++k){//over agents
           
            agent *cur = k->second;
            if(cur->age >= 15*365 && cur->age < 50*365){
                int index = int((int(cur->age/365))/5);
                double prob = 1 - exp(-birth_rate[index]*dt);
                
                if(random_real() < prob) ++total_births; 
            }
            
        }
        //now assigning births
        while (total_births > 0) {
            agent *bb = new agent(next_aid++, agg_param, 0); //have birth!
            bb->ngp = grp;
            bb->dgp = grp;
            grp->add_member(bb); //assigning baby to group
            grp->day_population.insert(pair<int, agent*>(bb->aid, bb));//assume baby stays within group during day
            
            --total_births;
        }
    }
}