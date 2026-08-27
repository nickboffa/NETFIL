#include <iostream>
#include <fstream>
#include "network.h"
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

void Region::load_config() {
    string default_path = string(DATADIR) + "default_params.json";
    ifstream df(default_path);
    if (!df) {
        cerr << "Could not open default params: " << default_path << "\n";
        exit(1);
    }
    json j = json::parse(df);

    string country_path = string(DATADIR) + COUNTRY_CONFIG;
    ifstream f(country_path);
    if (!f) {
        cerr << "Could not open country config: " << country_path << "\n";
        exit(1);
    }
    j.merge_patch(json::parse(f));  // country.json wins on any shared key

    init_year = j.at("init_year");
    end_year  = j.at("end_year");
    sim_years = end_year - init_year + 1;  // end_year is inclusive: the last calendar year simulated
    end_day   = j.value("end_day", 363);

    const auto& s = j.at("seeding");
    init_ant_prev     = s.at("init_ant_prev");
    init_ant_prev_min = s.at("init_ant_prev_min");
    init_ant_prev_max = s.at("init_ant_prev_max");

    if (s.contains("init_mf_prev_min")) {
        use_mf_prev_mode  = true;
        init_mf_prev      = s.at("init_mf_prev");
        init_mf_prev_min  = s.at("init_mf_prev_min");
        init_mf_prev_max  = s.at("init_mf_prev_max");
    } else {
        use_mf_prev_mode       = false;
        init_mf_ant_ratio_min  = s.at("init_mf_ant_ratio_min");
        init_mf_ant_ratio_max  = s.at("init_mf_ant_ratio_max");
    }

    for (const auto& r : j.value("mda_rounds", json::array())) {
        MDARound round;
        round.year     = r.at("year");
        round.month    = r.value("month", 0);
        round.drug     = r.at("drug");
        round.coverage = r.at("coverage");
        round.min_age  = r.at("min_age");
        mda_rounds.push_back(round);
    }

    const auto& ip = j.at("init_params");
    init_beta_b           = ip.at("init_beta_b");
    init_poisson          = ip.at("init_poisson");
    immature_to_antigen   = ip.at("immature_to_antigen");
    immature_and_ant      = ip.at("immature_and_ant");

    proportion_male_worm     = j.at("proportion_male_worm");
    proportion_male_agent    = j.at("proportion_male_agent");
    daily_prob_lose_ant      = j.at("daily_prob_lose_ant");
    immature_period_mean     = j.at("immature_period_mean");
    immature_period_mean_std = j.at("immature_period_mean_std");
    mature_period_mean       = j.at("mature_period_mean");
    mature_period_mean_std   = j.at("mature_period_mean_std");
    commuting_prop           = j.at("commuting_prop");

    auto read_array16 = [&](const char* key, array<double, 16>& arr) {
        const auto& v = j.at(key);
        for (int i = 0; i < 16; ++i) arr[i] = v[i];
    };
    read_array16("mortality_rates",  mortality_rates);
    read_array16("birth_rates",      birth_rates);
    read_array16("exposure_kernel",  exposure_kernel);

    const auto& ad = j.at("age_dist");
    for (int i = 0; i < 17; ++i) age_dist[i] = ad[i];

    cout << "Loaded config: "
         << init_year << "–" << end_year
         << " (" << sim_years << " years, "
         << mda_rounds.size() << " MDA rounds)\n";
}

MDAEvent Region::make_mda_event(int calendar_year) {
    for (const auto& r : mda_rounds) {
        if (r.year == calendar_year) {
            MDAEvent evt;
            evt.apply      = true;
            evt.day        = (r.month > 0) ? (r.month - 1) * 30 : 28;
            evt.drug_name  = r.drug;
            evt.drug_params = lookup_drug(r.drug, ster_to_kill, mda_total_effect);
            evt.coverage   = r.coverage * p_mda;  // p_mda discounts nominal coverage down to effective coverage
            evt.min_age    = r.min_age;
            return evt;
        }
    }
    return MDAEvent{};
}
