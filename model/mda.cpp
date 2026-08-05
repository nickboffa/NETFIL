#include <iostream>
#include <fstream>
#include "mda.h"
#include "json.hpp"

using namespace std;
using json = nlohmann::json;

CountryConfig load_country_config(const string& path) {
    ifstream f(path);
    if (!f) {
        cerr << "Could not open country config: " << path << "\n";
        exit(1);
    }
    json j = json::parse(f);

    CountryConfig cfg;
    cfg.init_year = j.at("init_year");
    cfg.end_year  = j.at("end_year");

    const auto& s = j.at("seeding");
    cfg.init_ant_prev     = s.at("init_ant_prev");
    cfg.init_ant_prev_min = s.at("init_ant_prev_min");
    cfg.init_ant_prev_max = s.at("init_ant_prev_max");

    if (s.contains("init_mf_prev_min")) {
        cfg.use_mf_prev_mode  = true;
        cfg.init_mf_prev_min  = s.at("init_mf_prev_min");
        cfg.init_mf_prev_max  = s.at("init_mf_prev_max");
    } else {
        cfg.use_mf_prev_mode       = false;
        cfg.init_mf_ant_ratio_min  = s.at("init_mf_ant_ratio_min");
        cfg.init_mf_ant_ratio_max  = s.at("init_mf_ant_ratio_max");
    }

    for (const auto& r : j.value("mda_rounds", json::array())) {
        MDARound round;
        round.year     = r.at("year");
        round.month    = r.value("month", 0);
        round.drug     = r.at("drug");
        round.coverage = r.at("coverage");
        round.min_age  = r.value("min_age", 2);
        cfg.mda_rounds.push_back(round);
    }

    cout << "Loaded country config: "
         << cfg.init_year << "–" << cfg.end_year
         << " (" << cfg.sim_years() << " years, "
         << cfg.mda_rounds.size() << " MDA rounds)\n";

    return cfg;
}

MDAEvent make_mda_event(const CountryConfig& cfg, int calendar_year) {
    for (const auto& r : cfg.mda_rounds) {
        if (r.year == calendar_year) {
            MDAEvent evt;
            evt.apply      = true;
            evt.day        = (r.month > 0) ? (r.month - 1) * 30 : 28;
            evt.drug_name  = r.drug;
            evt.drug_params = DRUG_DICT.at(r.drug);
            evt.coverage   = r.coverage;
            evt.min_age    = r.min_age;
            return evt;
        }
    }
    return MDAEvent{};
}
