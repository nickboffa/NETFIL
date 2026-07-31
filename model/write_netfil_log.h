#include <ctime>
#include "network.h"

void write_netfil(
    const string& filename,
    time_t start_time,
    time_t end_time,
    clock_t cpu_start,
    clock_t cpu_end,
    Region *rgn,
    string mda_data
);