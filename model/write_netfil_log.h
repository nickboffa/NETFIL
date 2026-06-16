#include <ctime>
#include "network.h"

void write_netfil(
    const string& filename,
    time_t start_time,
    time_t end_time,
    Region *rgn,
    string mda_data
);