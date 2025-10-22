#ifndef RANDOM_GEN_H
#define RANDOM_GEN_H
#include <random>
using namespace std;

extern mt19937 gen; // Declaration of the generator

// Reseed the global RNG so repeated runs are independent
void reseed_rng(unsigned int seed);

#endif // RANDOM_GEN_H