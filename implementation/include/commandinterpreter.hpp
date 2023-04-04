#pragma once

#include <array>
#include <cstdint>

#include "solve/definitions.hpp"

bool findSeed(std::array<uint_fast32_t, SEED_LENGTH>& seed, char* argv[], int& i);
void interpretCommandLine(const int argc, char* argv[]);

