#pragma once

#include <chrono>
#include <vector>
#include <string>
#include <G4ios.hh>
#include <iostream>
#include <fstream>
#include <thread>
#include "Settings.hh"

#include "Randomize.hh"

#include "sys/types.h"
#include "sys/sysinfo.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include <uuid/uuid.h>

namespace myUtils {

    long generate_a_unique_ID();

    int generate_a_unique_ID_int32();

    std::vector<double> logspace(const double a, const double b, const int n);

    double get_wall_time();

    double check_available_RAM();

    double check_USED_RAM();

    int parseLine(char *line);

}