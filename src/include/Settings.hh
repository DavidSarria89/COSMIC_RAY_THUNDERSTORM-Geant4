#pragma once

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <vector>
#include "G4ParticleDefinition.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
struct geant4_initial_cosmic_ray {
    G4ThreeVector momentum_ini;
    G4ThreeVector position_ini;
    double time;
    double energy;
    G4ParticleDefinition *g4_particle;
};
typedef unsigned int uint;

class Settings {
private:

    Settings() = default; // Private so that it can not be called

    Settings(Settings const &) {}

    // copy constructor is private
    // assignment operator is private
    static Settings *instance;

public:

    static Settings *getInstance();

public:

    G4String MODE = "run";
//    G4String MODE="visu";

    double POTENTIAL_VALUE = 100.0;   // MV, can be overwritten by input argument

    uint NB_PARTICLES_TO_GET = 50;

    const double MAX_POSSIBLE_TIME = 10.0 * second;

    const double RECORD_ALTITUDE = 4.3; // km
//    const double RECORD_ALTITUDE = 0.141; // km

    // Y axis is direction of altitude (positive towards increasing altitude)
    const double EFIELD_REGION_Y_LEN = 2.0;   // km
    // const double EFIELD_REGION_Y_LEN = 0.75;
    const double EFIELD_REGION_Y_CENTER = 5.3;   // km
//     const double EFIELD_REGION_Y_CENTER = 1.125;
    const double EFIELD_XZ_HALF_SIZE = 18.0; // km
    const double CR_SAMPLING_XZ_HALF_SIZE = 20.0; // km
    const double CR_GENERATION_ALT_MIN = 6.4; // km, used only by PARMA
    const double CR_GENERATION_ALT_MAX = 6.5; // km
    const double RECORD_XZ_HALF_SIZE = 8.0; // km
    const double SOIL_ALT_MAX = 4.299; // km
//    double SOIL_ALT_MAX = 0.140;

    const double RAM_FRACTION_LIMIT = 0.5; // between 0.0 and 1.0

    const double drOverR = 0.1;

    enum CR_GEN_ENGINE : int {
        PARMA, CRY
    };

    const CR_GEN_ENGINE chosen_CR_GEN_ENGINE = PARMA;

    enum PDG_nb : int {
        pdg_phot = 22,
        pdg_elec = 11,
        pdg_posi = -11,
        pdg_muP = -13,
        pdg_muN = 13,
        pdg_neut = 2112,
        pdg_prot = 2212
    };

    const double LOW_ENERGY_THRES = 8.0 * keV; // electrons with lower energies do not contribute to RREA

    enum Efield_State {
        efield_ON, efield_OFF
    };

    // list of PDG number of particles that we want to generate from parma and be recorded
    const std::vector<int> PDG_LIST = {pdg_phot, pdg_elec, pdg_posi, pdg_muP, pdg_muN};

    const std::vector<int> PDG_LIST_ALL = {pdg_phot, pdg_elec, pdg_posi, pdg_muP, pdg_muN, pdg_neut, pdg_prot};

    const double WORLD_MAX_ALT = 10;     // km

    const double GLOBAL_MAX_STEP = 5.0 * meter;
    const bool USE_GLOBAL_MAX_STEP = true;
    const double STEP_MAX_VAL_RECORD_REGION = 5.0 * meter;
    const bool USE_RECORD_REGION_STEP_MAX = false;

    const bool CR_GENRATOR_write_output_FOR_TEST = false; // test or not the CR particle generator (i.e. output the list of particles)
    const bool ATMOS_LAYERS_OUTPUT_TO_FILE = false; // write info about atmosphere layers on a file (for test / debug)
    const bool WRITE_MOM_OUTPUT_FOR_TEST = false; // write info about momentum (for test / debug)

    bool USE_STACKING_ACTION = true; // will mimic time oriented simulation if set to true
    // can be a bad idea to set it on, because it can use a lot of memory compared to default G4 behaviour
    // always turned off if potential is 0

    const double ENERGY_MIN_RECORD = 50. * keV;
    const double ENERGY_MAX_RECORD = 100. * MeV;

    ulong RANDOM_SEED = 250; // just for initialization, will be replace at beginning of main

    Efield_State current_efield_status = efield_ON;
    Efield_State initial_efield_status = efield_ON;

    //    const double longitude = 130.5; // ILDAS positron event coordinates
    //    (Australia)
    //    const double latitude = -13.5;  // ILDAS positron event coordinates
    //    (Australia)
    //    const double longitude = -103.5; // deg, FEGS glow coordinates (Colorado)
    //    const double latitude = 39.5;   // deg, FEGS glow coordinates (Colorado)

    const double longitude = 91.172; // Tibet
    const double latitude = 29.65;   // Tibet

//    const double longitude = 136.7; // Kanazawa
    //   const double latitude = 36.55;   // Kanazawa

    const double CR_GENERATION_ENER_MIN = 0.04;   // MeV
    const double CR_GENERATION_ENER_MAX = 0.90e6; // MeV
    const int year = 2016;
    const int month = 1;
    const int day = 20;


    uint RAM_USAGE_LIMIT_HAS_BEEN_REACHED = 0; // a flag indication if the limit of particles for RREA has been reached
    // meaning that real multiplciation factor is higher than the one obtained
    // and maybe that the feedback process is dominating

    std::vector<bool> event_lead_to_detection = {false, false, false, false, false};
    std::vector<uint> NB_EVENTS_WITH_DETECTION = {0, 0, 0, 0, 0};

    ulong NB_EVENT = 0; // initialisation is important here
    double VARIABLE_TIME_LIMIT_GLOBAL = 0;

    double T0 = 0.0;

    bool CPU_TIME_LIMIT_PER_EVENT_HAS_BEEN_REACHED = false;
    bool CPU_TIME_LIMIT_PER_EVENT_HAS_BEEN_REACHED_ONCE = false;
    uint NB_LONG_CPU_TIME = 0;

    const uint MAX_CPU_TICKS = 100;
    uint EVENT_DURATION_in_CPU_TICKS = 0;
};