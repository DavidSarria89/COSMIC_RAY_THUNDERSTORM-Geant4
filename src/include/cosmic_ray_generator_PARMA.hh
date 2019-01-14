#pragma once

#include "AnalysisManager.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "fortran.hh"
#include "globals.hh"
#include <locale.h>
#include <numeric>
#include <stdlib.h>
#include "Settings.hh"
#include <random>
#include "myUtils.hh"
#include "Randomize.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct cosmic_ray_parma_output {
    int type;
    double cos_zenith_angle;
    double altitude;
    double energy;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Cosmic_Ray_Generator_PARMA {
public:

    Cosmic_Ray_Generator_PARMA();

    ~Cosmic_Ray_Generator_PARMA();

    geant4_initial_cosmic_ray generate_Cosmic_ray();

private:

    double rand_double();

    Settings *settings = Settings::getInstance();

    enum Parma_ID : int {
        parma_phot = 33,
        parma_elec = 31,
        parma_posi = 32,
        parma_muP = 29,
        parma_muN = 30,
        parma_neut = 0,
        parma_prot = 1
    };

    ////////////////////////////////////////////////

    const std::vector<int> PDG_LIST = settings->PDG_LIST_ALL;      // the list of particles types we want

    int nb_part_type_wanted = 0;                                  // initialization, will have the size of ID_list_wanted and particles_wanted,
    // determined at runtime
    std::vector<int> parmaID_list_wanted;                                     // the particles (parma ID) we want
    std::vector<G4ParticleDefinition *> particles_type_wanted; // the particles types we want

    double min_cr_ener = settings->CR_GENERATION_ENER_MIN;      // MeV
    double max_cr_ener = settings->CR_GENERATION_ENER_MAX;      // MeV
    double min_alt = settings->CR_GENERATION_ALT_MIN;       // km
    double max_alt = settings->CR_GENERATION_ALT_MAX;       // km
    double min_cosAng = -1.0;
    double max_cosAng = 1.0;

    int iyear = settings->year;
    int imonth = settings->month;
    int iday = settings->day;
    double glat = settings->latitude;
    double glong = settings->longitude;

    int nebin = 512; // size of energy mesh (log)
    int nabin = 128; // size of angle mesh (linear)
    int naltbin = 2;   // size of altitude mesh (linear)
    // the bigger the numbers, more precise will be the sampling, but the more memory it will take

    int counter = 0;      // indexing of the sampled cosmic rays
    int seed_cr_smpl = 123456; // dummy value, random seed of CR generator

    const static int nb_int = 1000000; // number of cosmic to generate. If all are used, it will generate again
    double output_energies[nb_int];
    double output_cosangles[nb_int];
    double output_altitudes[nb_int];
    int output_types[nb_int];

    // Functions

    G4ThreeVector CR_direction_rand_sample(const double &cos_Sampled);

    G4ThreeVector sample_CR_secondary_position(const double &altitude);

    cosmic_ray_parma_output sample_One_CR_from_PARMA();

    void generate_output_for_test();

    void generate_output_for_test_momentum(const G4ThreeVector &mom);

    void Generate_CR_samples_list_from_PARMA();

    Parma_ID find_parma_ID_from_PDG(const int PDG_in);

    ///
    G4String name_outFile_mom = "./tests/cr_sampl_test_mom.txt";
    G4String filename_cr_sampling_test = "./tests/cr_sampl_test_";
};
