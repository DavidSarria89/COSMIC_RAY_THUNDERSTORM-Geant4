#include "cosmic_ray_generator_PARMA.hh"

extern "C" {
void gen_parma_cr_(const int *,
                   double[],       // MeV
                   double[],       // cosine of zenith angle (e.g. ang=1.0 for vertical direction, ang=0.0 for holizontal direction)
                   double[],       // km
                   int[],
                   const int *,    // Particle ID (Particle ID, 31:e-, 32:e+, 33:photon)
                   const int *,    // size of energy mesh (will be log)
                   const int *,    // size of angle mesh (linear)
                   const int *,    // size of altitude mesh (linear)
                   const double *, // minimum altitude (km)
                   const double *, // maximum altitude (km)
                   const double *, // emin MeV
                   const double *, // emax MeV
                   const int *,    // iyear
                   const int *,    // imonth
                   const int *,    // iday
                   const double *, // glat deg -90 =< glat =< 90
                   const double *, // glong deg -180 =< glat =< 180
                   int[],          // wanted Particle ID list (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
                   const int *);   // number of wanted Particle ID list

// See custom_subroutines.f90 for more info

void get_parma_particles_weights_(double[],
                                  int[],
                                  const int *,
                                  const double *,
                                  const double *,
                                  const double *,
                                  const double *,
                                  const int *,    // iyear
                                  const int *,    // imonth
                                  const int *,    // iday
                                  const double *, // glat deg -90 =< glat =< 90
                                  const double *);   // number of wanted Particle ID list
}

Cosmic_Ray_Generator_PARMA::Cosmic_Ray_Generator_PARMA() {
    setlocale(LC_ALL, "C"); // just in case

    // set up the list of indexes of particles we want

    int parmaID = find_parma_ID_from_PDG(settings->INITIAL_SAMPLE_TYPE);
    parmaID_wanted = parmaID;
    particles_type_wanted = G4ParticleTable::GetParticleTable()->FindParticle(settings->INITIAL_SAMPLE_TYPE);

    // first call to PARMA to generate the list of cosmic rays
    Generate_CR_samples_list_from_PARMA();

    settings->WEIGHTS.clear();
    settings->WEIGHTS = Calculate_Weights_from_PARMA();

    G4cout << "\nWeights: " << G4endl;

    int ii = 0;
    for (const double &val : settings->WEIGHTS) {
        G4cout << "  " << settings->NAMES[ii] << ": " << val << G4endl;
        ii++;
    }
    G4cout << G4endl;

    if (settings->INITIAL_SAMPLE_TYPE == 22) settings->CURRENT_WEIGHT = settings->WEIGHTS[i_phot];
    else if (settings->INITIAL_SAMPLE_TYPE == 11) settings->CURRENT_WEIGHT = settings->WEIGHTS[i_elec];
    else if (settings->INITIAL_SAMPLE_TYPE == -11) settings->CURRENT_WEIGHT = settings->WEIGHTS[i_posi];
    else if (settings->INITIAL_SAMPLE_TYPE == 13) settings->CURRENT_WEIGHT = settings->WEIGHTS[i_muN];
    else if (settings->INITIAL_SAMPLE_TYPE == -13) settings->CURRENT_WEIGHT = settings->WEIGHTS[i_muP];
    else if (settings->INITIAL_SAMPLE_TYPE == 2112) settings->CURRENT_WEIGHT = settings->WEIGHTS[i_neut];
    else if (settings->INITIAL_SAMPLE_TYPE == 2212) settings->CURRENT_WEIGHT = settings->WEIGHTS[i_prot];
    else std::abort();

#ifndef NDEBUG
    if (settings->WRITE_MOM_OUTPUT_FOR_TEST) {
        std::ofstream asciiFile7;
        asciiFile7.open(name_outFile_mom, std::ios::trunc);
        asciiFile7.close();
    }

    if (settings->CR_GENRATOR_write_output_FOR_TEST) {
        generate_output_for_test();
    }
#endif
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int Cosmic_Ray_Generator_PARMA::find_parma_ID_from_PDG(const int &PDG_in) {

    const int pdg_phot = 22;
    const int pdg_elec = 11;
    const int pdg_posi = -11;
    const int pdg_muN = 13;
    const int pdg_muP = -13;
    const int pdg_neut = 2112;
    const int pdg_prot = 2212;

    const int parma_phot = 33;
    const int parma_elec = 31;
    const int parma_posi = 32;
    const int parma_muN = 30;
    const int parma_muP = 29;
    const int parma_neut = 0;
    const int parma_prot = 1;

    switch (PDG_in) {
        case pdg_phot:
            return parma_phot;

        case pdg_elec:
            return parma_elec;

        case pdg_posi:
            return parma_posi;

        case pdg_muN:
            return parma_muN;

        case pdg_muP:
            return parma_muP;

        case pdg_neut:
            return parma_neut;

        case pdg_prot:
            return parma_prot;

        default:
            std::cout << "Error: not a valid PDG number in find_parma_ID_from_PDG in cosmic_ray_generator_PARMA.cc" << std::endl;
            std::abort();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

geant4_initial_cosmic_ray Cosmic_Ray_Generator_PARMA::generate_Cosmic_ray() {

    geant4_initial_cosmic_ray g4_cosmic = {};

    //            sampled_set = read_particles[index_sampling_part];
    cosmic_ray_parma_output sampled_set = sample_One_CR_from_PARMA();

    G4ThreeVector position_ini = sample_CR_secondary_position(sampled_set.altitude);

    // from PARMA OUTPUT:
    //   cos(theta)=1,  indicates  the  vertical  downward  direction,
    //   while  90  degree,  i.e.  cos(theta)=0, indicates the horizontal direction.
    G4ThreeVector momentum_ini = CR_direction_rand_sample(-1.0 * sampled_set.cos_zenith_angle);
    // multiplication by -1 is important, to make sure that when sampled_set.cos_zenith_angle is 1, the particle is sampled vertical downward

    const double time = 0.0;
    const double energy = sampled_set.energy * MeV;

    g4_cosmic.energy = energy;
    g4_cosmic.time = time;
    g4_cosmic.momentum_ini = momentum_ini;
    g4_cosmic.position_ini = position_ini;
    g4_cosmic.g4_particle = particles_type_wanted;

#ifndef NDEBUG
    if (settings->WRITE_MOM_OUTPUT_FOR_TEST) {
        generate_output_for_test_momentum(momentum_ini);
    }
#endif
    return g4_cosmic;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Following angular distributions compuer from PARMA code
G4ThreeVector Cosmic_Ray_Generator_PARMA::CR_direction_rand_sample(const double &cos_Sampled) {
    G4ThreeVector momentum_ini{};

    // if cos_Sampled == 1 (zenith) then direction should be (0,1,0)

    double uu = cos_Sampled;
    double theta = rand_double() * 2.0 * CLHEP::pi;

    momentum_ini.setX(sqrt(1.0 - uu * uu) * cos(theta));
    momentum_ini.setY(uu);
    momentum_ini.setZ(sqrt(1.0 - uu * uu) * sin(theta));

    return momentum_ini;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector Cosmic_Ray_Generator_PARMA::sample_CR_secondary_position(const double &altitude) {
    double r1 = rand_double();
    double r2 = rand_double();

    G4ThreeVector position = {(r1 - 0.5) * settings->CR_SAMPLING_XZ_HALF_SIZE * 2.0 * CLHEP::km,
                              altitude * CLHEP::km,
                              (r2 - 0.5) * settings->CR_SAMPLING_XZ_HALF_SIZE * 2.0 * CLHEP::km};

    return position;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// uses cumulative distribution sampling
cosmic_ray_parma_output Cosmic_Ray_Generator_PARMA::sample_One_CR_from_PARMA() {

    if (counter >= nb_int - 3) // if all the particles generated from Parma have been already used, generate new ones
    {
        Generate_CR_samples_list_from_PARMA();

#ifndef NDEBUG
        if (settings->CR_GENRATOR_write_output_FOR_TEST) {
            generate_output_for_test();
        }
#endif
        G4cout << "Generated " << nb_int << " new random cosmic ray particles." << G4endl;
    }

    const double eRand = output_energies[counter];
    const double cos_angRand = output_cosangles[counter];
    const double alt_Rand = output_altitudes[counter];
    counter++;

#ifndef NDEBUG // if debug mode, some sanity checks

    if ((alt_Rand < min_alt) || (alt_Rand > max_alt)) {
        G4cout << "Sampled altitude is not in the fixed altitude range. Aborting." << G4endl;
        std::abort();
    }

    if ((cos_angRand < min_cosAng) || (cos_angRand > max_cosAng)) {
        G4cout << "Sampled cosine of angle is not between -1 and 1. Aborting." << G4endl;
        std::abort();
    }

    if ((eRand < min_cr_ener) || (eRand > max_cr_ener)) {
        G4cout << "Energy is out of range. Aborting." << G4endl;
        std::abort();
    }

    if (idx_particle_sampled > nb_part_type_wanted) {
        G4cout << "ERROR : sampled type is not in sampled ID list. Aborting." << G4endl;
        std::abort();
    }

#endif // ifndef NDEBUG

    cosmic_ray_parma_output spld_set{cos_angRand, alt_Rand, eRand};

    return spld_set;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Cosmic_Ray_Generator_PARMA::Generate_CR_samples_list_from_PARMA() {

    int nb = nb_int;

    int *parmaID_list_wanted2 = &parmaID_wanted;

    const int the_seed = myUtils::generate_a_unique_ID_int32();

    gen_parma_cr_(&the_seed,
                  output_energies,
                  output_cosangles,
                  output_altitudes,
                  output_types,
                  &nb,
                  &nebin,
                  &nabin,
                  &naltbin,
                  &min_alt,
                  &max_alt,
                  &min_cr_ener,
                  &max_cr_ener,
                  &iyear,
                  &imonth,
                  &iday,
                  &glat,
                  &glong,
                  parmaID_list_wanted2,
                  &nb_part_type_wanted);

    counter = 0;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<double> Cosmic_Ray_Generator_PARMA::Calculate_Weights_from_PARMA() {

    double output_w[n_types_for_weights];

    int *parmaID_list_wanted2 = &parmaID_list_ALL[0]; // trick to transform a C++ vector into a C array

    get_parma_particles_weights_(output_w,
                                 parmaID_list_wanted2,
                                 &n_types_for_weights,
                                 &min_alt,
                                 &max_alt,
                                 &min_cr_ener,
                                 &max_cr_ener,
                                 &iyear,
                                 &imonth,
                                 &iday,
                                 &glat,
                                 &glong);

    std::vector<double> output;
    output.clear();
    output.reserve(7);
    for (int ii = 0; ii < n_types_for_weights; ++ii) {
        output.push_back(output_w[ii]);
    }

    return output;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Cosmic_Ray_Generator_PARMA::generate_output_for_test() {
    std::ofstream asciiFile6;
    asciiFile6.open(filename_cr_sampling_test + std::to_string(seed_cr_smpl) + ".txt", std::ios::trunc);

    for (int ii = 0; ii < nb_int; ++ii) {
        asciiFile6 << output_types[ii] << " " << output_energies[ii] << " " << output_altitudes[ii] << " " << output_cosangles[ii] << G4endl;
    }

    asciiFile6.close();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Cosmic_Ray_Generator_PARMA::generate_output_for_test_momentum(const G4ThreeVector &mom) {
    std::ofstream asciiFile7;
    asciiFile7.open(name_outFile_mom, std::ios::app);
    asciiFile7 << mom[0] << " " << mom[1] << " " << mom[2] << G4endl;
    asciiFile7.close();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double Cosmic_Ray_Generator_PARMA::rand_double() {
    return G4UniformRand();
}