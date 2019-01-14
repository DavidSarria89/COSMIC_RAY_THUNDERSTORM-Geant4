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
}

Cosmic_Ray_Generator_PARMA::Cosmic_Ray_Generator_PARMA() {
    setlocale(LC_ALL, "C"); // just in case

    if (settings->WRITE_MOM_OUTPUT_FOR_TEST) {
        std::ofstream asciiFile7;
        asciiFile7.open(name_outFile_mom, std::ios::trunc);
        asciiFile7.close();
    }

    nb_part_type_wanted = PDG_LIST.size();

    // set up the list of indexes of particles we want

    for (int ii = 0; ii < nb_part_type_wanted; ++ii) {
        Parma_ID parmaID = find_parma_ID_from_PDG(PDG_LIST[ii]);
        parmaID_list_wanted.push_back(static_cast<int>(parmaID));
        particles_type_wanted.push_back(G4ParticleTable::GetParticleTable()->FindParticle(PDG_LIST[ii]));
    }

    // first call to PARMA to generate the list of cosmic rays
    Generate_CR_samples_list_from_PARMA();

    if (settings->CR_GENRATOR_write_output_FOR_TEST) {
        generate_output_for_test();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Cosmic_Ray_Generator_PARMA::Parma_ID Cosmic_Ray_Generator_PARMA::find_parma_ID_from_PDG(const int PDG_in) {
    switch (PDG_in) {
        case Settings::pdg_phot:
            return parma_phot;

        case Settings::pdg_elec:
            return parma_elec;

        case Settings::pdg_posi:
            return parma_posi;

        case Settings::pdg_muP:
            return parma_muP;

        case Settings::pdg_muN:
            return parma_muN;

        case Settings::pdg_neut:
            return parma_neut;

        case Settings::pdg_prot:
            return parma_prot;

        default:
            std::cout << "Error: not a valid PDG number in find_parma_ID_from_PDG in cosmic_ray_generator_PARMA.cc" << std::endl;
            std::abort();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

geant4_initial_cosmic_ray Cosmic_Ray_Generator_PARMA::generate_Cosmic_ray() {
    geant4_initial_cosmic_ray g4_cosmic = {};

    G4ThreeVector position_ini(0., 0., 0.);
    G4ThreeVector momentum_ini(0., 0., 0.);

    G4ParticleDefinition *particle = nullptr;

    cosmic_ray_parma_output sampled_set{0, 0, 0, 0};

    //            sampled_set = read_particles[index_sampling_part];
    sampled_set = sample_One_CR_from_PARMA();

    position_ini = sample_CR_secondary_position(sampled_set.altitude);

    // from PARMA OUTPUT:
    //   cos(theta)=1,  indicates  the  vertical  downward  direction,
    //   while  90  degree,  i.e.  cos(theta)=0, indicates the horizontal direction.
    momentum_ini = CR_direction_rand_sample(-1.0 * sampled_set.cos_zenith_angle);
    // multiplication by -1 is important, to make sure that when sampled_set.cos_zenith_angle is 1, the particle is sampled vertical downward
    //            G4cout << momentum_ini[1] << G4endl;

    if (settings->WRITE_MOM_OUTPUT_FOR_TEST) {
        generate_output_for_test_momentum(momentum_ini);
    }

    const double time = 0.0;
    const double energy = sampled_set.energy * MeV;
    particle = particles_type_wanted[sampled_set.type - 1];

    g4_cosmic.energy = energy;
    g4_cosmic.time = time;
    g4_cosmic.momentum_ini = momentum_ini;
    g4_cosmic.position_ini = position_ini;
    g4_cosmic.g4_particle = particle;

    return g4_cosmic;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Following angular distributions compuer from PARMA code
G4ThreeVector Cosmic_Ray_Generator_PARMA::CR_direction_rand_sample(const double &cos_Sampled) {
    G4ThreeVector momentum_ini;

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
                              (r2 - 0.5) * settings->CR_SAMPLING_XZ_HALF_SIZE * 2.0 * CLHEP::km    };

    return position;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// uses cumulative distribution sampling
cosmic_ray_parma_output Cosmic_Ray_Generator_PARMA::sample_One_CR_from_PARMA() {

    if (counter >= nb_int - 3) // if all the particles generated from Parma have been already used, generate new ones
    {
        Generate_CR_samples_list_from_PARMA();

        if (settings->CR_GENRATOR_write_output_FOR_TEST) {
            generate_output_for_test();
        }

        G4cout << "Generated " << nb_int << " new random cosmic ray particles." << G4endl;
    }

    const double eRand = output_energies[counter];
    const double cos_angRand = output_cosangles[counter];
    const double alt_Rand = output_altitudes[counter];
    const int idx_particle_sampled = output_types[counter];
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

    cosmic_ray_parma_output spld_set{idx_particle_sampled, cos_angRand, alt_Rand, eRand};

    // equivalently:
    //    spld_set.energy = eRand;
    //    spld_set.cos_zenith_angle = cos_angRand;
    //    spld_set.altitude = alt_Rand;
    //    spld_set.type = idx_particle_sampled;

    return spld_set;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Cosmic_Ray_Generator_PARMA::Generate_CR_samples_list_from_PARMA() {


    int nb = nb_int;

    int *parmaID_list_wanted2 = &parmaID_list_wanted[0]; // trick to transform a C++ vector into a C array

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