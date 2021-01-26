// ********************************************************************

// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************

#include "AnalysisManager.hh"

// class following singleton pattern

AnalysisManager *AnalysisManager::instance = nullptr;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager *AnalysisManager::getInstance() // singleton lazy initialization
{
    if (instance == nullptr) {
        instance = new AnalysisManager;
    }

    return instance;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
AnalysisManager::AnalysisManager() // constructor
{
    // initialization of spectrums in energy and momentum direction
    memset(PART_SPEC, 0, sizeof(PART_SPEC));
    memset(PART_MOM_X, 0, sizeof(PART_MOM_X));
    memset(PART_MOM_Y, 0, sizeof(PART_MOM_Y));
    memset(PART_MOM_Z, 0, sizeof(PART_MOM_Z));

    memset(counter_up, 0, sizeof(counter_up));
    memset(counter_down, 0, sizeof(counter_down));

    memset(counter_total, 0, sizeof(counter_total));

    // sanity checks


    if (ENER_GRID.size() != size_grid_ener) {
        G4cout << "Error in AnalysisManager.hh/cc : energy grid vector does not have 256 elements." << G4endl;
        std::abort();
    }

    if (MOM_GRID.size() != size_grid_mom) {
        G4cout << "Error in AnalysisManager.hh/cc : momentum-direction grid vector does not have 128 elements." << G4endl;
        std::abort();
    }

    //////

    RECORDED_OUTPUT_STRINGS.clear();

    G4String second_part;
    G4String first_part;

    G4String pot = std::to_string(int(settings->POTENTIAL_VALUE)) + "MV";

    G4String dir_pot = "mkdir -p ./output/" + pot;
    G4String dir_seed = "/" + std::to_string(settings->RANDOM_SEED) + "/";
    system(dir_pot);

    if (settings->INITIAL_SAMPLE_TYPE == 22) {
        G4String mkdir_str = "mkdir -p ./output/" + pot + "/initial_photon";
        system(mkdir_str);
        G4String mkdir_str2 = "mkdir -p ./output/" + pot + "/initial_photon" + dir_seed;
        system(mkdir_str2);
        first_part = "./output/" + pot + "/initial_photon/" + dir_seed;
    } else if (settings->INITIAL_SAMPLE_TYPE == 11) {
        G4String mkdir_str = "mkdir -p ./output/" + pot + "/initial_electron";
        system(mkdir_str);
        G4String mkdir_str2 = "mkdir -p ./output/" + pot + "/initial_electron" + dir_seed;
        system(mkdir_str2);
        first_part = "./output/" + pot + "/initial_electron/" + dir_seed;
    } else if (settings->INITIAL_SAMPLE_TYPE == -11) {
        G4String mkdir_str = "mkdir -p ./output/" + pot + "/initial_positron";
        system(mkdir_str);
        G4String mkdir_str2 = "mkdir -p ./output/" + pot + "/initial_positron" + dir_seed;
        system(mkdir_str2);
        first_part = "./output/" + pot + "/initial_positron/" + dir_seed;
    } else if (settings->INITIAL_SAMPLE_TYPE == 13) {
        G4String mkdir_str = "mkdir -p ./output/" + pot + "/initial_muonN";
        system(mkdir_str);
        G4String mkdir_str2 = "mkdir -p ./output/" + pot + "/initial_muonN" + dir_seed;
        system(mkdir_str2);
        first_part = "./output/" + pot + "/initial_muonN/" + dir_seed;
    } else if (settings->INITIAL_SAMPLE_TYPE == -13) {
        G4String mkdir_str = "mkdir -p ./output/" + pot + "/initial_muonP";
        system(mkdir_str);
        G4String mkdir_str2 = "mkdir -p ./output/" + pot + "/initial_muonP" + dir_seed;
        system(mkdir_str2);
        first_part = "./output/" + pot + "/initial_muonP/" + dir_seed;
    } else if (settings->INITIAL_SAMPLE_TYPE == 2112) {
        G4String mkdir_str = "mkdir -p ./output/" + pot + "/initial_neutron";
        system(mkdir_str);
        G4String mkdir_str2 = "mkdir -p ./output/" + pot + "/initial_neutron" + dir_seed;
        system(mkdir_str2);
        first_part = "./output/" + pot + "/initial_neutron/" + dir_seed;
    } else if (settings->INITIAL_SAMPLE_TYPE == 2212) {
        G4String mkdir_str = "mkdir -p ./output/" + pot + "/initial_proton";
        system(mkdir_str);
        G4String mkdir_str2 = "mkdir -p ./output/" + pot + "/initial_proton" + dir_seed;
        system(mkdir_str2);
        first_part = "./output/" + pot + "/initial_proton/" + dir_seed;
    }

    second_part = first_part + "photon_ener_mom_dists_";

    G4String potential_str = "_" + std::to_string(int(settings->POTENTIAL_VALUE));
    G4String efield_alt_str = "_" + std::to_string(int(settings->EFIELD_REGION_Y_CENTER * 10.));
    G4String efield_length_str = "_" + std::to_string(int(settings->EFIELD_REGION_Y_LEN * 100.));
    G4String rec_alt_part = "_" + std::to_string(int(settings->RECORD_ALTITUDE * 10.));
    asciiFileName_phot = second_part + std::to_string(settings->RANDOM_SEED) + rec_alt_part + potential_str + efield_alt_str + efield_length_str + ".out";

    second_part = first_part + "electron_ener_mom_dists_";
    asciiFileName_elec = second_part + std::to_string(settings->RANDOM_SEED) + rec_alt_part + potential_str + efield_alt_str + efield_length_str + ".out";

    second_part = first_part + "positron_ener_mom_dists_";
    asciiFileName_posi = second_part + std::to_string(settings->RANDOM_SEED) + rec_alt_part + potential_str + efield_alt_str + efield_length_str + ".out";

    second_part = first_part + "muonn_ener_mom_dists_";
    asciiFileName_mn = second_part + std::to_string(settings->RANDOM_SEED) + rec_alt_part + potential_str + efield_alt_str + efield_length_str + ".out";

    second_part = first_part + "muonp_ener_mom_dists_";
    asciiFileName_mp = second_part + std::to_string(settings->RANDOM_SEED) + rec_alt_part + potential_str + efield_alt_str + efield_length_str + ".out";

    G4cout << "Reserving output file : " << asciiFileName_phot << G4endl;
    G4cout << "Reserving output file : " << asciiFileName_elec << G4endl;
    G4cout << "Reserving output file : " << asciiFileName_posi << G4endl;
    G4cout << "Reserving output file : " << asciiFileName_mp << G4endl;
    G4cout << "Reserving output file : " << asciiFileName_mn << G4endl;

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager::~AnalysisManager() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::add_NB_OUTPUT() {

    NB_OUTPUT++;

#ifndef NDEBUG // debug mode
    G4cout << "Number of outputs : " << NB_OUTPUT << G4endl;
#endif // ifndef NDEBUG
//    if (NB_OUTPUT % 100 == 0) {
//        G4cout << "Number of outputs : " << NB_OUTPUT << G4endl;
//    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void AnalysisManager::write_output_files() {

    clean_output_files();

    // Photons
    asciiFile_analysis_phot.open(asciiFileName_phot, std::ios::out | std::ios::app);

    if (asciiFile_analysis_phot.is_open()) {
        asciiFile_analysis_phot << std::scientific << std::setprecision(5)
                                << settings->RANDOM_SEED
                                << " " << settings->PDG_LIST[settings->i_p]
                                << " " << settings->NB_EVENT
                                << " " << settings->EFIELD_REGION_Y_CENTER
                                << " " << settings->EFIELD_REGION_Y_LEN
                                << " " << settings->POTENTIAL_VALUE
                                << " " << settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED
                                << " " << settings->NB_EVENTS_WITH_DETECTION[settings->i_p]
                                << " " << settings->NB_LONG_CPU_TIME
                                << " " << settings->INITIAL_SAMPLE_TYPE
                                << " " << settings->CURRENT_WEIGHT;

        asciiFile_analysis_phot << G4endl << G4endl;

        asciiFile_analysis_phot << counter_total[settings->i_p] << " ";

        asciiFile_analysis_phot << G4endl << G4endl;

        for (uint ii = 0; ii < ngride; ++ii) {
            asciiFile_analysis_phot << PART_SPEC[settings->i_p][ii] << " ";
        }

        asciiFile_analysis_phot << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_phot << PART_MOM_X[settings->i_p][ii] << " ";
        }

        asciiFile_analysis_phot << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_phot << PART_MOM_Y[settings->i_p][ii] << " ";
        }

        asciiFile_analysis_phot << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_phot << PART_MOM_Z[settings->i_p][ii] << " ";
        }

        asciiFile_analysis_phot << G4endl;
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting." << G4endl;
        std::abort();
    }

    asciiFile_analysis_phot.close();

    // Electrons
    asciiFile_analysis_elec.open(asciiFileName_elec, std::ios::out | std::ios::app);

    if (asciiFile_analysis_elec.is_open()) {
        asciiFile_analysis_elec << std::scientific << std::setprecision(5)
                                << settings->RANDOM_SEED
                                << " " << settings->PDG_LIST[settings->i_e]
                                << " " << settings->NB_EVENT
                                << " " << settings->EFIELD_REGION_Y_CENTER
                                << " " << settings->EFIELD_REGION_Y_LEN
                                << " " << settings->POTENTIAL_VALUE
                                << " " << settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED
                                << " " << settings->NB_EVENTS_WITH_DETECTION[settings->i_e]
                                << " " << settings->NB_LONG_CPU_TIME
                                << " " << settings->INITIAL_SAMPLE_TYPE
                                << " " << settings->CURRENT_WEIGHT;

        asciiFile_analysis_elec << G4endl << G4endl;

        asciiFile_analysis_elec << counter_total[settings->i_e] << " ";

        asciiFile_analysis_elec << G4endl << G4endl;

        for (uint ii = 0; ii < ngride; ++ii) {
            asciiFile_analysis_elec << PART_SPEC[settings->i_e][ii] << " ";
        }

        asciiFile_analysis_elec << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_elec << PART_MOM_X[settings->i_e][ii] << " ";
        }

        asciiFile_analysis_elec << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_elec << PART_MOM_Y[settings->i_e][ii] << " ";
        }

        asciiFile_analysis_elec << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_elec << PART_MOM_Z[settings->i_e][ii] << " ";
        }

        asciiFile_analysis_elec << G4endl;
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting." << G4endl;
        std::abort();
    }

    asciiFile_analysis_elec.close();

    // Positrons
    asciiFile_analysis_posi.open(asciiFileName_posi, std::ios::out | std::ios::app);

    if (asciiFile_analysis_posi.is_open()) {
        asciiFile_analysis_posi << std::scientific << std::setprecision(5)
                                << settings->RANDOM_SEED
                                << " " << settings->PDG_LIST[settings->i_ep]
                                << " " << settings->NB_EVENT
                                << " " << settings->EFIELD_REGION_Y_CENTER
                                << " " << settings->EFIELD_REGION_Y_LEN
                                << " " << settings->POTENTIAL_VALUE
                                << " " << settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED
                                << " " << settings->NB_EVENTS_WITH_DETECTION[settings->i_ep]
                                << " " << settings->NB_LONG_CPU_TIME
                                << " " << settings->INITIAL_SAMPLE_TYPE
                                << " " << settings->CURRENT_WEIGHT;

        asciiFile_analysis_posi << G4endl << G4endl;

        asciiFile_analysis_posi << counter_total[settings->i_ep] << " ";

        asciiFile_analysis_posi << G4endl << G4endl;

        for (uint ii = 0; ii < ngride; ++ii) {
            asciiFile_analysis_posi << PART_SPEC[settings->i_ep][ii] << " ";
        }

        asciiFile_analysis_posi << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_posi << PART_MOM_X[settings->i_ep][ii] << " ";
        }

        asciiFile_analysis_posi << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_posi << PART_MOM_Y[settings->i_ep][ii] << " ";
        }

        asciiFile_analysis_posi << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_posi << PART_MOM_Z[settings->i_ep][ii] << " ";
        }

        asciiFile_analysis_posi << G4endl;
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting." << G4endl;
        std::abort();
    }

    asciiFile_analysis_posi.close();

    // muons -
    asciiFile_analysis_mn.open(asciiFileName_mn, std::ios::out | std::ios::app);

    if (asciiFile_analysis_mn.is_open()) {
        asciiFile_analysis_mn << std::scientific << std::setprecision(5)
                              << settings->RANDOM_SEED
                              << " " << settings->PDG_LIST[settings->i_mn]
                              << " " << settings->NB_EVENT
                              << " " << settings->EFIELD_REGION_Y_CENTER
                              << " " << settings->EFIELD_REGION_Y_LEN
                              << " " << settings->POTENTIAL_VALUE
                              << " " << settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED
                              << " " << settings->NB_EVENTS_WITH_DETECTION[settings->i_mn]
                              << " " << settings->NB_LONG_CPU_TIME
                              << " " << settings->INITIAL_SAMPLE_TYPE
                              << " " << settings->CURRENT_WEIGHT;

        asciiFile_analysis_mn << G4endl << G4endl;

        asciiFile_analysis_mn << counter_total[settings->i_mn] << " ";

        asciiFile_analysis_mn << G4endl << G4endl;

        for (uint ii = 0; ii < ngride; ++ii) {
            asciiFile_analysis_mn << PART_SPEC[settings->i_mn][ii] << " ";
        }

        asciiFile_analysis_mn << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_mn << PART_MOM_X[settings->i_mn][ii] << " ";
        }

        asciiFile_analysis_mn << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_mn << PART_MOM_Y[settings->i_mn][ii] << " ";
        }

        asciiFile_analysis_mn << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_mn << PART_MOM_Z[settings->i_mn][ii] << " ";
        }


        asciiFile_analysis_mn << G4endl;
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting." << G4endl;
        std::abort();
    }

    asciiFile_analysis_mn.close();

    if (settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED) {
        G4cout << G4endl << "Aborting run since the limit ram usage has been exceeded." << G4endl;
        std::abort();
    }

    // muons +
    asciiFile_analysis_mp.open(asciiFileName_mp, std::ios::out | std::ios::app);

    if (asciiFile_analysis_mp.is_open()) {
        asciiFile_analysis_mp << std::scientific << std::setprecision(5)
                              << settings->RANDOM_SEED
                              << " " << settings->PDG_LIST[settings->i_mp]
                              << " " << settings->NB_EVENT
                              << " " << settings->EFIELD_REGION_Y_CENTER
                              << " " << settings->EFIELD_REGION_Y_LEN
                              << " " << settings->POTENTIAL_VALUE
                              << " " << settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED
                              << " " << settings->NB_EVENTS_WITH_DETECTION[settings->i_mp]
                              << " " << settings->NB_LONG_CPU_TIME
                              << " " << settings->INITIAL_SAMPLE_TYPE
                              << " " << settings->CURRENT_WEIGHT;

        asciiFile_analysis_mp << G4endl << G4endl;

        asciiFile_analysis_mp << counter_total[settings->i_mp] << " ";

        asciiFile_analysis_mp << G4endl << G4endl;

        for (uint ii = 0; ii < ngride; ++ii) {
            asciiFile_analysis_mp << PART_SPEC[settings->i_mp][ii] << " ";
        }

        asciiFile_analysis_mp << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_mp << PART_MOM_X[settings->i_mp][ii] << " ";
        }

        asciiFile_analysis_mp << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_mp << PART_MOM_Y[settings->i_mp][ii] << " ";
        }

        asciiFile_analysis_mp << G4endl << G4endl;

        for (uint ii = 0; ii < ngridm; ++ii) {
            asciiFile_analysis_mp << PART_MOM_Z[settings->i_mp][ii] << " ";
        }

        asciiFile_analysis_mp << G4endl;
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting." << G4endl;
        std::abort();
    }

    asciiFile_analysis_mp.close();

}

// ======================================================================
// Returns interpolated value at x from parallel arrays ( xData, yData )
//   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
//   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
G4double AnalysisManager::interpolate(std::vector<G4double> &xData, std::vector<G4double> &yData, G4double x, bool extrapolate) {
    int size = static_cast<int>(xData.size());

    int i = 0;                // find left end of interval for interpolation

    if (x >= xData[size - 2]) // special case: beyond right end
    {
        i = size - 2;
    } else {
        while (x > xData[i + 1]) i++;
    }

    double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1]; // points on either side (unless beyond ends)

    if (!extrapolate)                                                          // if beyond ends of array and not extrapolating
    {
        if (x < xL) {
            yR = yL;
        }

        if (x > xR) {
            yL = yR;
        }
    }

    double dydx = (yR - yL) / (xR - xL); // gradient

    return yL + dydx * (x - xL);         // linear interpolation
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double AnalysisManager::get_scale(const G4double &alt) {
    // returns the atmospheric relative scale compared to sea level (>1)
    // input is altitude in km
    std::vector<G4double> alt_list =
            {0,
             0.500000000000000,
             1.000000000000000,
             1.500000000000000,
             2.000000000000000,
             2.500000000000000,
             3.000000000000000,
             3.500000000000000,
             4.000000000000000,
             4.500000000000000,
             5.000000000000000,
             5.500000000000000,
             6.000000000000000,
             6.500000000000000,
             7.000000000000000,
             7.500000000000000,
             8.000000000000000,
             8.500000000000000,
             9.000000000000000,
             9.500000000000000,
             10.000000000000000,
             10.500000000000000,
             11.000000000000000,
             11.500000000000000,
             12.000000000000000,
             12.500000000000000,
             13.000000000000000,
             13.500000000000000,
             14.000000000000000,
             14.500000000000000,
             15.000000000000000,
             15.500000000000000,
             16.000000000000000,
             16.500000000000000,
             17.000000000000000,
             17.500000000000000,
             18.000000000000000,
             18.500000000000000,
             19.000000000000000,
             19.500000000000000,
             20.000000000000000,};

    std::vector<G4double> scale_list =
            {1.000000000000000,
             1.059301380991064,
             1.121238177128117,
             1.184377838328792,
             1.249042145593870,
             1.317304778260430,
             1.388415672913118,
             1.463688404983724,
             1.543377914546100,
             1.628371628371629,
             1.719862833025587,
             1.818435364663227,
             1.925291598996014,
             2.041647095663066,
             2.168634624979211,
             2.307556184746062,
             2.460377358490566,
             2.628502318081032,
             2.813376483279396,
             3.018518518518519,
             3.245395719263315,
             3.496916063287745,
             3.774240231548481,
             4.080100125156446,
             4.414353419092755,
             4.781811514484781,
             5.180770758839889,
             5.613430908308223,
             6.082089552238807,
             6.589186457806973,
             7.133479212253830,
             7.720544701006513,
             8.348271446862997,
             9.024221453287199,
             9.753178758414361,
             10.541632983023444,
             11.398601398601400,
             12.325141776937620,
             13.334696799263730,
             14.431164231961047,
             15.627996164908916,};

    if ((alt > 20.) || (alt < 0.)) {
        G4cout << "ERROR in get_scale : altitude of E field is not between 0 and 20. Aborting.";
        std::abort();
    }

    return interpolate(alt_list, scale_list, alt, false);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisManager::fill_histogram_E(const G4int idx_part, const G4double &value) {
    for (uint ii = 0; ii < ngride - 1; ++ii) {
        if ((value >= ENER_GRID[ii]) && (value < ENER_GRID[ii + 1])) {
            PART_SPEC[idx_part][ii]++;
            return;
        }

        if (value >= 100.0) {
            PART_SPEC[idx_part][ngride - 1]++;
            return;
        }
    }
}

void AnalysisManager::fill_histogram_mX(const G4int idx_part, const G4double &value) {
    for (uint ii = 0; ii < ngridm - 1; ++ii) {
        if ((value >= MOM_GRID[ii]) && (value < MOM_GRID[ii + 1])) {
            PART_MOM_X[idx_part][ii]++;
            return;
        }
    }
}

void AnalysisManager::fill_histogram_mY(const G4int idx_part, const G4double &value) {
    for (uint ii = 0; ii < ngridm - 1; ++ii) {
        if ((value >= MOM_GRID[ii]) && (value < MOM_GRID[ii + 1])) {
            PART_MOM_Y[idx_part][ii]++;
            return;
        }
    }
}

void AnalysisManager::fill_histogram_mZ(const G4int idx_part, const G4double &value) {
    for (uint ii = 0; ii < ngridm - 1; ++ii) {
        if ((value >= MOM_GRID[ii]) && (value < MOM_GRID[ii + 1])) {
            PART_MOM_Z[idx_part][ii]++;
            return;
        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// open / close with trunc to clean the files (if it does exist) or to create it
void AnalysisManager::clean_output_files() {
    asciiFile_analysis_phot.open(asciiFileName_phot, std::ios::trunc);
    if (asciiFile_analysis_phot.is_open()) {
        asciiFile_analysis_phot.close();
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting" << G4endl;
        std::abort();
    }

    asciiFile_analysis_elec.open(asciiFileName_elec, std::ios::trunc);
    if (asciiFile_analysis_elec.is_open()) {
        asciiFile_analysis_elec.close();
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting" << G4endl;
        std::abort();
    }

    asciiFile_analysis_posi.open(asciiFileName_posi, std::ios::trunc);
    if (asciiFile_analysis_posi.is_open()) {
        asciiFile_analysis_posi.close();
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting" << G4endl;
        std::abort();
    }

    asciiFile_analysis_mp.open(asciiFileName_mp, std::ios::trunc);
    if (asciiFile_analysis_mp.is_open()) {
        asciiFile_analysis_mp.close();
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting" << G4endl;
        std::abort();
    }

    asciiFile_analysis_mn.open(asciiFileName_mn, std::ios::trunc);
    if (asciiFile_analysis_mn.is_open()) {
        asciiFile_analysis_mn.close();
    } else {
        G4cout << G4endl << "ERROR : cannot open output file. Aborting" << G4endl;
        std::abort();
    }
}

// ------------------------------------------------------------------------
// check to make sure particle not recorded more than 2 times (one upward and one downward)
// returns true if the particle can be recorded, false otherwise
bool AnalysisManager::check_record(const uint &id, const double &momy) {

    if (record_list.empty()) return true;

    bool up = momy > 0.0;
    bool down = !up;

    bool found = false;

    for (auto &a_record : record_list) {
        if (id == a_record.ID) {
            found = true;
            if (a_record.downward && a_record.upward) {
                G4cout << "SKIPPED" << G4endl;
                return false;
            } else if (!a_record.downward && down) {
                a_record.downward = true;
                return true;
            } else if (!a_record.upward && up) {
                a_record.upward = true;
                return true;
            }
        }
    }

    if (!found) {
        record_list.push_back({id, up, down});
        return true;
    }

    return true;

}

// ------------------------------------------------------------------------

double AnalysisManager::get_fraction_RREA_thres(const double potential_MV, const double center_alt_km, const double length_km) {

    const double RREA_thes_see_level = 284.0; // MV/km
    const double scale = get_scale(center_alt_km + length_km / 2.0);

    const double RREA_thres = RREA_thes_see_level / scale;

    const double fraction_RREA_thres = potential_MV / length_km / RREA_thres;

    return fraction_RREA_thres;
}

// ------------------------------------------------------------------------

double AnalysisManager::get_avalanchge_time(const double potential_MV, const double center_alt_km, const double length_km) {

    const double RREA_thes_see_level = 284.0; // MV/km
    const double scale = get_scale(center_alt_km + length_km / 2.0);

    const double RREA_thres = RREA_thes_see_level / scale;

    const double fraction_RREA_thres = potential_MV / length_km / RREA_thres;

    return fraction_RREA_thres;
}