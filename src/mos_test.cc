
#include "G4RunManager.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "G4UImanager.hh"
#include "PhysicsList.hh"

#ifdef G4VIS_USE

# include "G4VisExecutive.hh"

#endif // ifdef G4VIS_USE

#ifdef G4UI_USE

# include "G4UIExecutive.hh"

#endif // ifdef G4UI_USE

#include <chrono>
#include "G4PhysListFactory.hh"
#include "myUtils.hh"

double INITIAL_ENERGY_MEV = 0;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv) {

    const char *env_var[2] = {"CC", "CXX"};
    char *env_val[2];
    for (int i = 0; i < 2; ++i) {
        env_val[i] = getenv(env_var[i]);
        if (env_val[i] != NULL)
            G4cout << "Variable = " << env_var[i] << ", Value= " << env_val[i] << G4endl;
        else {
            G4cout << env_var[i] << " doesn't exist" << G4endl;
        }
    }

    const double WT1 = myUtils::get_wall_time();

    G4cout << "Total available RAM : " << myUtils::check_available_RAM() << " GB" << G4endl;

    Settings *settings = Settings::getInstance();

    setlocale(LC_ALL, "C");  // just in case, to avoid potential bug for PARMA ("," <-> ".")

    G4String NB_PARTICLES_TO_SHOOT = "10000";

    if (argc == 4) {
        settings->POTENTIAL_VALUE = std::stod(argv[1]);
        NB_PARTICLES_TO_SHOOT = argv[2];
        settings->INITIAL_SAMPLE_TYPE = std::stoi(argv[3]);
    } else if (argc > 4) {
        G4cout << "ERROR. Executable should have 0 or 3 input arguments." << G4endl;
        std::abort();
    } else if (argc > 1 && argc < 4) {
        G4cout << "ERROR. Executable should have 0 or 3 input arguments." << G4endl;
        std::abort();
    }

    if (settings->POTENTIAL_VALUE != 0.) {
        settings->initial_efield_status = settings->efield_ON;
        settings->current_efield_status = settings->efield_ON;
    } else if (settings->POTENTIAL_VALUE == 0.) {
        settings->initial_efield_status = settings->efield_OFF;
        settings->current_efield_status = settings->efield_OFF;
        settings->USE_STACKING_ACTION = false;
    }

    settings->RANDOM_SEED = myUtils::generate_a_unique_ID();

    G4cout << "First random seed for GEANT4: " << settings->RANDOM_SEED << G4endl;

    if (settings->EFIELD_REGION_Y_CENTER > 20.0) {
        G4cout << "ERROR. Efield should be set with a center below 20 km" << G4endl;
        std::abort();
    }

    // choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::MixMaxRng);
    long seeds[2];
    seeds[0] = settings->RANDOM_SEED;
    seeds[1] = myUtils::generate_a_unique_ID();
    CLHEP::HepRandom::setTheSeeds(seeds, 2);

    AnalysisManager *analysis = AnalysisManager::getInstance();

    const double f_rrea_t = analysis->get_fraction_RREA_thres(settings->POTENTIAL_VALUE, settings->EFIELD_REGION_Y_CENTER, settings->EFIELD_REGION_Y_LEN);

//    if (f_rrea_t < 0.3) settings->USE_STACKING_ACTION = false;

    auto *runManager = new G4RunManager;

    //    myExceptionHandler *ExceptionHandler = new myExceptionHandler(0); // G4ExceptionHandler class with verbosity option (0 -> silent)

    runManager->SetVerboseLevel(0);

    // set mandatory initialization classes
    auto *det = new DetectorConstruction;
    runManager->SetUserInitialization(det);

    G4bool EM_ONLY = false;
    G4VUserPhysicsList *phys = new PhysicsList(EM_ONLY);
    phys->SetVerboseLevel(0);
    runManager->SetUserInitialization(phys);

    runManager->SetUserInitialization(new ActionInitialization(det));

    // get the pointer to the User Interface manager
    G4UImanager *UI = G4UImanager::GetUIpointer();
    UI->SetVerboseLevel(0);

    if (settings->MODE == "run") {
        // Initialize G4 kernel
        runManager->Initialize();

        //     UI->ApplyCommand("/run/printProgress 1000000");
        //            AnalysisManager *analysis = AnalysisManager::getInstance();

        UI->ApplyCommand("/run/beamOn " + NB_PARTICLES_TO_SHOOT);

        analysis->write_output_files();

//        return 0;

    } else if (settings->MODE == "visu") // define visualization and UI terminal for interactive mode
    {
        runManager->Initialize();
#ifdef G4VIS_USE
        G4VisManager *visManager = new G4VisExecutive;
        visManager->SetVerboseLevel(0);
        visManager->Initialize();
#endif // ifdef G4VIS_USE
#ifdef G4UI_USE
        G4UIExecutive *ui_ex = new G4UIExecutive(argc, argv);
        UI->ApplyCommand("/control/execute vis.mac");
        UI->SetVerboseLevel(0);
        ui_ex->SessionStart();
        delete ui_ex;
#endif // ifdef G4UI_USE
#ifdef G4VIS_USE
        delete visManager;
#endif // ifdef G4VIS_USE
    } else {
        G4cout << G4endl << "ERROR : Mode should be set to 'visu' or 'run'. Aborting. " << G4endl;
        std::abort();
    }

    // job termination

    delete runManager;

    G4double WT2 = myUtils::get_wall_time();

    G4cout << "GEANT4 run finished" << G4endl;
    G4cout << "Potential: " << settings->POTENTIAL_VALUE << G4endl;
    G4cout << "Time taken: " << (WT2 - WT1) / 1.0e6 << " seconds" << G4endl;

    ///
    int out = system("mkdir -p ./logs/");
    std::ofstream asciiFile_log;
    asciiFile_log.open("./logs/log_" + std::to_string(settings->RANDOM_SEED) + "_" + std::to_string(int(settings->POTENTIAL_VALUE)) + ".txt",
                       std::ios::out | std::ios::app);

    double duration_s = (WT2 - WT1) / 1.0e6;
    asciiFile_log << "Run ID (rng seed): " << settings->RANDOM_SEED << G4endl;
    asciiFile_log << "Time taken: " << duration_s << " seconds" << G4endl;
    asciiFile_log << "Number of events: " << settings->NB_EVENT << G4endl;
    asciiFile_log << "Number of outputs: " << analysis->NB_OUTPUT << G4endl;

    asciiFile_log << "Potential: " << settings->POTENTIAL_VALUE << G4endl;
    asciiFile_log << "E-field mid altitude: " << settings->EFIELD_REGION_Y_CENTER << G4endl;
    asciiFile_log << "E-field length: " << settings->EFIELD_REGION_Y_LEN << G4endl;

    asciiFile_log << "Record altitude: " << settings->RECORD_ALTITUDE << G4endl;

    asciiFile_log << "Average time per 10 events: " << duration_s / double(settings->NB_EVENT) * 10.0 << " seconds" << G4endl;
    asciiFile_log << "Average time per 10 records: " << duration_s / double(analysis->NB_OUTPUT) * 10.0 << " seconds" << G4endl;

    asciiFile_log << "Has reached CPU time per event limit: " << settings->CPU_TIME_LIMIT_PER_EVENT_HAS_BEEN_REACHED_ONCE << G4endl;
//    asciiFile_log << "Has reached RAM limit: " << settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED << G4endl;
    asciiFile_log.close();

    return 0;
}
