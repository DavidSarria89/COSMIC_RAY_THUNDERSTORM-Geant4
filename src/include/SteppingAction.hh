#pragma once

#include "AnalysisManager.hh"
#include "G4ThreeVector.hh"
#include "G4UserSteppingAction.hh"
#include "RegionInformation.hh"
#include "globals.hh"
#include <vector>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#include "Settings.hh"

#include <sys/time.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include "myUtils.hh"

class DetectorConstruction;

class EventAction;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

struct index_found {
    uint index;
    bool found;
};

class SteppingAction : public G4UserSteppingAction {
public:

    SteppingAction(DetectorConstruction *, EventAction *);

    ~SteppingAction() override;

    void UserSteppingAction(const G4Step *aStep) override;

private:

    G4String get_name(int PDG);

    index_found find_particle_index(const int PDG_in);

    uint ALL_TICKS_COUNTER = 0;

    uint CHECK_COUNTER = 0;
    const uint RAM_CHECK_COUNTER_MAX = 3 * 500000;
    const double TIME_PRINT_SEC = 0.5 * 60.0;

    double WT1 = 0;
    double WT2 = 0;

    double WT3 = 0;
    double WT4 = 0;

    double INITIAL_TIME = 0;

    Settings *settings = Settings::getInstance();

    AnalysisManager *analysis = AnalysisManager::getInstance();

    DetectorConstruction *fDetector = nullptr;
    EventAction *fEventAction = nullptr;

    G4StepPoint *thePrePoint = nullptr;
    G4StepPoint *thePostPoint = nullptr;
    G4Track *theTrack = nullptr;


    bool is_inside_eField_region(const G4double &alt,
                                 const G4double &xx,
                                 const G4double &zz);

    const double EFIELD_alt_min = settings->EFIELD_REGION_Y_CENTER - settings->EFIELD_REGION_Y_LEN / 2.0; // km
    const double EFIELD_alt_max = settings->EFIELD_REGION_Y_CENTER + settings->EFIELD_REGION_Y_LEN / 2.0; // km

    std::vector<int> PDG_LST = settings->PDG_LIST;

    bool first_time = false;

};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
