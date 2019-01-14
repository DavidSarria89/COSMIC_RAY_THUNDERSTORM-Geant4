#pragma once

#include "AnalysisManager.hh"
#include "G4UserEventAction.hh"
#include "Settings.hh"
#include "globals.hh"

#include <sys/time.h>
#include <time.h>

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction {
public:

    EventAction();

    ~EventAction() override;

public:

    void BeginOfEventAction(const G4Event *) override;

    void EndOfEventAction(const G4Event *) override;

private:

    Settings *settings = Settings::getInstance();

    ulong print_nb = 1; // initialisation

    AnalysisManager *analysis = AnalysisManager::getInstance();

};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
