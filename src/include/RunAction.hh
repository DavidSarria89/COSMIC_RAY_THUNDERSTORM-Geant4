#pragma once

#include "G4UserRunAction.hh"
#include "globals.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;

class PrimaryGeneratorAction;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction: public G4UserRunAction
{
public:

  RunAction();

  ~RunAction() override;

public:

  void BeginOfRunAction(const G4Run *) override;

  void EndOfRunAction(const G4Run *) override;

private:
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
