#pragma once

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>

class DetectorConstruction;

class G4ParticleDefinition;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run: public G4Run
{
public:

  Run(DetectorConstruction *);

  ~Run();

public:

  void EndOfRun();

private:

private:

  DetectorConstruction *fDetector;
  G4ParticleDefinition *fParticle;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
