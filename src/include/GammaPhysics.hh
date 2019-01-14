// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#pragma once

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class GammaPhysics: public G4VPhysicsConstructor
{
public:

  explicit GammaPhysics(const G4String& name = "gamma");

  ~GammaPhysics() override;

public:

  void ConstructParticle() override
  {}

  void ConstructProcess() override;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
