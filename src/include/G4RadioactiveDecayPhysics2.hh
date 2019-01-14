// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#pragma once

#include "G4VPhysicsConstructor.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4RadioactiveDecay;

class G4RadioactiveDecayPhysics2: public G4VPhysicsConstructor
{
public:

  G4RadioactiveDecayPhysics2(G4int verbose = 0);

  G4RadioactiveDecayPhysics2(const G4String &name);

  virtual ~G4RadioactiveDecayPhysics2();

public:

  // This method is dummy for physics
  virtual void ConstructParticle();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

private:

  // G4RadioactiveDecay*  theRadioactiveDecay;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
