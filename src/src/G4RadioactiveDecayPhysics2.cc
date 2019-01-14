// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RadioactiveDecayPhysics2.hh"

#include "G4GenericIon.hh"
#include "G4RadioactiveDecay.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"

//
G4_DECLARE_PHYSCONSTR_FACTORY(G4RadioactiveDecayPhysics2);

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RadioactiveDecayPhysics2::G4RadioactiveDecayPhysics2(G4int) : G4VPhysicsConstructor("G4RadioactiveDecay") // , theRadioactiveDecay(0)
{}

G4RadioactiveDecayPhysics2::G4RadioactiveDecayPhysics2(const G4String& name) : G4VPhysicsConstructor(name)  // , theRadioactiveDecay(0)
{}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RadioactiveDecayPhysics2::~G4RadioactiveDecayPhysics2()
{
  // delete theRadioactiveDecay;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RadioactiveDecayPhysics2::ConstructParticle()
{
  G4GenericIon::GenericIon();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4RadioactiveDecayPhysics2::ConstructProcess()
{
  G4RadioactiveDecay *radioactiveDecay = new G4RadioactiveDecay();

  radioactiveDecay->SetICM(true); // Internal Conversion
  radioactiveDecay->SetARM(true); // Atomic Rearangement
  G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
