#include "G4SystemOfUnits.hh"
#include "FieldSetup.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FieldSetup::FieldSetup()
        : fFieldManager(0),
          fLocalFieldManager(0),
          fChordFinder(0),
          fLocalChordFinder(0),
          fEquation(0),
          fElectricField(0),
          fLocalElectricField(0),
          fStepper(0),
          fLocalStepper(0) {

    EFIELD_mag = settings->POTENTIAL_VALUE * megavolt / (settings->EFIELD_REGION_Y_LEN * km);

    fElectricField =
            new G4UniformElectricField_custom(G4ThreeVector(0.0,
                                                     -1.0 * EFIELD_mag,
                                                     0.0));

    fEquation = new G4EqMagElectricField(fElectricField);

    fStepper = new G4ClassicalRK4(fEquation, nvar);

    fFieldManager = GetGlobalFieldManager();

    fFieldManager->SetDetectorField(fElectricField);

    UpdateField();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FieldSetup::~FieldSetup() {
    delete fElectricField;
    delete fChordFinder;
    delete fStepper;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FieldSetup::UpdateField() {
    // It must be possible to call 'again' - e.g. to choose an alternative stepper
    //   has been chosen, or in case other changes have been made.


    // 1. First clean up previous state.
    delete fChordFinder;
    fChordFinder = nullptr;

    // 2. Create the steppers ( Note: this also deletes the previous ones. )
    SetStepper();

    double minEps = 1.0e-6; //   Minimum & value for smallest steps
    double maxEps = 1.0e-5; //   Maximum & value for largest steps

    // 3. Create the chord finder(s)
    EIntgrDriver = new G4MagInt_Driver(fMinStep,
                                       fStepper,
                                       fStepper->GetNumberOfVariables());

    fChordFinder = new G4ChordFinder(EIntgrDriver);

    fFieldManager->SetChordFinder(fChordFinder);

    fFieldManager->SetMinimumEpsilonStep(minEps);
    fFieldManager->SetMaximumEpsilonStep(maxEps);
    fFieldManager->SetDeltaOneStep(10.0 * cm);
    fFieldManager->SetDeltaIntersection(1.0 * cm);

    // 4. Ensure that the field is updated (in Field manager & equation)
    fFieldManager->SetDetectorField(fElectricField);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FieldSetup::SetStepper() {
    delete fStepper;
    fStepper = nullptr;

    fStepper = new G4ClassicalRK4(fEquation, nvar);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4FieldManager *FieldSetup::GetGlobalFieldManager() {
    return G4TransportationManager::GetTransportationManager()
            ->GetFieldManager();
}
