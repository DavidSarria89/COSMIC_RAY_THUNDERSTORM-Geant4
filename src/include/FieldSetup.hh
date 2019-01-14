
#pragma once

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4DormandPrince745.hh"

#include "G4CachedMagneticField.hh"
#include "G4UniformElectricField.hh"
#include "Settings.hh"
#include "G4EqMagElectricField.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4UniformElectricField_custom.hh"

class G4FieldManager;

class G4ChordFinder;

class G4Mag_UsualEqRhs;

class G4MagIntegratorStepper;


class FieldSetup {
public:
    explicit FieldSetup();           //  A zero field
    virtual ~FieldSetup();

    void SetStepper();

    void SetMinStep(G4double sss) { fMinStep = sss; }

    void UpdateField();

    G4FieldManager *GetLocalFieldManager() { return fLocalFieldManager; }

private:
    Settings *settings = Settings::getInstance();

protected:

    double EFIELD_mag = 0;

    // Find the global Field Manager

    G4FieldManager *GetGlobalFieldManager();

    G4MagInt_Driver *EIntgrDriver;

    G4FieldManager *fFieldManager;
    G4FieldManager *fLocalFieldManager;
    G4ChordFinder *fChordFinder;
    G4ChordFinder *fLocalChordFinder;
    G4EqMagElectricField *fEquation;
    G4ElectroMagneticField *fElectricField;
    G4ElectroMagneticField *fLocalElectricField;

    G4MagIntegratorStepper *fStepper;
    G4MagIntegratorStepper *fLocalStepper;

    G4int nvar = 8; // to integrate time and energy

    G4double fMinStep = 1. * micrometer; // minimal step, micrometer

};
