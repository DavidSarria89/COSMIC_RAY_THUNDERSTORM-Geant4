#pragma once

#include "AnalysisManager.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "fortran.hh"
#include "globals.hh"
#include <locale.h>
#include <numeric>
#include <stdlib.h>
#include "Settings.hh"
#include <random>
#include "CRYSetup.h"
#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"
#include "RNGWrapper.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Cosmic_Ray_Generator_CRY {
public:

    Cosmic_Ray_Generator_CRY(const char *inputfile);

    ~Cosmic_Ray_Generator_CRY();

    geant4_initial_cosmic_ray generate_Cosmic_ray();

private:
    G4ParticleTable *particleTable;
    G4int InputState;

    Settings *settings = Settings::getInstance();

    // all the particle types that it is possible to have
    G4ParticleDefinition *Gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
    G4ParticleDefinition *Electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    G4ParticleDefinition *Positron = G4ParticleTable::GetParticleTable()->FindParticle("e+");
    G4ParticleDefinition *Neutron = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    G4ParticleDefinition *Proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    G4ParticleDefinition *MuonP = G4ParticleTable::GetParticleTable()->FindParticle("mu+");
    G4ParticleDefinition *MuonN = G4ParticleTable::GetParticleTable()->FindParticle("mu-");

    void InputCRY();

    void UpdateCRY(std::string *MessInput);

    void CRYFromFile(G4String newValue);

    std::vector<CRYParticle *> *vect;
    CRYGenerator *gen;

    double altitude_generator = 11.3*km;

};
