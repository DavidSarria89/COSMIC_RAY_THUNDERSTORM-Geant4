// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#pragma once

#include "AnalysisManager.hh"
#include "G4Box.hh"
#include "G4Navigator.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "fortran.hh"
#include "globals.hh"
#include <locale.h>
#include <numeric>
#include <stdlib.h>
#include "cosmic_ray_generator_PARMA.hh"
#include "Settings.hh"

class G4Event;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:

    PrimaryGeneratorAction();

    ~PrimaryGeneratorAction() override;

public:

    void GeneratePrimaries(G4Event *) override;

    G4ParticleGun *GetParticleGun() {
        return fParticleGun;
    }

private:

    geant4_initial_cosmic_ray out_cosmic;

    Settings *settings = Settings::getInstance();

    Cosmic_Ray_Generator_PARMA *cosmic_ray_gene_p = nullptr;

    G4ParticleGun *fParticleGun; // pointer a to G4 service class

};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
