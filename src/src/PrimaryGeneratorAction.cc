// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(nullptr) {

    if (settings->chosen_CR_GEN_ENGINE == settings->CR_GEN_ENGINE::PARMA) {
        cosmic_ray_gene_p = new Cosmic_Ray_Generator_PARMA();
    } else if (settings->chosen_CR_GEN_ENGINE == settings->CR_GEN_ENGINE::CRY) {

        if (settings->CR_GENERATION_ALT_MAX > 11.3) {
            G4cout << "ERROR: altitude of cosmic ray generation cannot bet more than 11.3 km for the CRY model. Aborting." << G4endl;
            std::abort();
        }

        cosmic_ray_gene_c = new Cosmic_Ray_Generator_CRY("setup_CRY.file");
    }

    G4int n_particle = 1;

    fParticleGun = new G4ParticleGun(n_particle);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {


    if (settings->chosen_CR_GEN_ENGINE == settings->CR_GEN_ENGINE::PARMA) {
        out_cosmic = cosmic_ray_gene_p->generate_Cosmic_ray();
    } else if (settings->chosen_CR_GEN_ENGINE == settings->CR_GEN_ENGINE::CRY) {
        out_cosmic = cosmic_ray_gene_c->generate_Cosmic_ray();
    }

    fParticleGun->SetParticleEnergy(out_cosmic.energy);
    fParticleGun->SetParticlePosition(out_cosmic.position_ini);
    fParticleGun->SetParticleMomentumDirection(out_cosmic.momentum_ini);
    fParticleGun->SetParticleTime(out_cosmic.time);
    fParticleGun->SetParticleDefinition(out_cosmic.g4_particle);
    fParticleGun->GeneratePrimaryVertex(anEvent);
} // PrimaryGeneratorAction::GeneratePrimaries
