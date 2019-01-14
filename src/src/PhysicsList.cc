//

// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
// $Id: PhysicsList.cc 70268 2013-05-28 14:17:50Z maire $

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option1_dr.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "GammaPhysics.hh"

// particles

#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(const bool em_only) : G4VModularPhysicsList() {
    G4int verb = 0;

    SetVerboseLevel(verb);

    // EM physics
    //    RegisterPhysics(new G4EmStandardPhysics_option1());
    RegisterPhysics(new G4EmStandardPhysics_option1_dr());

    //    RegisterPhysics(new G4EmStandardPhysics_option3());
    //    RegisterPhysics(new G4EmStandardPhysicsSS());

    if (!em_only) {
        // Decay
        RegisterPhysics(new G4DecayPhysics());

        // Radioactive decay
        RegisterPhysics(new G4RadioactiveDecayPhysics());

        // Hadron Elastic scattering
        RegisterPhysics(new G4HadronElasticPhysics(verb));

        //            RegisterPhysics(new G4HadronElasticPhysicsHP(verb));

        // Hadron Inelastic physics
        RegisterPhysics(new G4HadronPhysicsFTFP_BERT(verb));

        //            RegisterPhysics(new G4HadronPhysicsFTFP_BERT_HP(verb));
        //            RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP(verb));
        //            RegisterPhysics(new G4HadronInelasticQBBC(verb));
        //            RegisterPhysics( new G4HadronPhysicsINCLXX(verb));

        // Ion Elastic scattering
        RegisterPhysics(new G4IonElasticPhysics(verb));

        // Ion Inelastic physics
        RegisterPhysics(new G4IonPhysics(verb));

        //            RegisterPhysics(new G4IonINCLXXPhysics(verb));

        // Gamma-Nuclear Physics
        RegisterPhysics(new GammaPhysics("gamma"));
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle() {
    G4BosonConstructor pBosonConstructor;

    G4BosonConstructor::ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    G4LeptonConstructor::ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    G4MesonConstructor::ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    G4BaryonConstructor::ConstructParticle();

    G4IonConstructor pIonConstructor;
    G4IonConstructor::ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    G4ShortLivedConstructor::ConstructParticle();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts() {
    //// Production thresholds for world (default) ; !!!! : That is actually for the zone inside the plane with STP air
    G4double cutvalp = 0.0;
    G4double cutvalr = 1.0 * mm;

    //    G4double cutval=1. * cm ;
    SetCutValue(cutvalp, G4ProductionCuts::GetIndex("proton"));
    SetCutValue(cutvalr, G4ProductionCuts::GetIndex("e-"));
    SetCutValue(cutvalr, G4ProductionCuts::GetIndex("e+"));
    SetCutValue(cutvalr, G4ProductionCuts::GetIndex("gamma"));
    SetCutValue(cutvalr, G4ProductionCuts::GetIndex("mu+"));
    SetCutValue(cutvalr, G4ProductionCuts::GetIndex("mu-"));

    G4EmParameters *EM_params = G4EmParameters::Instance();
    EM_params->SetNumberOfBinsPerDecade(16);

    //// Production thresholds for E-field region
    G4Region *EFIELD_region = G4RegionStore::GetInstance()->GetRegion("EFIELD");
    auto *cuts_EFIELD = new G4ProductionCuts;

    G4double cut_value = 1.0 / 100. * mm;
    cuts_EFIELD->SetProductionCut(cutvalp, G4ProductionCuts::GetIndex("neutron"));
    cuts_EFIELD->SetProductionCut(cutvalp, G4ProductionCuts::GetIndex("proton"));
    cuts_EFIELD->SetProductionCut(cut_value, G4ProductionCuts::GetIndex("e-"));
    cuts_EFIELD->SetProductionCut(cut_value, G4ProductionCuts::GetIndex("e+"));
    cuts_EFIELD->SetProductionCut(cut_value, G4ProductionCuts::GetIndex("gamma"));
    cuts_EFIELD->SetProductionCut(cut_value, G4ProductionCuts::GetIndex("mu-"));
    cuts_EFIELD->SetProductionCut(cut_value, G4ProductionCuts::GetIndex("mu+"));
    EFIELD_region->SetProductionCuts(cuts_EFIELD);

//    G4double lowlimit = 8 * CLHEP::keV;
//    G4ProductionCutsTable *aPCTable = G4ProductionCutsTable::GetProductionCutsTable();
//    aPCTable->SetEnergyRange(lowlimit, 100 * CLHEP::GeV);

    if (settings->USE_GLOBAL_MAX_STEP) {
        Add_global_StepMax();
    }
    if (settings->USE_RECORD_REGION_STEP_MAX) {
        Add_StepMax_for_record_regions();
    }
} // PhysicsList::SetCuts

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "StepMax.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiter.hh"

void PhysicsList::Add_global_StepMax() {
    // Step limitation seen as a process
    StepMax *stepMaxProcess = new StepMax();

    stepMaxProcess->SetMaxStep(GLOBAL_stepMaxVal);

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)()) {
        G4ParticleDefinition *particle = particleIterator->value();
        G4ProcessManager *pmanager = particle->GetProcessManager();

        if (stepMaxProcess->IsApplicable(*particle)) {
            pmanager->AddDiscreteProcess(stepMaxProcess);
        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4StepLimiter.hh"
#include "G4ProcessManager.hh"

void PhysicsList::Add_StepMax_for_record_regions() {
    // Step limitation seen as a process
    G4StepLimiter *stepLimiter = new G4StepLimiter();

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)())  // for all particles
    {
        G4ParticleDefinition *particle = particleIterator->value();
        G4ProcessManager *pmanager = particle->GetProcessManager();
        pmanager->AddDiscreteProcess(stepLimiter);
    }
}