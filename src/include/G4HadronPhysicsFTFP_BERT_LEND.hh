#pragma once

#include "G4ios.hh"
#include "globals.hh"

#include "G4VPhysicsConstructor.hh"

#include "G4BertiniPiKBuilder.hh"
#include "G4FTFPPiKBuilder.hh"
#include "G4PiKBuilder.hh"

#include "G4BertiniProtonBuilder.hh"
#include "G4FTFPNeutronBuilder.hh"
#include "G4FTFPProtonBuilder.hh"
#include "G4ProtonBuilder.hh"

#include "G4NeutronBuilder.hh"

#include "G4BertiniNeutronBuilder.hh"

#include "G4FTFPNeutronBuilder.hh"

#include "G4NeutronLENDBuilder.hh"
#include "G4NeutronPHPBuilder.hh"

#include "G4AntiBarionBuilder.hh"
#include "G4FTFPAntiBarionBuilder.hh"
#include "G4HyperonFTFPBuilder.hh"

#include "G4LENDCapture.hh"
#include "G4LENDCaptureCrossSection.hh"
#include "G4LENDFission.hh"
#include "G4LENDFissionCrossSection.hh"
#include "G4LENDInelastic.hh"

#include "G4LENDInelasticCrossSection.hh"

class G4ComponentGGHadronNucleusXsc;

class G4HadronPhysicsFTFP_BERT_LEND: public G4VPhysicsConstructor
{
public:

  G4HadronPhysicsFTFP_BERT_LEND(G4int verbose = 1);

  G4HadronPhysicsFTFP_BERT_LEND(const G4String &name, G4bool quasiElastic = false);

  virtual ~G4HadronPhysicsFTFP_BERT_LEND();

public:

  virtual void ConstructParticle();

  virtual void ConstructProcess();

private:

  void               CreateModels();

  G4HadronicProcess* FindInelasticProcess(const G4ParticleDefinition *);

  struct ThreadPrivate
  {
    G4NeutronBuilder        *theNeutrons;
    G4BertiniNeutronBuilder *theBertiniNeutron;
    G4FTFPNeutronBuilder    *theFTFPNeutron;
    G4NeutronLENDBuilder    *theHPNeutron;

    G4PiKBuilder        *thePiK;
    G4BertiniPiKBuilder *theBertiniPiK;
    G4FTFPPiKBuilder    *theFTFPPiK;

    G4ProtonBuilder        *thePro;
    G4BertiniProtonBuilder *theBertiniPro;
    G4FTFPProtonBuilder    *theFTFPPro;

    G4HyperonFTFPBuilder *theHyperon;

    G4AntiBarionBuilder     *theAntiBaryon;
    G4FTFPAntiBarionBuilder *theFTFPAntiBaryon;

    G4ComponentGGHadronNucleusXsc *xsKaon;
    G4VCrossSectionDataSet        *xsNeutronCaptureXS;

    G4LENDFissionCrossSection *theLENDFissionCrossSection;

    G4LENDCaptureCrossSection *theLENDCaptureCrossSection;

    G4LENDInelastic             *theLENDInelastic;
    G4LENDInelasticCrossSection *theLENDInelasticCrossSection;
  };
  static G4ThreadLocal ThreadPrivate *tpdata;

  // G4VCrossSectionDataSet * BGGProton;
  // G4VCrossSectionDataSet * BGGNeutron;
  G4bool QuasiElastic;
};
