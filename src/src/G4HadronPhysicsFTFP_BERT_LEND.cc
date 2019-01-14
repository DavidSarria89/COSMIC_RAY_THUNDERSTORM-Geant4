#include "G4HadronPhysicsFTFP_BERT_LEND.hh"

#include "G4SystemOfUnits.hh"

#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4LFission.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4NeutronRadCapture.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"

//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_LEND);

G4ThreadLocal G4HadronPhysicsFTFP_BERT_LEND::ThreadPrivate *G4HadronPhysicsFTFP_BERT_LEND::tpdata = 0;

G4HadronPhysicsFTFP_BERT_LEND::G4HadronPhysicsFTFP_BERT_LEND(G4int) : G4VPhysicsConstructor("hInelastic FTFP_BERT_HP")

/*    , theNeutrons(0)
 * , theBertiniNeutron(0)
 * , theFTFPNeutron(0)
 * , theHPNeutron(0)
 * , thePiK(0)
 * , theBertiniPiK(0)
 * , theFTFPPiK(0)
 * , thePro(0)
 * , theBertiniPro(0)
 * , theFTFPPro(0)
 * , theHyperon(0)
 * , theAntiBaryon(0)
 * , theFTFPAntiBaryon(0)*/
  , QuasiElastic(false)

/*, xsKaon(0)
 * , xsNeutronCaptureXS(0)*/
{}

G4HadronPhysicsFTFP_BERT_LEND::G4HadronPhysicsFTFP_BERT_LEND(const G4String& name, G4bool quasiElastic) : G4VPhysicsConstructor(name)

/*    , theNeutrons(0)
 * , theBertiniNeutron(0)
 * , theFTFPNeutron(0)
 * , theHPNeutron(0)
 * , thePiK(0)
 * , theBertiniPiK(0)
 * , theFTFPPiK(0)
 * , thePro(0)
 * , theBertiniPro(0)
 * , theFTFPPro(0)
 * , theHyperon(0)
 * , theAntiBaryon(0)
 * , theFTFPAntiBaryon(0)*/
  , QuasiElastic(quasiElastic)

/*, xsKaon(0)
 * , xsNeutronCaptureXS(0)*/
{}

void G4HadronPhysicsFTFP_BERT_LEND::CreateModels()
{
  tpdata->theNeutrons    = new G4NeutronBuilder(true); // Fission on
  tpdata->theFTFPNeutron = new G4FTFPNeutronBuilder(QuasiElastic);
  tpdata->theNeutrons->RegisterMe(tpdata->theFTFPNeutron);
  tpdata->theNeutrons->RegisterMe(tpdata->theBertiniNeutron = new G4BertiniNeutronBuilder);
  tpdata->theBertiniNeutron->SetMinEnergy(19.9 * MeV);
  tpdata->theBertiniNeutron->SetMaxEnergy(5 * GeV);
  tpdata->theNeutrons->RegisterMe(tpdata->theHPNeutron = new G4NeutronLENDBuilder);
  tpdata->thePro     = new G4ProtonBuilder;
  tpdata->theFTFPPro = new G4FTFPProtonBuilder(QuasiElastic);
  tpdata->thePro->RegisterMe(tpdata->theFTFPPro);
  tpdata->thePro->RegisterMe(tpdata->theBertiniPro = new G4BertiniProtonBuilder);
  tpdata->theBertiniPro->SetMaxEnergy(5 * GeV);
  tpdata->thePiK     = new G4PiKBuilder;
  tpdata->theFTFPPiK = new G4FTFPPiKBuilder(QuasiElastic);
  tpdata->thePiK->RegisterMe(tpdata->theFTFPPiK);
  tpdata->thePiK->RegisterMe(tpdata->theBertiniPiK = new G4BertiniPiKBuilder);
  tpdata->theBertiniPiK->SetMaxEnergy(5 * GeV);
  tpdata->theHyperon    = new G4HyperonFTFPBuilder;
  tpdata->theAntiBaryon = new G4AntiBarionBuilder;
  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon = new G4FTFPAntiBarionBuilder(QuasiElastic));
}

G4HadronPhysicsFTFP_BERT_LEND::~G4HadronPhysicsFTFP_BERT_LEND()
{
  if (!tpdata)
  {
    return;
  }

  delete tpdata->theNeutrons;
  delete tpdata->theBertiniNeutron;
  delete tpdata->theFTFPNeutron;
  delete tpdata->theHPNeutron;
  delete tpdata->thePiK;
  delete tpdata->theBertiniPiK;
  delete tpdata->theFTFPPiK;
  delete tpdata->thePro;
  delete tpdata->theBertiniPro;
  delete tpdata->theFTFPPro;
  delete tpdata->theHyperon;
  delete tpdata->theAntiBaryon;
  delete tpdata->theFTFPAntiBaryon;
  delete tpdata;
  tpdata = nullptr;
}

void G4HadronPhysicsFTFP_BERT_LEND::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;

  G4MesonConstructor::ConstructParticle();
  G4BaryonConstructor pBaryonConstructor;
  G4BaryonConstructor::ConstructParticle();
  G4ShortLivedConstructor pShortLivedConstructor;
  G4ShortLivedConstructor::ConstructParticle();
}

#include "G4ProcessManager.hh"

void G4HadronPhysicsFTFP_BERT_LEND::ConstructProcess()
{
  if (tpdata == nullptr)
  {
    tpdata = new ThreadPrivate;
  }

  CreateModels();
  tpdata->theNeutrons->Build();
  tpdata->thePro->Build();
  tpdata->thePiK->Build();

  // --- Kaons ---
  tpdata->xsKaon = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet *kaonxs = new G4CrossSectionInelastic(tpdata->xsKaon);
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(kaonxs);
  tpdata->theHyperon->Build();
  tpdata->theAntiBaryon->Build();

  // --- Neutrons ---
  G4HadronicProcess *capture  = nullptr;
  G4HadronicProcess *fission  = nullptr;
  G4ProcessManager  *pmanager = G4Neutron::Neutron()->GetProcessManager();
  G4ProcessVector   *pv       = pmanager->GetProcessList();

  for (size_t i = 0; i < static_cast<size_t>(pv->size()); ++i)
  {
    if (fCapture == ((*pv)[i])->GetProcessSubType())
    {
      capture = dynamic_cast<G4HadronicProcess *>((*pv)[i]);
    }
    else if (fFission == ((*pv)[i])->GetProcessSubType())
    {
      fission = dynamic_cast<G4HadronicProcess *>((*pv)[i]);
    }
  }

  if (!capture)
  {
    capture = new G4HadronCaptureProcess("nCapture");
    pmanager->AddDiscreteProcess(capture);
  }

  tpdata->xsNeutronCaptureXS = (G4NeutronCaptureXS *)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name());
  capture->AddDataSet(tpdata->xsNeutronCaptureXS);
  capture->AddDataSet(new G4ParticleHPCaptureData);
  auto *theNeutronRadCapture = new G4NeutronRadCapture();
  theNeutronRadCapture->SetMinEnergy(19.9 * MeV);
  capture->RegisterMe(theNeutronRadCapture);

  if (!fission)
  {
    fission = new G4HadronFissionProcess("nFission");
    pmanager->AddDiscreteProcess(fission);
  }

  G4LFission *theNeutronLEPFission = new G4LFission();
  theNeutronLEPFission->SetMinEnergy(19.9 * MeV);
  fission->RegisterMe(theNeutronLEPFission);
} // G4HadronPhysicsFTFP_BERT_LEND::ConstructProcess
