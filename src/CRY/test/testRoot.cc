/* 

Copyright (c) 2007-2012, The Regents of the University of California. 
Produced at the Lawrence Livermore National Laboratory 
UCRL-CODE-227323. 
All rights reserved. 
 
For details, see http://nuclear.llnl.gov/simulations
Please also read this http://nuclear.llnl.gov/simulations/additional_bsd.html
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
 
1.  Redistributions of source code must retain the above copyright
notice, this list of conditions and the disclaimer below.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the disclaimer (as noted below) in
the documentation and/or other materials provided with the
distribution.

3. Neither the name of the UC/LLNL nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OF
THE UNIVERSITY OF CALIFORNIA, THE U.S. DEPARTMENT OF ENERGY OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include "CRYGenerator.h"
#include "CRYSetup.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include <math.h>
#include <stdlib.h>  // For Ubuntu Linux

int main( int argc, const char *argv[]) {

  if (!TROOT::Initialized()) {
    static TROOT root("RooTuple", "RooTuple ROOT God in CosmicRaY simulation");
  }


  if ( argc < 2 ) {
    std::cout << "usage " << argv[0] << " <setup file name> <N events>\n";
    std::cout << "N events = 10k by default\n";
    return 0;
  }

  TFile *outputFile=new TFile("test.root","RECREATE");
  
  TH1I* multiplicity = new TH1I("mult",       "mult",      100,     0,  100);
  TH1F* latLoc       = new TH1F("latLoc",     "latLoc",    100, -150.0, 150.0);
  TH1F* keMuon       = new TH1F("keMuon",     "keMuon",     5*16, -9.0, 7.0);
  TH1F* keNeutron    = new TH1F("keNeutron",  "keNeutron",  5*16, -9.0, 7.0);
  TH1F* kePion       = new TH1F("kePion",     "kePion",     5*16, -9.0, 7.0);
  TH1F* keKaon       = new TH1F("keKaon",     "keKaon",     5*16, -9.0, 7.0);
  TH1F* keGamma      = new TH1F("keGamma",    "keGamma",    5*16, -9.0, 7.0);
  TH1F* keElectron   = new TH1F("keElectron", "keElectron", 5*16, -9.0, 7.0);
  TH1F* keProton     = new TH1F("keProton",   "keProton",   5*16, -9.0, 7.0);
  TH1F* kePrimary    = new TH1F("kePrimary",  "kePrimary",  5*16, -9.0, 8.0);
  TH1F* chargeHist   = new TH1F("chargeHist", "chargeHist",   10, -0.5, 9.5);
  TH1F* multHist     = new TH1F("multHist",   "multHist",     10, -0.5, 9.5);

  float box_size = 50.;
  int xybins = 50;
  TH1F* xall         = new TH1F("xall",   "xall",   xybins, -box_size, box_size);
  TH1F* yall         = new TH1F("yall",   "yall",   xybins, -box_size, box_size);
  TH1F* xmuon        = new TH1F("xmuon",  "xmuon",  xybins, -box_size, box_size);
  TH1F* ymuon        = new TH1F("ymuon",  "ymuon",  xybins, -box_size, box_size);
  TH2F* xyall        = new TH2F("xyall",  "xyall",  xybins, -box_size, box_size, xybins, -box_size, box_size);
  TH2F* xymuons      = new TH2F("xymuon", "xymuon", xybins, -box_size, box_size, xybins, -box_size, box_size);

  TH1F* costhe       = new TH1F("costhe",  "costhe plot;cos(theta)",   200, -1., 0. );

  std::map<CRYParticle::CRYId,TH1F*> keHistos;

  keHistos[CRYParticle::Neutron]=keNeutron;
  keHistos[CRYParticle::Muon]=keMuon;
  keHistos[CRYParticle::Proton]=keProton;
  keHistos[CRYParticle::Electron]=keElectron;
  keHistos[CRYParticle::Gamma]=keGamma;
  keHistos[CRYParticle::Kaon]=keKaon;
  keHistos[CRYParticle::Pion]=kePion;

  int nEv = 100000;
  if (argc > 2 ) nEv = atoi(argv[2]);

  // Read the setup file into setupString
  std::ifstream inputFile;
  inputFile.open(argv[1],std::ios::in);
  char buffer[1000];

  std::string setupString("");
  while (!inputFile.getline(buffer,1000).eof()) {
    setupString.append(buffer);
    setupString.append(" ");
  }

  // Parse the contents of the setup file
  CRYSetup *setup=new CRYSetup(setupString,"./data");

  // Setup the CRY event generator
  CRYGenerator gen(setup);

  // Generate the events
  int nMuon = 0;
  std::vector<CRYParticle*> *ev=new std::vector<CRYParticle*>;
  for (int i = 0; i < nEv; i++) {
    ev->clear();
    gen.genEvent(ev);

    if (i % 1000 == 0) std::cout << "Event: " << i << std::endl;

    // Fill the histograms
    multiplicity->Fill(ev->size());
    kePrimary->Fill(log10(gen.primaryParticle()->ke()));
    for (unsigned j = 0; j < ev->size(); j++) {
      CRYParticle *part=(*ev)[j];
      
      //....printout all secondaries every 1000 events just for fun
      if (i % 1000 == 0) {
        std::cout << "Secondary " << j << " " <<
        CRYUtils::partName(part->id()) << " ke=" << part->ke() << "\n";
      }
      
      keHistos[part->id()]->Fill(log10(part->ke()));
      latLoc->Fill(sqrt(part->x()*part->x()+part->y()*part->y())); 
      chargeHist->Fill(part->id(),part->charge());
      multHist->Fill(part->id());
      xall->Fill(part->x());
      yall->Fill(part->y());
      xyall->Fill(part->x(), part->y());
      costhe->Fill(part->w());
      
      if (part->id() == CRYParticle::Muon) {
        xmuon->Fill(part->x());
        ymuon->Fill(part->y());      
        xymuons->Fill(part->x(), part->y());
        nMuon++;
      }
      
      delete (*ev)[j];
    }
  }

  chargeHist->Divide(multHist);
  std::cout << "Run completed.\n";
  std::cout << "Total time simulated: " << gen.timeSimulated() << " seconds\n";
  double muonsPerSecondPerm2=nMuon/(300.0*300.0*gen.timeSimulated());
  std::cout << "Muons per second per m2 " << muonsPerSecondPerm2 << std::endl;

  // Write the histogram file
  outputFile->Write();
  outputFile->Close();

  delete setup;

  return 0;
}
