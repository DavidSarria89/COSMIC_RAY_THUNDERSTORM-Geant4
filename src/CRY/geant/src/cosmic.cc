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


//******************************************************************************
// cosmic.cc  GEANT4 user application for testing handling of
//            cosmic-rays
//
// 1.00 JMV, LLNL, JAN-2007:  First version.
//******************************************************************************
//

// misc includes
//
#include <fstream>
#include <math.h>
#include "G4ios.hh"

// package includes
//
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

// geant4 includes
//
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

//------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Run manager
  //------------
  G4RunManager * theRunManager = new G4RunManager;

  // UserInitialization classes
  //---------------------------
  DetectorConstruction* theDetector = new DetectorConstruction;
  PhysicsList* thePhysicsList = new PhysicsList;

  // UserAction classes
  //-------------------

  theRunManager->SetUserInitialization(theDetector);
  theRunManager->SetUserInitialization(thePhysicsList);
  theRunManager->SetUserAction(new PrimaryGeneratorAction(""));

  // Initialize G4 kernel
  //---------------------
  theRunManager->Initialize();

  // User interactions
  //------------------
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  if (argc > 1){  //....geant command file specified on command line
    UI->ApplyCommand("/control/execute "+G4String(argv[1]));

  }else{           //....no command file specified, Start interactive session 
    G4UIsession * theUIsession = new G4UIterminal(new G4UItcsh);
    theUIsession->SessionStart();
    delete theUIsession;
  }

  delete theRunManager;
  return EXIT_SUCCESS;
}
