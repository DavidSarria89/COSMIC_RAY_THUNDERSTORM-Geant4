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


//
// Compare old and new version of same histogram stored in differenct files
// Loads root files auotmatically so, do:
//
//   root compare_root.cc
//
// Method for looping over keys copied from
// http://www.hep.wisc.edu/~cmsprod/farmoutCmsJobs/mergeFiles.C
//
// Doug Wright

{
  TFile* fref = new TFile("test_ref.root");
  TFile* fnew = new TFile("test.root");

  TList* keys = new TList();
  

  // Make a list of all the 1D histograms that appear in both files
  TKey* keyref;
  TKey* keynew;
  TIter nextref(fref->GetListOfKeys());
  while((keyref = (TKey*)nextref())) {
    TIter nextnew(fnew->GetListOfKeys());
    while((keynew = (TKey*)nextnew())) {
      //      if (strcmp(keyref->GetName(), keynew->GetName()) == 0) {
      if( TString(keyref->GetName()) == TString(keynew->GetName()) ) {  // The keys match
        if( TString(keynew->GetClassName()).BeginsWith("TH1") ) {      // Only include 1D histograms
          keys->AddLast(keynew);
          break;
        }
      }
    }
  }

  
  if(keys->GetSize() < 1) {
    std:cout << "There are no matching histograms " << endl;
    return:
  }
  
  //std:cout << "number of plots: " << keys->GetSize() << endl;
  
  //....make enough zones to plot all histograms
  int n = keys->GetSize();
  int j = ceil(sqrt(float(n)));
  int k = (n+j-1)/j;  //....adding j-1 to the numerator makes this round up which is what I want
                      //    n/j would round down, and thus might not be enough zones

  //....make large canvas with jxk plots
  TCanvas* c = new TCanvas();
  c->SetWindowPosition(0,0); //....upper left-hand of screen
  c->SetWindowSize(1600.,860.);
  c->Divide(k,j,0.001,0.001);
  c->cd(1);

  // loop over all keys in the reference file and plot all histograms
  TIter nextkey( keys );
  TKey* key;
  while ((key = (TKey*)nextkey())){
    char* name = key->GetName();
    fref->Get(name)->Draw();
    fnew->Get(name)->Draw("Esame");
    c->cd(gPad->GetNumber()+1);
  }

  //....set first plot to have log y axis
  c->cd(1); 
  gPad->SetLogy() ;

  gPad->Update();
  c->Print("compare_root.pdf");
}
