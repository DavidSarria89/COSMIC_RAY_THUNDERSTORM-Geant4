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

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

int main( int argc, const char *argv[]) {

  int nEv=1000000;

  std::vector<double> dates;
  dates.push_back(2008.0);
  dates.push_back(2008.0+5.5);

  std::vector<double> alt;
  alt.push_back(0.);
  alt.push_back(2100.0);
  alt.push_back(11300.0);

  std::string setupString("");
  std::vector<CRYParticle*> *ev=new std::vector<CRYParticle*>;
  CRYSetup *setup=new CRYSetup(setupString,"./data");

  // Iterate over the dates
  for ( unsigned int d=0; d< dates.size(); d++) {
    setup->setParam(CRYSetup::date,dates[d]);

    // Iterate over the three altitudes
    for ( unsigned int r=0; r< alt.size(); r++) {
      std::cout << "Starting test for altitude: " << alt[r] << "m ";
      std::cout << "and year: " << dates[d] << std::endl;
      setup->setParam(CRYSetup::altitude,alt[r]);
      setup->setParam(CRYSetup::subboxLength,300.0);
      setup->setParam(CRYSetup::latitude,90.0);
      
      // Set all the particle counts to zero
      std::map<CRYParticle::CRYId, int> counts;
      for ( CRYParticle::CRYId i=CRYParticle::CRYIdMin; i<=CRYParticle::CRYIdMax; i=CRYParticle::CRYId(i+1)) {
        counts[i]=0;
      }
      
      // Set up the generator
      CRYGenerator gen(setup);
      
      // Generate the events
      for ( int i=0; i<nEv; i++) {
        ev->clear();
        gen.genEvent(ev);
        
        // Count the particles of each type
        for ( unsigned j=0; j<ev->size(); j++) {
          CRYParticle *part=(*ev)[j];
          counts[part->id()]+=1;
          delete (*ev)[j];
        }
      }
      
      double bsize = setup->param(CRYSetup::subboxLength);
      std::cout << "Box size: " << bsize << "m X " << bsize << "m\n";
      std::cout << "Total time simulated: " << gen.timeSimulated() << " seconds (" << nEv << " events)\n";
      for ( CRYParticle::CRYId i=CRYParticle::CRYIdMin; i<=CRYParticle::CRYIdMax; i=CRYParticle::CRYId(i+1)) {
        if( i != CRYParticle::Kaon)  // skip kaons in output, because we do not have data for them yet
          std::cout << CRYUtils::partName(i) << "s per second per m2: "<< 
            counts[i]/(setup->param(CRYSetup::subboxLength)*setup->param(CRYSetup::subboxLength)*gen.timeSimulated()) << std::endl;
      }
      std::cout << "========================================================\n";
    }
  }
  return 0;
}
