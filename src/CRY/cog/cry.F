C /* 
C 
C Copyright (c) 2007-2012, The Regents of the University of California. 
C Produced at the Lawrence Livermore National Laboratory 
C UCRL-CODE-227323. 
C All rights reserved. 
C  
C For details, see http://nuclear.llnl.gov/simulations
C Please also read this http://nuclear.llnl.gov/simulations/additional_bsd.html
C  
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C  
C 1.  Redistributions of source code must retain the above copyright
C notice, this list of conditions and the disclaimer below.
C 
C 2. Redistributions in binary form must reproduce the above copyright
C notice, this list of conditions and the disclaimer (as noted below) in
C the documentation and/or other materials provided with the
C distribution.
C 
C 3. Neither the name of the UC/LLNL nor the names of its contributors
C may be used to endorse or promote products derived from this software
C without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OF
C THE UNIVERSITY OF CALIFORNIA, THE U.S. DEPARTMENT OF ENERGY OR
C CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
C EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
C PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
C PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
C LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
C NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
C SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C */
C 
C 
      subroutine cry()
c
c  USRSOR subroutine for CRY cosmic ray shower library
c
c  Chris Hagmann Jan 2007 
c
c Outputs:
c User must generate values for these source parameters, which are
c part of the COG source include file Sorcm.inc.
c--integer*4 values----
c   Sntype = particle ID.
c--Sntype Values----
c Sntype     Particle
c  1          neutron
c  2          proton
c  3          deuteron
c  4          triton
c  5          he-3
c  6          alpha
c  7          photon
c  8          electron
c  9          positron
c
c--real*8 values----
c   Sx,Sy,Sz = Source particle location.
c   Su,Sv,Sw = Direction cosines of particle's velocity vector.
c   Seng = Particle kinetic energy.
c   Svel = Particle velocity.
c   Sage = Particle age at time of creation.
c   Swate = Particle weight.
c   xnors = Source Strength.
c
c COG routine invoked:
c  RNG (random number generator)
c
c *******BEGIN DECLARATIONS*******************************************
      implicit none
c
c COG Include files needed.
#include "basic.inc"
#include "source.inc"
#include "spav.inc"
#include "partyp.inc"
#include "Sorcm.inc"
#include "sorpar.inc"
#include "filenames.inc"
#include "FlagsEtc.inc"
#include "PhysConst.inc"
c
c Include Statement for MPI Information
       include '/usr/gapps/dmc/work/verbeke2/cog_new/include/mpif.h'
c
c
c     --- external functions ---
      real*4 RNG
      real*8 drng

      integer*4 crytodo 
      integer*4 ierr
      integer*4 cry_init
      integer*4 cry_smp

      integer*4 npmax

      parameter (npmax=1000)

      double precision cry_e(npmax)
      double precision cry_c(npmax)
      double precision cry_x(npmax)
      double precision cry_y(npmax)
      double precision cry_z(npmax)
      double precision cry_u(npmax)
      double precision cry_v(npmax)
      double precision cry_w(npmax)
      double precision cry_t(npmax)

      integer*4 cry_pid(npmax)
c     
c     
c     --- local variables ---
      integer*4 i
      integer*4 firstflag       ! == 1 on fisrt call
 
      data firstflag /1/
      save firstflag


      data crytodo/0/
      save crytodo

      save cry_e,cry_c,cry_u,cry_v,cry_w,cry_x,cry_y,cry_z,cry_t,cry_pid

      external rng
      external drng
c
c *******END DECLARATIONS*******************************************
c
      if (firstflag .eq. 1) then
         firstflag = 0
         ierr= cry_init(drng)   
      endif

      Swate = 1.  ! particle weight


111   if (crytodo.eq.0) then
c new shower from cry
      crytodo = cry_smp(cry_e,cry_x,cry_y,cry_z,cry_u,cry_v,cry_w,
     &               cry_t,cry_pid,npmax)

      endif

      if(cry_pid(crytodo).eq.1.and.kp(1).ne.0) then  ! neutron tracked ?
           Sntype = 1
           Seng = cry_e(crytodo)
           Svel = cspeed * 
     &     dsqrt(Seng*(Seng+2.*rme(Sntype)))/(Seng+rme(Sntype))
      else if(cry_pid(crytodo).eq.7.and.kp(7).ne.0) then   ! photon tracked ?
           Sntype = 7
           Seng = cry_e(crytodo)
           Svel = cspeed
      else     ! other (Not tracked)
           crytodo = crytodo-1 
           goto 111
      endif
      
      Su = cry_u(crytodo)
      Sv = cry_v(crytodo)
      Sw = cry_w(crytodo)
      Sx = cry_x(crytodo)
      Sy = cry_y(crytodo)
      Sx = cry_z(crytodo)
      Sage = cry_t(crytodo)
      crytodo = crytodo - 1


      return
      end



      real*8 function drng(rng)

      real*4 rng
      external rng
      
      drng = dble(rng())

      return
      end
      
