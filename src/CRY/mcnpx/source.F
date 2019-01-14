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
      subroutine source
!
! Call LLNL Cosmic-Ray Physics (CRY) library to generate secondaries
! produced by cosmic-ray showers.
!

       use global_data







       implicit real(kindr) (a-h,o-z)
       implicit integer(kindi) (i-n)





      integer*4 npmax
      parameter (npmax=1000)

      integer*4 pidmap(0:6)

      real*8 cry_e(npmax) ! kin. E (MeV)
      real*8 cry_x(npmax) ! x (m)
      real*8 cry_y(npmax)
      real*8 cry_z(npmax)
      real*8 cry_u(npmax)
      real*8 cry_v(npmax)
      real*8 cry_w(npmax)
      real*8 cry_t(npmax)

      integer*4 cry_init
      integer*4 cry_smp
      integer*4 cry_pid(npmax) ! particle id
      integer*4 cry_charge(npmax) ! particle charge
      integer*4 itracked(npmax) 
      integer*4 ierr
      integer*4 ifirst
      integer*4 npar
      integer*4 npartracked
      integer*4 ip,kp

! arg= CRY pid: neutron=0,proton=1,pion=2,kaon=3,muon=4,electron=5,gamma=6
! map= mcnpx pid

      data pidmap /1,9,20,22,4,3,2/

      data ifirst /1/


      external rang

      save ifirst

! source cell # = idum(1)
! z pos of cosmic ray = rdum(1)

      if(ifirst.eq.1)then
        ifirst = 0
        ierr = cry_init(rang)
      endif

      icl = namchg(1,idum(1)) ! convert from user cell to internal cell #
      zzz = rdum(1) ! z start position in cm

      jsu = 0
      wgt = 1.0d0
      npa = 1

!....Copy cutoffs and importance settings from original source cell to the
! current cell
      do i=1,mipt
         elc(i)=elp(i,icl) ! elc=e-cutoffs in current cell, elp = cell dep cutoffs
         fiml(i)=fim(i,icl) ! fiml= importances in current cell, fim = cell dep importances
      enddo


    1 npar = cry_smp(cry_e,cry_x,cry_y,cry_z,cry_u,cry_v,cry_w,
     & cry_t,cry_pid,cry_charge,npmax)





      npartracked = 0

      do 10 kp=1,npar

        erg = cry_e(kp)
        ipt = pidmap(cry_pid(kp))

! check if particle is tracked
        if(ipt .gt. mipt ) goto 10 ! MCNP can only handle n,e,photon
        if(kpt(ipt).eq.0) goto 10 ! discard particles not tracked
        if(erg.gt.emx(ipt)) goto 10 ! discard particles above max E

        npartracked = npartracked+1
        itracked(npartracked) = kp
10    continue

      if(npartracked.eq.0) goto 1
 
 
      do 20 ip=1,npartracked
        kp = itracked(ip)
 
        ipt = pidmap(cry_pid(kp))
        erg = cry_e(kp)
        jan = cry_charge(kp)

! special cases:

! neutron
        if(ipt.eq.1) jan = 1

! electron
        if(ipt.eq.3) jan = -jan

! proton
        if(ipt.eq.9) jan = -jan



        if (ipt.eq.2) then
          vel = slite
        else
          vel = slite*dsqrt(erg*(erg+2.*gpt(ipt)))/(erg+gpt(ipt))
        endif

        xxx= 100.*cry_x(kp) ! m -> cm
        yyy= 100.*cry_y(kp)

        tme= 1.e8*cry_t(kp) ! sec -> shake

        uuu= cry_u(kp)
        vvv= cry_v(kp)
        www= cry_w(kp)

        if (ip.lt.npartracked) then
           nx = igain_all_pls_source
           paxtc(1,nx,ipt) = paxtc(1,nx,ipt) + 1
           paxtc(2,nx,ipt) = paxtc(2,nx,ipt) + wgt
           paxtc(3,nx,ipt) = paxtc(3,nx,ipt) + (wgt*erg)
           pac(kpac+ipt,1,icl) = pac(kpac+ipt,1,icl)+one
           call bankit(5)
         endif

20    continue

      return
      end
