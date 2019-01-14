subroutine gen_parma_cr(seed,nb,energies,cosangles,altitude_list,types)
!  Generate cosmic-ray based on PARMA model
USE parma, ONLY : getr,getd,getspecangfinal,getspec,gethp
implicit none
integer,parameter :: npart=33 ! number of applicable particle
integer,parameter :: nebin=1024 ! number of energy mesh (divided by log)
integer,parameter :: nabin=128 ! number of angle mesh (divided by linear)
integer,parameter :: naltbin=256
real(8) :: ehigh(0:nebin),emid(nebin) ! higher and middle point of energy bin
real(8) :: ahigh(0:nabin),amid(nabin) ! higher and middle point of angular bin
real(8) :: althigh(0:naltbin),altmid(naltbin)
real(8) :: etable(0:nebin)            ! probability table (0.0 for 0, 1.0 for nebin)
real(8) :: atable(0:nabin,0:nebin)    ! probability table (0.0 for 0, 1.0 for nabin)

real(8) :: alt_table(0:naltbin)            ! probability table
real(8) :: e_table(0:naltbin,0:nebin)            ! probability table
real(8) :: a_table(0:naltbin,0:nabin,0:nebin)    ! probability table

integer :: IangPart(0:33)=(/1,2,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,5,5,6/) ! Angular particle ID

real(8) :: PPII = acos(-1.0)
real(8) :: spec_val = 0.0,spec_angle=0.0

integer,intent(in) :: nb
integer,intent(in) :: seed
real(8),intent(out) :: energies(nb)
real(8),intent(out) :: cosangles(nb)
real(8),intent(out) :: altitude_list(nb)
integer,intent(out) :: types(nb)

integer :: ip,nevent,iyear,imonth,iday,i,i_alt,ia,ie,idummy
real(8) :: Alti,glat,glong,g,emin,emax,amin,amax,altmax,altmin,amid_ia,emid_ie
real(8) :: elog,estep,astep,alt_step,s,rr,d,e,w,totalflux,sampled_alt
real(8) :: getgeneration

call srand(seed)

! Set number of particle to be generated, and initial random seed
nevent=nb ! number of particles to be generated

! Set Particle ID
ip=33 ! Particle ID (Particle ID, 0:neutron, 1-28:H-Ni, 29-30:muon+-, 31:e-, 32:e+, 33:photon)
Alti=12.0   ! Altitude (km)

! Set Conditions (location, time, and local geometry)
idummy = 0
iyear=2016  ! Year
imonth=3    ! Month
iday=20      ! Date
glat=-13.0   ! Latitude (deg), -90 =< glat =< 90
glong=130.0 ! Longitude (deg), -180 =< glong =< 180
g=10.0      ! Local geometry parameter, 0=< g =< 1: water weight fraction, 10:no-earth, 100:blackhole, -10< g < 0: pilot, g < -10: cabin



! Set energy and angle ranges for generation
emin = 2.0e-2  ! Minimum energy of particle
emax = 1.0e6  ! Maximum energy of particle
amin = -1.0   ! Minimum cosine of particle
amax =  1.0   ! Maximum cosine of particle
altmax = 40.0 ! km
altmin = 7.0 ! km
if(ip.eq.0.and.emin.lt.1.0e-8) emin=1.0e-8 ! Minimum energy for neutron is 10 meV
if(ip.ne.0.and.emin.lt.1.0e-2) emin=1.0e-2 ! Minimum energy for other particle is 10 keV

! Make energy and angle mesh and altitude mesh
elog=log10(emin)
estep=(log10(emax)-log10(emin))/nebin
do ie=0,nebin
 ehigh(ie)=10d0**elog
 if(ie.ne.0) emid(ie)=sqrt(ehigh(ie)*ehigh(ie-1))
 elog=elog+estep
enddo

astep=(amax-amin)/nabin
do ia=0,nabin
 ahigh(ia)=amin+astep*ia
 if(ia.ne.0) amid(ia)=(ahigh(ia)+ahigh(ia-1))*0.5
enddo

alt_step=(altmax-altmin)/naltbin
do i_alt=0,naltbin
 althigh(i_alt)=altmin+alt_step*i_alt
 if(i_alt.ne.0) altmid(i_alt)=(althigh(i_alt)+althigh(i_alt-1))*0.5
enddo

! Make probability table for particle type

a_table(:,:,:)=0.0d0 ! initialization
e_table(:,:)=0.0d0
alt_table(:)=0.0d0

do i_alt=1,naltbin
 s=41.0 ! Solar activity (W index) of the day
 rr=14.0               ! Vertical cut-off rigidity (GV)
 d=getd(altmid(i_alt),glat)                 ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

 do ie=1,nebin
  do ia=1,nabin
   emid_ie = emid(ie)
   spec_val = getSpec(ip,s,rr,d,emid_ie,g)
   write(*,*) spec_val
   amid_ia = amid(ia)
   if (amid_ia.gt.1.0) amid_ia = 0.999999
   if (amid_ia.lt.-1.0) amid_ia = -0.999999

   spec_angle = getSpecAngFinal(iangpart(ip),s,rr,d,emid_ie,g,amid_ia)

   a_table(i_alt,ia,ie) = a_table(i_alt,ia-1,ie) + spec_val * spec_angle &
   * (2.0*PPII) * (ahigh(ia)-ahigh(ia-1))

  enddo
  e_table(i_alt,ie)=e_table(i_alt,ie-1)+a_table(i_alt,nabin,ie)*(ehigh(ie)-ehigh(ie-1))
 enddo

alt_table(i_alt)=alt_table(i_alt-1)+e_table(i_alt,nebin)*(althigh(i_alt)-althigh(i_alt-1))

enddo




sampled_alt = 12.0

s=getHP(iyear,imonth,iday,idummy) ! Solar activity (W index) of the day
rr=getr(glat,glong)                ! Vertical cut-off rigidity (GV)
d=getd(sampled_alt,glat)                 ! Atmospheric depth (g/cm2), set glat = 100 for use US Standard Atmosphere 1976.

! Make probability table for energy and angle (absolute value)
atable(:,:)=0.0d0 ! initialization
etable(:)=0.0d0

do ie=1,nebin
 do ia=1,nabin
  spec_val = getSpec(ip,s,rr,d,emid(ie),g)
  spec_angle = getSpecAngFinal(iangpart(ip),s,rr,d,emid(ie),g,amid(ia))
  atable(ia,ie) = atable(ia-1,ie) + spec_val * spec_angle * (2.0*PPII) * (ahigh(ia)-ahigh(ia-1)) ! angular integrated value
 enddo
enddo

do ie=1,nebin
 etable(ie)=etable(ie-1)+atable(nabin,ie)*(ehigh(ie)-ehigh(ie-1)) ! energy integrated value
enddo

TotalFlux=etable(nebin) ! Total Flux (/cm2/s), used for normalization

! Make probability table (normalized to 1)
do ie=1,nebin
 etable(ie)=etable(ie)/etable(nebin)
 do ia=1,nabin
  atable(ia,ie)=atable(ia,ie)/atable(nabin,ie)
 enddo
enddo

! Particle Generation
do i=1,nevent
 e=getGeneration(ie,nebin,ehigh(0),etable(0))    ! energy
 energies(i) = e

 w=getGeneration(ia,nabin,ahigh(0),atable(0,ie)) ! z direction, -1.0:upward, 0.0:horizontal, 1.0:downward
 cosangles(i) = w;

types(i) = 33
altitude_list(i) = Alti

enddo

!open(2,file='CheckGeneration.out')
!write(2,*) 'Total Flux (/cm2/s) =',TotalFlux
!write(2,'(102i12)') (ia,ia=1,nabin+2)
!write(2,'(''   Emid(MeV)  /cm2/s/MeV'',100f12.4)') (amid(ia),ia=1,nabin)
!do ie=1,nebin
! write(2,'(102es12.4)') emid(ie),(flux(ia,ie)/nevent,ia=0,nabin)
!enddo

end subroutine gen_parma_cr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function getGeneration(ibin,nbin,high,table)
implicit real*8 (a-h, o-z)
dimension high(0:nbin)
dimension table(0:nbin)

randd=rand() ! random number
do i=1,nbin-1
 if(randd.le.table(i)) exit
enddo
ibin=i ! bin ID

randd=rand() ! random number
getGeneration = high(ibin-1)*randd + high(ibin)*(1.0d0-randd)

return

end
