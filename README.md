COSMIC_RAY_THUNDERSTORM-Gean4
=======
Geant4-based code for simulation effect of thunderstorms electric fields on cosmic ray fluxes.
=======

## Generalities
* Code based on Geant4 (particle propagation in arbitrary materials, with arbitrary electric-fields), NRL-MSISE-00 (atmosphere model) and PARMA or CRY (cosmic-ray model) to simulate the effect of thunderstorms electric fields on cosmic ray fluxes (electrons, positrons, muons, photons, etc.).
* Integrates the [NRL-MSISE-00 model](https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/) for the atmosphere. The cosmic-ray particles can be generated using [PARMA](https://phits.jaea.go.jp/expacs/) or [CRY](https://nuclear.llnl.gov/simulation/main.html) (this last cannot generate particules at an altitude higher than 11.3 km).
* The code itself mostly uses Geant4 features. See [documentation](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/ForApplicationDeveloper/html/index.html "Geant4 documentation").
* This output is written as an ASCII file in the folder `build/output` as histograms of energies and momentums, for photons, electrons and positrons. The file `build/OUTPUT_FILE_DESCRIPTION.md` explains the output files' content.
* Feel free to clone, fork and suggest improvements for the code.

## Requirements
* Linux-based OS
* GCC compiler suite (C, C++, Fortran)
* Requires a Geant4 installation (library and headers), see [the website](http://geant4.web.cern.ch/) and [installation instructions](http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html/index.html). Recommended version is 10.6.3. For Linux and Ubuntu users, easy installation scripts are provided in the [following repository](https://github.com/DavidSarria89/GEANT4-easy-install-script).

## Compilation
* The source code is localted in `src/` and the build should be done in the folder `build/`.
* First, go to the folder `src/CRY/`, open a terminal and run the command `make` to compile CRY.
* Then, open a terminal in the `build/` directory and type the commands `cmake ../` and then `make` to compile the code. It produces the executable `build/mos_test`

## Usage
* The executable `mos_test` runs with parameters specified inside the code (see **Simulation Settings**).
* The python script `build/python_job_local.py` makes it possible to run the code on multiple threads (CPU cores) by running several times the executable. It requires `mpi4py`, `numpy`, and possibly other python libraries. 

## Simulation Settings:
* Most of settings can be adjusted in `src/include/Settings.hh`. In particular:
  * `std::vector < int > PDG_LIST` : list of particles to be simulated (available are photons, electrons, positrons, positive muons, negative muons, neutrons, protons)
  * `WORLD_MAX_ALT` : the maximum altitude of the simulation world
  * `POTENTIAL_VALUE` : the thunderstorm potential difference (MV)
  * `EFIELD_REGION_Y_LEN` : the altitude extend of the E-field region
  * `EFIELD_REGION_Y_CENTER` : the altitude center of the E-field region
  * `EFIELD_XZ_HALF_SIZE` : the half size of the E-field region in the XZ plane (perpendicular to altitude direction)
  * `longitude`, `latitude` : coordinates were the simulation takes places
  * `year`, `month`, `day` : time when the simulation takes place
  * `CR_GENERATION_ALT_MIN`, `CR_GENERATION_ALT_MAX` : minimum and maximum altitude to generate the initial cosmic ray particles
  * `CR_SAMPLING_XZ_HALF_SIZE` : half size in the XZ plane (perpendicular to altitude) where the initial cosmic rays are generated
  * `RECORD_ALTITUDE` : altitude to record particles.
  * `NB_PARTICLES_TO_GET` : number of particles to get before stopping the simulation.
  * `RECORD_XZ_HALF_SIZE` : half size in the XZ plane of the area where particles are recorded
  * `MAX_POSSIBLE_TIME` : maximum allowed time (particles are killed if their time is higher)

## Other info
* By default, it uses the `G4EmStandardPhysics_option1` physics list (can be changed in `src/src/PhysicsList.cc`) with a change in the "Dr over R value". A maximum step can also be used. Not using one of these two tweaks can lead to incorrect results for large electric fields. See [this article](https://www.geosci-model-dev.net/11/4515/2018/) for more information.
* Each run has its own random number seed that guarantees independent simulations. The `CLHEP::MixMaxRng` random number generator is preferred, as it makes sure sequences of random numbers are independent if two seed numbers are close.
* It is possible to set the simulation to "visualization" mode using the setting `MODE` (=`"run"` or `"visu"`). It is just for debugging and to check that everything is working as expected, because it makes the simulation slow and uses a lot of memory.
* Particles are allowed to be recorded several times if they cross the detection altitude several times. This can be changed in the file `src/src/SteppingAction.cc`.
* Particles are recorded with energies starting from 50 keV. Only energies between 50 keV and 100 MeV are used to build the energy spectrums. Counts above 100 MeV are put in the last energy bin (overflow).