#include <G4ThreeVector.hh>
#include "cosmic_ray_generator_CRY.hh"
#include "RNGWrapper.hh"
#include "Randomize.hh"


Cosmic_Ray_Generator_CRY::Cosmic_Ray_Generator_CRY(const char *inputfile) {

    std::ifstream inputFile;
    inputFile.open(inputfile, std::ios::in);
    char buffer[1000];

    if (inputFile.fail()) {
        if (*inputfile != 0) //....only complain if a filename was given
            G4cout << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << G4endl;

        InputState = -1;
    } else {
        std::string setupString("");

        while (!inputFile.getline(buffer, 1000).eof()) {
            setupString.append(buffer);
            setupString.append(" ");
        }

        CRYSetup *setup = new CRYSetup(setupString, "./CRY_data/");

        gen = new CRYGenerator(setup);

        // set random number generator
        RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);

        setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
        InputState = 0;
    }

    // create a vector to store the CRY particle properties
    vect = new std::vector<CRYParticle *>;

    // Create the table containing all particle names
    particleTable = G4ParticleTable::GetParticleTable();

}

//----------------------------------------------------------------------------//
void Cosmic_Ray_Generator_CRY::InputCRY() {
    InputState = 1;
}

//----------------------------------------------------------------------------//
void Cosmic_Ray_Generator_CRY::UpdateCRY(std::string *MessInput) {
    CRYSetup *setup = new CRYSetup(*MessInput, "../data");

    gen = new CRYGenerator(setup);

    // set random number generator
    RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
    setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
    InputState = 0;

}

void Cosmic_Ray_Generator_CRY::CRYFromFile(G4String newValue) {
    // Read the cry input file
    std::ifstream inputFile;
    inputFile.open(newValue, std::ios::in);
    char buffer[1000];

    if (inputFile.fail()) {
        G4cout << "Failed to open input file " << newValue << G4endl;
        G4cout << "Make sure to define the cry library on the command line" << G4endl;
        InputState = -1;
    } else {
        std::string setupString("");

        while (!inputFile.getline(buffer, 1000).eof()) {
            setupString.append(buffer);
            setupString.append(" ");
        }

        CRYSetup *setup = new CRYSetup(setupString, "../data");

        gen = new CRYGenerator(setup);

        // set random number generator

        RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);


        setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
        InputState = 0;
    }
}



// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

geant4_initial_cosmic_ray Cosmic_Ray_Generator_CRY::generate_Cosmic_ray() {

    std::vector<geant4_initial_cosmic_ray> g4_cosmic_list = {};

    if (InputState != 0) {
        G4String *str = new G4String("CRY library was not successfully initialized");
        //G4Exception(*str);
        G4Exception("PrimaryGeneratorAction", "1",
                    RunMustBeAborted, *str);
    }

//    G4String particleName;
    vect->clear();
    gen->genEvent(vect);

    int jj = 0;

    const double CRY_max_possible_XY_size = 300.0 * meter;
    const double scale_factor = (settings->CR_SAMPLING_XZ_HALF_SIZE * km) / CRY_max_possible_XY_size * 2.0;

    G4ThreeVector mom_dir = G4ThreeVector((*vect)[jj]->u(), (*vect)[jj]->w(), (*vect)[jj]->v());
    G4ThreeVector position = G4ThreeVector((*vect)[jj]->x() * m * scale_factor, (*vect)[jj]->z() * m + altitude_generator, (*vect)[jj]->y() * m * scale_factor);
    double energy = (*vect)[jj]->ke() * MeV;
    double time = 0.0;

    G4String particleName = CRYUtils::partName((*vect)[jj]->id());

    if (particleName == "electron") {
        particleName = "e-";
    } else if (particleName == "muon") {
        particleName = "mu-";
    }

    G4ParticleDefinition *particle = particleTable->FindParticle(particleName);

//    int PDG_nb = particle->GetPDGEncoding();

    geant4_initial_cosmic_ray g4_cosmic = {};

    g4_cosmic.energy = energy;
    g4_cosmic.time = time;
    g4_cosmic.momentum_ini = mom_dir;
    g4_cosmic.position_ini = position;
    g4_cosmic.g4_particle = particle;

    return g4_cosmic;

}