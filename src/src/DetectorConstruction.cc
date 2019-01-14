// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include <G4GDMLParser.hh>
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "RegionInformation.hh"

using namespace std;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// C     INPUT VARIABLES:
// C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
// C              (Year ignored in current model)
// C        SEC - UT(SEC)
// C        ALT - ALTITUDE(KM)
// C        GLAT - GEODETIC LATITUDE(DEG)
// C        GLONG - GEODETIC LONGITUDE(DEG)
// C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
// C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
// C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
// C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
// C           - ARRAY CONTAINING:
// C             (1) DAILY AP
// C             (2) 3 HR AP INDEX FOR CURRENT TIME
// C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
// C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
// C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
// C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
// C                    TO CURRENT TIME
// C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
// C                    TO CURRENT TIME
// C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
// C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
// C                 MASS 17 IS Anomalous O ONLY.)
// C
// C     NOTES ON INPUT VARIABLES:
// C        UT, Local Time, and Longitude are used independently in the
// C        model and are not of equal importance for every situation.
// C        For the most physically realistic calculation these three
// C        variables should be consistent (STL=SEC/3600+GLONG/15).
// C        The Equation of Time departures from the above formula
// C        for apparent local time can be included if available but
// C        are of minor importance.
// c
// C        F107 and F107A values used to generate the model correspond
// C        to the 10.7 cm radio flux at the actual distance of the Earth
// C        from the Sun rather than the radio flux at 1 AU. The following
// C        site provides both classes of values:
// C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
// C
// C        F107, F107A, and AP effects are neither large nor well
// C        established below 80 km and these parameters should be set to
// C        150., 150., and 4. respectively.
// C
// C     OUTPUT VARIABLES:
// C        D(1) - HE NUMBER DENSITY(CM-3)
// C        D(2) - O NUMBER DENSITY(CM-3)
// C        D(3) - N2 NUMBER DENSITY(CM-3)
// C        D(4) - O2 NUMBER DENSITY(CM-3)
// C        D(5) - AR NUMBER DENSITY(CM-3)
// C        D(6) - TOTAL MASS DENSITY(GM/CM3)
// C        D(7) - H NUMBER DENSITY(CM-3)
// C        D(8) - N NUMBER DENSITY(CM-3)
// C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
// C        T(1) - EXOSPHERIC TEMPERATURE
// C        T(2) - TEMPERATURE AT ALT

// IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T
#include <fortran.hh>

// extrernal fortran subroutine to get MSIS atmospheric densities
extern "C" {
void gtd7_(INTEGER &IYD,   // YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
           REAL &SEC,   // UT(SEC)
           REAL &ALT,   // ALTITUDE(KM)
           REAL &GLAT,  // GEODETIC LATITUDE(DEG)
           REAL &GLONG, // GEODETIC LONGITUDE(DEG)
           REAL &STL,   // LOCAL APPARENT SOLAR TIME
           REAL &F107A, // 81 day AVERAGE OF F10.7 FLUX (centered on day DDD
           REAL &F107,  // DAILY F10.7 FLUX FOR PREVIOUS DAY
           REAL &AP,    // MAGNETIC INDEX(DAILY)
           INTEGER &MASS,  // MASS NUMBER
           REAL *D,
           REAL *T);    // OUTPUT VARIABLES temperatures
}

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction() {
    sWorld = nullptr;
    lWorld = nullptr;
    pWorld = nullptr;

    if (settings->ATMOS_LAYERS_OUTPUT_TO_FILE) {
        asciiFile.open("./tests/alt_dens.txt", std::ios::trunc);

        asciiFile << "ID // altitude (km) // range (km) // density (g/cm2) // has_efield // is_record" << G4endl;
    }

    // sanity check
    if (World_XZ_half_size < settings->CR_SAMPLING_XZ_HALF_SIZE * km) {
        G4cout << "Error in building geometry (DetectorConstruction): CR_SAMPLING_XZ_HALF_SIZE is larger than World_XZ_half_size" << G4endl;
        std::abort();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
    if (settings->ATMOS_LAYERS_OUTPUT_TO_FILE) {
        asciiFile.close();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct() {
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();

    //
    // World Material : Vaccum
    G4NistManager *man = G4NistManager::Instance();
    G4Material *vac = man->FindOrBuildMaterial("G4_Galactic");

    ////////////////////////////////////////
    /// // assigning default region information
    auto *Reginfo_def = new RegionInformation();
    Reginfo_def->Set_World();
    G4Region *defaultRegion = (*(G4RegionStore::GetInstance()))[0]; // the default (world) region is index 0 in the region store
    defaultRegion->SetUserInformation(Reginfo_def);

    ////////////////////////////////////////
    // World Definition
    //
    sWorld = new G4Box("sWorld", World_XZ_half_size * 1.01, World_Y_size * 1.01, World_XZ_half_size * 1.01);
    lWorld = new G4LogicalVolume(sWorld, vac, "lWorld", globalfieldMgr);

    pWorld = new G4PVPlacement(nullptr,
                               G4ThreeVector(),
                               lWorld,
                               "pWorld",
                               nullptr,
                               false,
                               0);

    G4VisAttributes *VisAtt_vac = new G4VisAttributes(G4Colour(0.0, 0.8, 0.0, 0.2));
    lWorld->SetVisAttributes(VisAtt_vac);

    //    VisAtt_vac->SetVisibility(false);

//    build_soil();
    build_air_layers();

    if (settings->USE_RECORD_REGION_STEP_MAX && nb_layers_with_maxStep == 0) {
        G4cout << "ERROR : should have at least one layer with a maximum step set. Aborting" << G4endl;
        std::abort();
    }

    // always return the root volume
    //
    return pWorld;
} // DetectorConstruction::Construct

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Air layers
void DetectorConstruction::build_air_layers() {

    double EF_val = 0.0;

    std::vector<G4Material *> Airs;

    altitude_list = calculate_altitudes_list(0.01 * km, World_Y_size, nb_altitudes);
    Airs = Construct_Atmos_layers_Materials_simple(altitude_list);

    std::vector<G4VSolid *> solid_boxes;
    std::vector<G4LogicalVolume *> logic_boxes;
    std::vector<G4VPhysicalVolume *> phys_boxes;

    G4VisAttributes *visAtt_reg = new G4VisAttributes(G4Colour(0.8, 0.2, 0.8, 0.2));
    visAtt_reg->SetVisibility(true);

    G4VisAttributes *visAtt_reg2 = new G4VisAttributes(G4Colour(0.8, 0.3, 0.8, 0.2));
    visAtt_reg2->SetVisibility(true);

    G4VisAttributes *visAtt_eField_region = new G4VisAttributes(G4Colour(0.49, 0.98, 0., 0.2));
    visAtt_eField_region->SetVisibility(true);

    G4VisAttributes *visAtt_record_region = new G4VisAttributes(G4Colour(0.78, 0.29, 0., 0.2));
    visAtt_eField_region->SetVisibility(true);

    // building altitude layers

    double alt_Min_Efield = (settings->EFIELD_REGION_Y_CENTER - settings->EFIELD_REGION_Y_LEN / 2.000000000000) * km;
    double alt_Max_Efield = (settings->EFIELD_REGION_Y_CENTER + settings->EFIELD_REGION_Y_LEN / 2.000000000000) * km;

    for (unsigned int idx_alt = 0; idx_alt < altitude_list.size() - 1; ++idx_alt) {
        const double delta_alt = (altitude_list[idx_alt + 1] - altitude_list[idx_alt]);

        const double ymin = altitude_list[idx_alt];
        const double ymax = altitude_list[idx_alt + 1];

        const double mean_alt = (ymax + ymin) * 0.50000000000000;

        const double position_Y = mean_alt;

        bool has_efield = false;
        bool has_record = false;

        //            G4cout << Airs[index_air]->GetDensity() / (g / cm3) << G4endl ;

        solid_boxes.push_back(
                new G4Box("s_layer_" + std::to_string(idx_alt), World_XZ_half_size, delta_alt * 0.50000000000000, World_XZ_half_size));

        logic_boxes.push_back(new G4LogicalVolume(solid_boxes.back(), Airs[idx_alt], "l_layer_" + std::to_string(idx_alt)));


        const double densityval = Airs[idx_alt]->GetDensity() / (g / cm3);

        // small change of color to better visualize the layers
        if (idx_alt % 2 == 0) {
            logic_boxes.back()->SetVisAttributes(visAtt_reg);
        } else {
            logic_boxes.back()->SetVisAttributes(visAtt_reg2);
        }

        const double alt_record = settings->RECORD_ALTITUDE * km;

        if (ymin == alt_record) {
            logic_boxes.back()->SetVisAttributes(visAtt_record_region);
            has_record = true;
        }

        if (settings->USE_RECORD_REGION_STEP_MAX) {
            if (std::abs(ymin - alt_record) < 0.25 * km) {
                logic_boxes.back()->SetUserLimits(stepLimit_record);

                nb_layers_with_maxStep++;
            }
        }

        // setting E-field in some regions
        if ((ymin >= alt_Min_Efield) && (ymax <= alt_Max_Efield)) {

            logic_boxes.back()->SetRegion(EFIELD_Region);
            EFIELD_Region->AddRootLogicalVolume(logic_boxes.back());

            logic_boxes.back()->SetVisAttributes(visAtt_eField_region);

            has_efield = true;
            EF_val = (settings->POTENTIAL_VALUE * megavolt) / (settings->EFIELD_REGION_Y_LEN * km);

            logic_boxes.back()->SetFieldManager(localfieldMgr, allLocal);

        } else {
            has_efield = false;
            EF_val = 0.0;
        }

        G4ThreeVector position_box{0, position_Y, 0};

        phys_boxes.push_back(new G4PVPlacement(rotation_null, position_box, logic_boxes.back(),
                                               "atmosLayer_" + std::to_string(idx_alt) + "_" + std::to_string(mean_alt / km),
                                               lWorld, false, 0, false));


        if (settings->ATMOS_LAYERS_OUTPUT_TO_FILE) {
            asciiFile << idx_alt << " " << ymin / km << " " << ymax / km << " " << delta_alt * 0.50 / km << " " << densityval << " " << has_efield
                      << " " << EF_val << " " <<
                      has_record << G4endl;
        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<double> DetectorConstruction::calculate_altitudes_list(const double alt_min, const double alt_max_construction, const uint nb_altis)

// fills the vector altitude_list
{
    // defining the altitude vector
    for (uint idx_alt = 0; idx_alt < nb_altis; idx_alt++) // geocentric altitude_list
    {
        double div = nb_altis - 1.0;
        double idx_alt_ld = idx_alt;
        double alt = std::exp(std::log(alt_min) + (std::log(alt_max_construction) - std::log(alt_min)) * idx_alt_ld / div);
        altitude_list.push_back(alt);

        if (std::isnan(altitude_list.back())) {
            G4cout << "Error, one altitude is NaN" << G4endl;
            std::abort();
        }
    }

    // adding altitude_list of E-field region

    if (!contains(altitude_list, EFIELD_alt_Mid)) {
        altitude_list.push_back(EFIELD_alt_Mid);
    }
    double to_add = EFIELD_alt_Min - 0.01 * km;

    if (!contains(altitude_list, to_add)) {
        altitude_list.push_back(to_add);
    }

    to_add = EFIELD_alt_Max + 0.01 * km;

    if (!contains(altitude_list, to_add)) {
        altitude_list.push_back(to_add);
    }

    if (!contains(altitude_list, EFIELD_alt_Min)) {
        altitude_list.push_back(EFIELD_alt_Min);
    }

    if (!contains(altitude_list, EFIELD_alt_Max)) {
        altitude_list.push_back(EFIELD_alt_Max);
    }

    // adding altitude_list of record region (narrow slice)


    to_add = settings->RECORD_ALTITUDE * km;

    if (!contains(altitude_list, to_add)) {
        altitude_list.push_back(to_add);
    }

    to_add = settings->RECORD_ALTITUDE * km - 1.0 * km;

    if (!contains(altitude_list, to_add)) {
        altitude_list.push_back(to_add);
    }


    to_add = settings->RECORD_ALTITUDE * km + 1.0 * km;

    if (!contains(altitude_list, to_add)) {
        altitude_list.push_back(to_add);
    }

    // check if there are duplicates on the list (that should not happen), if so abort

    if (hasDuplicates(altitude_list)) {
        G4cout << "ERROR : There are duplicates values in the altitude list. Aborting." << G4endl;
        std::abort();
    }

    // sorting in increasing value
    std::sort(altitude_list.begin(), altitude_list.end());

    if (!is_increasing(altitude_list)) {
        G4cout << "ERROR: altitude list should be strictly increasing. Aborting." << G4endl;
        std::abort();
    }

    G4cout << "Successfully created altitude list" << G4endl;

    return altitude_list;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool DetectorConstruction::hasDuplicates(const std::vector<double> &arr) {
    for (uint i = 0; i < arr.size(); ++i) {
        for (uint j = i + 1; j < arr.size(); ++j) {
            if (arr[i] == arr[j]) {
                return true;
            }
        }
    }

    return false;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool DetectorConstruction::is_increasing(const std::vector<double> &arr) {
    for (uint i = 0; i < arr.size() - 1; ++i) {
        if (arr[i + 1] < arr[i]) return false;
    }

    return true;
}
// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Caltulating materials of the atmopsheric layers, based on the MSIS C++ model integrated to this code
// ref : https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html
// simplified version using only the total mass density and same ratio of N2 O2 and AR
// Caltulating materials of the atmopsheric layers, based on the MSIS C++ model integrated to this code
// ref : https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html
// simplified version using only the total mass density and same ratio of N2 O2
std::vector<G4Material *> DetectorConstruction::Construct_Atmos_layers_Materials_simple(const std::vector<double> &altitudes) {
    std::vector<G4Material *> Airs;

    // Vaccum
    G4NistManager *man = G4NistManager::Instance();
    G4Material *vaccum = man->FindOrBuildMaterial("G4_Galactic");

    //    elHe = new G4Element(name = "Helium", symbol = "He", z = 2., He_molarMass);
    //    elH  = new G4Element(name = "Hydrogen", symbol = "H", z = 1., H_molarMass);

    for (uint idx_alt = 0; idx_alt < altitudes.size() - 1; idx_alt++) {
        const double innerAlt = altitudes[idx_alt];
        const double outerAlt = altitudes[idx_alt + 1];
        const double mid_altitude_in_km = (innerAlt + outerAlt) / 2.0 / km; // geocentric altitude

        if (mid_altitude_in_km > World_Y_size / km) {
            Airs.push_back(vaccum);
        } else {
            INTEGER input_iyd = 172; // IYD - YEAR AND DAY AS YYDDD
            REAL input_sec = 29000.0;
            REAL input_alt = (REAL) mid_altitude_in_km;
            REAL input_g_lat = (REAL) settings->latitude;
            REAL input_g_long = (REAL) settings->longitude;
            REAL input_lst = 16.0;
            REAL input_f107A = 150.0;
            REAL input_f107 = 150.0;
            REAL input_ap = 4.0;
            INTEGER input_mass = 48;
            REAL output_D[9];
            REAL output_T[2];

            gtd7_(input_iyd, input_sec, input_alt, input_g_lat, input_g_long, input_lst, input_f107A, input_f107, input_ap, input_mass, output_D, output_T); //
            // MSIS,// call
            if (std::isnan(output_D[5]) || std::isinf(isnan(output_D[5]))) {
                G4cout << "ERROR : density from gtd7_ is NaN. Aborting" << G4endl;
                std::abort();
            }

            double density_air = output_D[5] * g / cm3; // getting density and converting it to the GEANT4 system of unit

            if (MAKE_EVERYTHING_VACCUM) {
                density_air = 1e-24 * (g / cm3);
            }

            G4String name;
            G4int ncomponents;
            double fractionmass;

            Airs.push_back(new G4Material(name = "Air_" + std::to_string(idx_alt), density_air, ncomponents = 2));

            Airs.back()->AddElement(elN, fractionmass = 0.7615);
            Airs.back()->AddElement(elO, fractionmass = 0.2385);
        }
    }

    G4cout << "Successfully created air materials list" << G4endl;
    return Airs;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Returns interpolated value at x from parallel arrays ( xData, yData )
//   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
//   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
double DetectorConstruction::interpolate(const vector<double> &xData, const vector<double> &yData, const double &x, const bool extrapolate) {
    int size = static_cast<int>(xData.size());

    int i = 0;                // find left end of interval for interpolation

    if (x >= xData[size - 2]) // special case: beyond right end
    {
        i = size - 2;
    } else {
        while (x > xData[i + 1]) i++;
    }

    double xL = xData[i], yL = yData[i], xR = xData[i + 1], yR = yData[i + 1]; // points on either side (unless beyond ends)

    if (!extrapolate)                                                          // if beyond ends of array and not extrapolating
    {
        if (x < xL) {
            yR = yL;
        }

        if (x > xR) {
            yL = yR;
        }
    }

    double dydx = (yR - yL) / (xR - xL); // gradient

    return yL + dydx * (x - xL);         // linear interpolation
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DetectorConstruction::contains(const vector<double> &v, const double &x) {
    return std::find(v.begin(), v.end(), x) != v.end();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double DetectorConstruction::get_scale(const double &alt) {
    // returns the atmospheric relative scale compared to sea level (>1)
    std::vector<double> alt_list =
            {0, 0.500000000000000, 1.000000000000000, 1.500000000000000,
             2.000000000000000,
             2.500000000000000,
             3.000000000000000, 3.500000000000000,
             4.000000000000000, 4.500000000000000, 5.000000000000000,
             5.500000000000000, 6.000000000000000, 6.500000000000000, 7.000000000000000,
             7.500000000000000,
             8.000000000000000,
             8.500000000000000, 9.000000000000000, 9.500000000000000, 10.000000000000000,
             10.500000000000000, 11.000000000000000, 11.500000000000000, 12.000000000000000,
             12.500000000000000,
             13.000000000000000,
             13.500000000000000, 14.000000000000000, 14.500000000000000, 15.000000000000000,
             15.500000000000000, 16.000000000000000, 16.500000000000000, 17.000000000000000,
             17.500000000000000,
             18.000000000000000,
             18.500000000000000, 19.000000000000000, 19.500000000000000, 20.000000000000000,};

    std::vector<double> scale_list =
            {1.000000000000000, 1.059301380991064, 1.121238177128117, 1.184377838328792, 1.249042145593870, 1.317304778260430,
             1.388415672913118, 1.463688404983724,
             1.543377914546100, 1.628371628371629,
             1.719862833025587, 1.818435364663227, 1.925291598996014, 2.041647095663066, 2.168634624979211, 2.307556184746062,
             2.460377358490566,
             2.628502318081032, 2.813376483279396, 3.018518518518519,
             3.245395719263315, 3.496916063287745, 3.774240231548481, 4.080100125156446, 4.414353419092755, 4.781811514484781,
             5.180770758839889,
             5.613430908308223, 6.082089552238807, 6.589186457806973,
             7.133479212253830, 7.720544701006513, 8.348271446862997, 9.024221453287199, 9.753178758414361, 10.541632983023444,
             11.398601398601400,
             12.325141776937620, 13.334696799263730, 14.431164231961047,
             15.627996164908916,};

    if ((alt > 20.) || (alt < 0.)) {
        G4cout << "ERROR in get_scale : altitude is not between 0 and 20 km. Aborting.";
        std::abort();
    }

    return interpolate(alt_list, scale_list, alt, false);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::generate_Soil_material() {
    G4NistManager *man = G4NistManager::Instance();
    G4Material *SILICON_DIOXIDE = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material *ALUMINUM_OXIDE = man->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    G4Material *FERRIC_OXIDE = man->FindOrBuildMaterial("G4_FERRIC_OXIDE");
    G4Material *CALCIUM_OXIDE = man->FindOrBuildMaterial("G4_CALCIUM_OXIDE");
    G4Material *MAGNESIUM_OXIDE = man->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE");
    G4Material *TITANIUM_DIOXIDE = man->FindOrBuildMaterial("G4_TITANIUM_DIOXIDE");
    //    G4_SILICON_DIOXIDE 61.3 %
    //    G4_ALUMINUM_OXIDE 13 %
    //    G4_FERRIC_OXIDE 2.5 %
    //    G4_CALCIUM_OXIDE 1.6 %
    //    G4_MAGNESIUM_OXIDE 0.7 %
    //    G4_TITANIUM_DIOXIDE 0.6 %
    //    organic material : 20.3 %; made of
    //    C 50  %
    //    O 42  %
    //    H 5 %
    //    N 3 %
    G4Element *elH = man->FindOrBuildElement("H");
    G4Element *elC = man->FindOrBuildElement("C");
    G4Element *elOO = man->FindOrBuildElement("O");
    G4Element *elNN = man->FindOrBuildElement("N");
    const double density = 0.8 * g / cm3;
    G4String name;
    G4int ncomponents, natoms;
    G4Material *organic_mat = new G4Material(name = "organic_mat", density, ncomponents = 4);
    organic_mat->AddElement(elC, natoms = 50);
    organic_mat->AddElement(elOO, natoms = 42);
    organic_mat->AddElement(elH, natoms = 5);
    organic_mat->AddElement(elNN, natoms = 3);
    double fraction;
    const double soil_density = 1.4 * g / cm3;
    SOIL = new G4Material(name = "Soil", soil_density, ncomponents = 7);
    SOIL->AddMaterial(SILICON_DIOXIDE, fraction = 61.3 * perCent);
    SOIL->AddMaterial(ALUMINUM_OXIDE, fraction = 13. * perCent);
    SOIL->AddMaterial(FERRIC_OXIDE, fraction = 2.5 * perCent);
    SOIL->AddMaterial(CALCIUM_OXIDE, fraction = 1.6 * perCent);
    SOIL->AddMaterial(MAGNESIUM_OXIDE, fraction = 0.7 * perCent);
    SOIL->AddMaterial(TITANIUM_DIOXIDE, fraction = 0.6 * perCent);
    SOIL->AddMaterial(organic_mat, fraction = 20.3 * perCent);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::build_soil() {
    generate_Soil_material();

    G4VisAttributes *visAtt_soil = new G4VisAttributes(G4Colour(0.6250, 0.3203, 0.1758, 0.2));
    visAtt_soil->SetVisibility(true);

    G4Box *Soil_solid = new G4Box("Soil_solid", World_XZ_half_size, 2.5 * km, World_XZ_half_size);

    G4LogicalVolume *Soil_logic = new G4LogicalVolume(Soil_solid, SOIL, "Soil_logic");

    Soil_logic->SetVisAttributes(visAtt_soil);

    G4VPhysicalVolume *Soil_phys = new G4PVPlacement(rotation_null, G4ThreeVector{0, (settings->SOIL_ALT_MAX * km - 2.5 * km), 0},
                                                     Soil_logic,
                                                     "Soil_phys",
                                                     lWorld, false, 0, true);

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField() {


    if (!fEmFieldSetup.Get()) {
        FieldSetup *emFieldSetup = new FieldSetup();

        fEmFieldSetup.Put(emFieldSetup);
        G4AutoDelete::Register(emFieldSetup); //Kernel will delete the messenger
    }

    bool allLocall = true;
    lWorld->SetFieldManager(fEmFieldSetup.Get()->GetLocalFieldManager(),
                            allLocall);
    // step limitation only for charged particles, initialization value, will be changed in the SteppingAction
//    G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(
//            10.0 * m);

}
