#pragma once

#include "G4SDManager.hh"
#include <G4AssemblyVolume.hh>
#include <G4String.hh>
#include <G4Types.hh>
#include <G4VUserDetectorConstruction.hh>
#include "G4VPhysicalVolume.hh"

#include <CLHEP/Units/SystemOfUnits.h>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

// #include "E_Field.hh"

#include "Settings.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4EqMagElectricField.hh"
#include "G4ElectroMagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4FieldManager.hh"
#include "G4ElectroMagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4ClassicalRK4.hh"
#include "G4UniformElectricField.hh"
#include "FieldSetup.hh"
#include "G4UserLimits.hh"
#include "G4DormandPrince745.hh"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "G4MagIntegratorDriver.hh"
#include <vector>

#include "G4Cache.hh"
#include "G4AutoDelete.hh"
// #include "E_Field.hh"

extern "C" {
#include <coordinates_conversions.h>
}
extern "C" {
#include <nrlmsise-00.h>
}

class G4LogicalVolume;

class G4Material;

class G4Box;

class G4Cons;

class G4Tubs;

class G4VPhysicalVolume;

class G4BooleanSolid;

class G4UnionSolid;

class G4UserLimits;


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction {
public:

    DetectorConstruction();

    ~DetectorConstruction() override;

    void ConstructSDandField() override;

public:

    G4VPhysicalVolume *Construct() override;

private:

    G4Cache<FieldSetup *> fEmFieldSetup;

    void generate_Soil_material();

    bool allLocal = true;

    Settings *settings = Settings::getInstance();

    std::vector<double> altitude_list;

    // usefull null position and rotation
    G4RotationMatrix *rotation_null = new G4RotationMatrix;
    G4ThreeVector translation_null = G4ThreeVector(0, 0, 0);

    // World volume
    const double World_XZ_half_size = 50. * CLHEP::km;
    const double World_Y_size = settings->WORLD_MAX_ALT * CLHEP::km;
    int nb_altitudes = 100;
    G4Box *sWorld;           // pointer to the solid envelope
    G4LogicalVolume *lWorld; // pointer to the logical envelope
    G4VPhysicalVolume *pWorld; // pointer to the physical envelope
    std::vector<G4Material *> Construct_Atmos_layers_Materials_simple(const std::vector<double> &altitudes);

    G4FieldManager *globalfieldMgr = nullptr;
    G4FieldManager *localfieldMgr = nullptr;
    double field_val_temp = 0.0; // just for debug

    bool MAKE_EVERYTHING_VACCUM = false; // flag to put all density to very low value, to be able to check/debug cosmic ray sampling

    double interpolate(const std::vector<double> &xData,
                       const std::vector<double> &yData,
                       const double &x,
                       const bool extrapolate);

    //        std::vector < G4Material * > Construct_Atmos_layers_Materials(const std::vector < double > altitudes_);
    std::vector<double> calculate_altitudes_list(const double alt_min, const double alt_max_construction, const uint nb_altis);

    void build_air_layers();

    void build_soil();


    G4Region *EFIELD_Region = new G4Region("EFIELD");

    std::ofstream asciiFile;

    G4bool contains(const std::vector<double> &v,
                    const double &x);

    bool hasDuplicates(const std::vector<double> &arr);

    bool is_increasing(const std::vector<double> &arr);

    double get_scale(const double &alt);

    G4Element *elN = new G4Element("Nitrogen",
                                   "N",
                                   7.,
                                   14.01 * g / mole);
    G4Element *elO = new G4Element("Oxygen",
                                   "O",
                                   8.,
                                   16.00 * g / mole);

    G4Material *SOIL = nullptr;

    double EFIELD_alt_Min = (settings->EFIELD_REGION_Y_CENTER - settings->EFIELD_REGION_Y_LEN / 2.00000000) * km;
    double EFIELD_alt_Max = (settings->EFIELD_REGION_Y_CENTER + settings->EFIELD_REGION_Y_LEN / 2.00000000) * km;
    double EFIELD_alt_Mid = settings->EFIELD_REGION_Y_CENTER * km;

    G4UserLimits *stepLimit_record = new G4UserLimits(settings->STEP_MAX_VAL_RECORD_REGION);

    uint nb_layers_with_maxStep = 0;
};

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace datetools {
    namespace details {
        constexpr unsigned int days_to_month[2][12] =
                {
                        // non-leap year
                        {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334},
                        // leap year
                        {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335},
                };
    }

    constexpr bool is_leap(int const year) noexcept {
        return year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    }

    constexpr unsigned int day_of_year(int const year, unsigned int const month, unsigned int const day) {
        return details::days_to_month[is_leap(year)][month - 1] + day;
    }
}