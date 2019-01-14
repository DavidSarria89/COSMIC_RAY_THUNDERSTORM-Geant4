#pragma once

#include <vector>

#include "globals.hh"

#include "G4ElectroMagneticField.hh"
#include "Settings.hh"


class G4UniformElectricField_custom : public G4ElectroMagneticField {
public:

    G4UniformElectricField_custom(const G4ThreeVector& x_y_z_values);

    ~G4UniformElectricField_custom();

    /// DoesFieldChangeEnergy() returns true.
    virtual G4bool DoesFieldChangeEnergy() const {
        return true;
    }

    void GetFieldValue(const G4double Point[4],
                       G4double *field) const;

private:

    Settings *settings = Settings::getInstance();

    const double total_potential_diff = settings->POTENTIAL_VALUE * megavolt;

    const double acceleration_region_length = settings->EFIELD_REGION_Y_LEN * km;

    double EFIELD_mag = 0.;

    const double time_limit = settings->MAX_POSSIBLE_TIME;

    const double EFIELD_XY_HALF_SIZE_mm = settings->EFIELD_XZ_HALF_SIZE * km;

    const double EFIELD_alt_min_mm = settings->EFIELD_REGION_Y_CENTER * km - settings->EFIELD_REGION_Y_LEN / 2.0 * km;
    const double EFIELD_alt_max_mm = settings->EFIELD_REGION_Y_CENTER * km + settings->EFIELD_REGION_Y_LEN / 2.0 * km;

    const double density_sea_level = 1.340E-03 * g / cm3;

    double EFIELD_x, EFIELD_y, EFIELD_z;

};