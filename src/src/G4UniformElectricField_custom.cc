#include "G4SystemOfUnits.hh"
#include "G4UniformElectricField_custom.hh"


// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UniformElectricField_custom::G4UniformElectricField_custom(const G4ThreeVector& x_y_z_values) {
    EFIELD_x = x_y_z_values[0];
    EFIELD_y = x_y_z_values[1];
    EFIELD_z = x_y_z_values[2];

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UniformElectricField_custom::~G4UniformElectricField_custom() {}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// global E-field is made with a component from the plane (1/r decrease) and a component with a constant E-field

void G4UniformElectricField_custom::GetFieldValue(const G4double Point[4], G4double *field) const {

    // field is really field[6]: Bx,By,Bz,Ex,Ey,Ez
    // Point[0],Point[1],Point[2] are x-, y-, z-cordinates, Point[3] is time

    field[0] = field[1] = field[2] = field[3] = field[4] = field[5] = 0.0;

    // protect against Geant4 bug that calls us with point[] NaN.
    if (Point[0] != Point[0]) return;

//    if (Point[3] > time_limit) return;

    if (settings->current_efield_status == settings->efield_OFF) return;

    if ((Point[1] > EFIELD_alt_min_mm) && (Point[1] < EFIELD_alt_max_mm)
        && (std::abs(Point[0]) < EFIELD_XY_HALF_SIZE_mm)
        && (std::abs(Point[2]) < EFIELD_XY_HALF_SIZE_mm)) {

        field[3] = EFIELD_x;
        field[4] = EFIELD_y;
        field[5] = EFIELD_z;
        return;

    } else {

        return;
    }

}