#include "Settings.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Settings *Settings::instance = nullptr;

Settings * Settings::getInstance()
{
  if (instance == nullptr)
  {
    instance = new Settings;
  }

  return instance;
}
