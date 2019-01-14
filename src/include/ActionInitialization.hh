#pragma once

#include "EventAction.hh"
#include "G4VUserActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"

class DetectorConstruction;

/// Action initialization class.
///

class ActionInitialization: public G4VUserActionInitialization
{
public:

  explicit ActionInitialization(DetectorConstruction *detector);

  ~ActionInitialization() override;

  void BuildForMaster() const override;

  void Build() const override;

private:

  Settings *settings = Settings::getInstance();

  DetectorConstruction *fDetector;
};
