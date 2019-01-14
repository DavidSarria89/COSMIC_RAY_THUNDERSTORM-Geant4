// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "Run.hh"
#include "myUtils.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction() {

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *) {
    settings->NB_EVENT++;
    settings->event_lead_to_detection = {false, false, false, false, false};

#ifndef NDEBUG // debug mode

    if (settings->NB_EVENT % print_nb == 0) {
        G4cout << "Begin of event : " << settings->NB_EVENT << G4endl;
    }

#endif // ifndef NDEBUG

    settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED = 0;

    analysis->record_list.clear();
    analysis->record_list.reserve(100);

    settings->current_efield_status = settings->initial_efield_status;

    settings->T0 = myUtils::get_wall_time() / 1.0e6;

    settings->EVENT_DURATION_in_CPU_TICKS = 0;

    settings->CPU_TIME_LIMIT_PER_EVENT_HAS_BEEN_REACHED = false;

} // EventAction::BeginOfEventAction

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *) {

    for (uint ii = 0; ii < settings->event_lead_to_detection.size(); ++ii) {
        if (settings->event_lead_to_detection[ii]) {
            settings->NB_EVENTS_WITH_DETECTION[ii]++;
        }
    }

} // EventAction::EndOfEventAction