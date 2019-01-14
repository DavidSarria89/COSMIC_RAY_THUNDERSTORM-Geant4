//

// ********************************************************************
#include "StackingAction.hh"
#include "G4RunManager.hh"

// Added to mimic time dependent simulation, and to count the number of
// energetic electron every 1 microseconds
// and turn OFF the electric field if there is a memory overflow

// can be a bad idea to set it on, because it can use a lot of memory compared
// to default G4 behaviour

BaseStackingAction::BaseStackingAction() {

}

BaseStackingAction::~BaseStackingAction() {}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BaseStackingAction::PrepareNewEvent() {

    if (check_RAM_USAGE_LIMIT()) {
        return;
    }

    LIST_ENERGETIC_PART_IDS.clear();
    LIST_ENERGETIC_PART_IDS.reserve(buffer_size);

    VARIABLE_TIME_LIMIT = TIME_STEP;

    NB_EVENT = settings->NB_EVENT;

    settings->current_efield_status = settings->initial_efield_status;

    settings->VARIABLE_TIME_LIMIT_GLOBAL = VARIABLE_TIME_LIMIT;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ClassificationOfNewTrack BaseStackingAction::ClassifyNewTrack(const G4Track *aTrack) {

    // ignore if efield is off
    if (settings->current_efield_status == settings->efield_OFF) return fUrgent;

    // ignore if we have exceeded the ram usage limit
    if (check_RAM_USAGE_LIMIT()) {
        return fUrgent; // default classification
    }

    // Warning : dot not use "aTrack->GetStep()->GetPreStepPoint()", this is not
    // the stepping action, step may not exist

    const G4int ID = aTrack->GetTrackID();

    // if energy > ENER_THRES and if the ID is not already saved

    if (aTrack->GetKineticEnergy() > ENER_THRES) {
        if (does_not_contain(ID, LIST_ENERGETIC_PART_IDS)) {
            LIST_ENERGETIC_PART_IDS.push_back(ID); // save the ID
        }
    }

    if (aTrack->GetGlobalTime() > VARIABLE_TIME_LIMIT) {
        return fWaiting;
    }

    // default classification
    return fUrgent;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BaseStackingAction::print_status() {
    G4cout << "current max time : " << VARIABLE_TIME_LIMIT / microsecond << " microsecond" << G4endl;
    G4cout << ">0.1 MeV electron count : " << LIST_ENERGETIC_PART_IDS.size() << G4endl;
    G4cout << "Event nb : " << NB_EVENT << G4endl;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BaseStackingAction::NewStage() // called when the "Urgent stack" is empty
// and the particles in the "waiting stack" are transfered to the "Urgent stack"
{
    if (check_RAM_USAGE_LIMIT()) return;

    if (settings->current_efield_status == settings->efield_OFF) return;

    if (!LIST_ENERGETIC_PART_IDS.empty()) LIST_ENERGETIC_PART_IDS.clear();


    VARIABLE_TIME_LIMIT += TIME_STEP;

    settings->VARIABLE_TIME_LIMIT_GLOBAL = VARIABLE_TIME_LIMIT;


    if (LIST_ENERGETIC_PART_IDS.size() > 10000) {
        G4cout << "WARNING: more than 10000 particles stacked for 1 single event." << G4endl;
    }

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool BaseStackingAction::check_RAM_USAGE_LIMIT() {
    if (settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED == 1) {
        settings->current_efield_status = settings->efield_OFF;

        if (!LIST_ENERGETIC_PART_IDS.empty()) LIST_ENERGETIC_PART_IDS.clear();

        return true;
    } else {
        return false;
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

bool BaseStackingAction::does_not_contain(const G4int ID, const std::vector<G4int> &LIST_IDS) {

    if (LIST_IDS.empty()) return true;

    return !(std::find(LIST_IDS.begin(), LIST_IDS.end(), ID) != LIST_IDS.end());
}

