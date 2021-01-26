#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "Run.hh"

extern double INITIAL_ENERGY_MEV;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction *det, EventAction *event) : G4UserSteppingAction(), fDetector(det), fEventAction(event) {

    INITIAL_TIME = myUtils::get_wall_time() / 1.0e6;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *aStep) {
    //////// AVOID PUTTING THIS LINE, it will produce incorrect energy record
    //        if (aStep->GetTrack()->GetTrackStatus() != fAlive) return;
    ////////////////////////////////////////////

    theTrack = aStep->GetTrack();
    thePrePoint = aStep->GetPreStepPoint();
    thePostPoint = aStep->GetPostStepPoint();

//    if (settings->CPU_TIME_LIMIT_PER_EVENT_HAS_BEEN_REACHED) {
//        theTrack->SetTrackStatus(fStopAndKill);
//        return;
//    }

    // check RAM usage

    CHECK_COUNTER++;

    // usage of "CPU ticks" is preferred to have consistent calculation amount between runs and usage of difference computers or load

    if (CHECK_COUNTER > RAM_CHECK_COUNTER_MAX) {

        CHECK_COUNTER = 0;

        settings->EVENT_DURATION_in_CPU_TICKS++;

        ALL_TICKS_COUNTER++;

        // check event CPU time duration
        if (!settings->CPU_TIME_LIMIT_PER_EVENT_HAS_BEEN_REACHED) {

            if (settings->EVENT_DURATION_in_CPU_TICKS > settings->MAX_CPU_TICKS) {

                G4cout << G4endl << G4endl << "Run ID (rng seed): " << settings->RANDOM_SEED << G4endl;

                G4cout << "WARNING: event seems to take too long in terms of computation time. Turning OFF electric field." << G4endl;

                settings->CPU_TIME_LIMIT_PER_EVENT_HAS_BEEN_REACHED = true;

                settings->CPU_TIME_LIMIT_PER_EVENT_HAS_BEEN_REACHED_ONCE = true;

                settings->current_efield_status = settings->efield_OFF;

                settings->NB_LONG_CPU_TIME++;

            }
        }

        // check RAM usage
        WT1 = myUtils::get_wall_time() / 1.0e6; // seconds

        if (std::abs(WT1 - WT2) > TIME_PRINT_SEC) {

            WT2 = WT1;

            const double USED_RAM = myUtils::check_USED_RAM(); // in GB

            const double total_ram = myUtils::check_available_RAM(); // in GB

            if (USED_RAM > 5.0) {
                G4cout << G4endl << G4endl << "WARNING: used RAM is more than 5 GB" << G4endl << G4endl;
            }

            if (USED_RAM > total_ram * settings->RAM_FRACTION_LIMIT) {
                G4cout << G4endl << G4endl << "WARNING: used RAM is more than " + std::to_string(total_ram * settings->RAM_FRACTION_LIMIT) + " GB."
                       << G4endl << G4endl;
//                    settings->current_efield_status = settings->efield_OFF;
                settings->RAM_USAGE_LIMIT_HAS_BEEN_REACHED = 1;
//                    settings->CHECK_RAM = false;
            }

            const double dt_since_run_start = (myUtils::get_wall_time() / 1.0e6 - INITIAL_TIME);

            const double dt_since_event_start = (myUtils::get_wall_time() / 1.0e6 - settings->T0);

            G4cout << G4endl << G4endl << "Run ID (rng seed): " << settings->RANDOM_SEED << G4endl;
            G4cout << "Use stacking action: " << settings->USE_STACKING_ACTION << G4endl;
            G4cout << "Initial sampled type (PARMA): " << get_name(settings->INITIAL_SAMPLE_TYPE) << G4endl;
            G4cout << "Energy of current initial particle: " << std::round(INITIAL_ENERGY_MEV * 100) / 100 << " MeV" << G4endl;
            G4cout << "Currently used RAM: " << std::round(USED_RAM * 1000) / 1000 << " GB" << G4endl;
            G4cout << "Available RAM: " << std::round(total_ram * 1000) / 1000 << " GB" << G4endl;
            G4cout << "Average CPU time to compute 10 events: " << std::round(dt_since_run_start / double(settings->NB_EVENT) * 10.0 * 1000.0) / 1000.0
                   << " seconds" << G4endl;
            G4cout << "Average CPU time to get 10 records: " << std::round(dt_since_run_start / double(analysis->NB_OUTPUT) * 10.0 * 1000.0) / 1000.0
                   << " seconds" << G4endl;
            G4cout << "Potential: " << settings->POTENTIAL_VALUE << G4endl;
            G4cout << "Number of events: " << settings->NB_EVENT << G4endl;
            G4cout << "Number of outputs: " << analysis->NB_OUTPUT << G4endl;
            G4cout << "Ratio: " << double(analysis->NB_OUTPUT) / double(settings->NB_EVENT) << G4endl;
            G4cout << "Current event duration: " << std::round(dt_since_event_start * 1000.0) / 1000.0 << " seconds" << G4endl;
            G4cout << "Current event duration: " << settings->EVENT_DURATION_in_CPU_TICKS << " tick(s)" << G4endl;
            G4cout << "Real duration of a CPU tick (seconds), whole run: " << std::round(dt_since_run_start / double(ALL_TICKS_COUNTER) * 1000.0) / 1000.0
                   << G4endl;

            G4cout << "Real duration of a CPU tick (seconds), this event only: "
                   << std::round(dt_since_event_start / double(settings->EVENT_DURATION_in_CPU_TICKS) * 1000.0) / 1000.0
                   << G4endl;

            G4cout << "E-field status: " << !settings->current_efield_status << G4endl;

        }
    }


    if (thePrePoint->GetGlobalTime() > settings->MAX_POSSIBLE_TIME) {
#ifndef NDEBUG
        G4cout << G4endl << "PARTICLE IS TERMINATED BECAUSE IT EXCEEDED THE TIME LIMIT. " << thePrePoint->GetGlobalTime() / us << G4endl << G4endl;
#endif // ifndef NDEBUG
        theTrack->SetTrackStatus(fStopAndKill);
        return;
    }


    if (thePrePoint->GetPosition().y() < double(settings->SOIL_ALT_MAX * km - 0.5 * km)) {
        theTrack->SetTrackStatus(fStopAndKill);
        return;
    }

    if (thePrePoint->GetPosition().y() > settings->WORLD_MAX_ALT * km) {
        theTrack->SetTrackStatus(fStopAndKill);
        return;
    }

    if (std::abs(thePrePoint->GetPosition().x()) > settings->RECORD_XZ_HALF_SIZE * km) {
        return;
    }

    if (std::abs(thePrePoint->GetPosition().z()) > settings->RECORD_XZ_HALF_SIZE * km) {
        return;
    }

    ///////////////////////
    const int PDG_num = theTrack->GetParticleDefinition()->GetPDGEncoding();
    ///////////////////////

//    if (PDG_num == 22 && thePrePoint->GetKineticEnergy() > 10.0 * GeV) {
//        G4cout << G4endl << "KILLING PARTICLE BECAUSE ENERGY IS TOO HIGH. " << G4endl << G4endl;
//        theTrack->SetTrackStatus(fStopAndKill);
//        return;
//    }

    if (PDG_num == -12 || PDG_num == 12) { // killing if neutrino
        theTrack->SetTrackStatus(fStopAndKill);
        return;
    }

    // cleaning particles below 8 keV to improve performance
    if (thePrePoint) {
        // avoid killing positrons below low energy threshold to make sure we still have annihilation

        if (thePrePoint->GetKineticEnergy() < settings->LOW_ENERGY_THRES) {
            if (PDG_num != settings->pdg_posi) {
                theTrack->SetTrackStatus(fStopAndKill);
                return;
            }
            // using fStopAndAlive for positrons makes a bug that saturates the RAM...
        }
    }

    // skipping the rest if current particle is not in list of wanted particles
    index_found indexFound = find_particle_index(PDG_num);

    if (!indexFound.found) {
        return;
    }

    /// stacking action check, if pseudo time oriented simulation is activated (slow)
    if (settings->USE_STACKING_ACTION) {
        if (thePrePoint) {
            if (settings->current_efield_status == settings->efield_ON) {
                if (theTrack->GetGlobalTime() > settings->VARIABLE_TIME_LIMIT_GLOBAL) {
                    theTrack->SetTrackStatus(fSuspend);
                    return;
                }
            }
        }
    }

    // check if particle should be recorded

    if (thePrePoint->GetStepStatus() == fGeomBoundary || thePostPoint->GetStepStatus() == fGeomBoundary) {

        const G4String &vol_name = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();

        if (vol_name.find("record_phys") != std::string::npos && thePrePoint->GetPosition().y() == settings->RECORD_ALTITUDE * km) {

            const double ener = thePrePoint->GetKineticEnergy() * MeV;

            const double momy = thePrePoint->GetMomentumDirection().y();

            // WARNING : PARTICLES ARE ALLOWED TO BE RECORED 2 TIMES if they cross the detection altitude with upwards and downwards momentum

            if ((ener >= settings->ENERGY_MIN_RECORD) && (ener <= settings->ENERGY_MAX_RECORD)) {

                if (analysis->check_record(theTrack->GetTrackID(), momy)) {

                    analysis->add_NB_OUTPUT();

                    uint i_part = indexFound.index;

                    analysis->fill_histogram_E(i_part, ener);
                    analysis->fill_histogram_mX(i_part, thePrePoint->GetMomentumDirection().x());
                    analysis->fill_histogram_mY(i_part, thePrePoint->GetMomentumDirection().y());
                    analysis->fill_histogram_mZ(i_part, thePrePoint->GetMomentumDirection().z());

                    analysis->counter_total[i_part]++;

                    settings->event_lead_to_detection[i_part] = true;

                    if (momy > 0.0) {
                        analysis->counter_up[i_part]++;
                    } else {
                        analysis->counter_down[i_part]++;
                    }
                }
            } else if (ener >= settings->ENERGY_MAX_RECORD) {

                if (analysis->check_record(theTrack->GetTrackID(), momy)) {

                    analysis->add_NB_OUTPUT();

                    uint i_part = indexFound.index;

                    analysis->fill_histogram_E(i_part, ener);
                    analysis->fill_histogram_mX(i_part, thePrePoint->GetMomentumDirection().x());
                    analysis->fill_histogram_mY(i_part, thePrePoint->GetMomentumDirection().y());
                    analysis->fill_histogram_mZ(i_part, thePrePoint->GetMomentumDirection().z());

                    analysis->counter_total[i_part]++;

                    settings->event_lead_to_detection[i_part] = true;

                    if (momy > 0.0) {
                        analysis->counter_up[i_part]++;
                    } else {
                        analysis->counter_down[i_part]++;
                    }
                }
            }
        }

    }
}

// ------------------------------------------------------------------------

bool SteppingAction::is_inside_eField_region(const double &alt, const double &xx, const double &zz)
// alt, xx, zz assumed in km
{
    return alt > EFIELD_alt_min && alt < EFIELD_alt_max && std::abs(xx) < settings->EFIELD_XZ_HALF_SIZE && std::abs(zz) < settings->EFIELD_XZ_HALF_SIZE;
}

// ------------------------------------------------------------------------

index_found SteppingAction::find_particle_index(const int PDG_in) {

    index_found my_index_found{1000, false};

    for (uint ii = 0; ii < settings->PDG_LIST.size(); ++ii) {
        if (PDG_in == settings->PDG_LIST[ii]) {
            my_index_found.index = ii;
            my_index_found.found = true;
            return my_index_found;
        }
    }

    return my_index_found;
}

// ------------------------------------------------------------------------
G4String SteppingAction::get_name(int PDG) {
    G4String name;
    if (PDG == 22) name = "photon";
    else if (PDG == 11) name = "electron";
    else if (PDG == -11) name = "positron";
    else if (PDG == -13) name = "MuonP";
    else if (PDG == 13) name = "MuonN";
    else if (PDG == 2212) name = "proton";
    else if (PDG == 2112) name = "neutron";
    else std::abort();

    return name;
}