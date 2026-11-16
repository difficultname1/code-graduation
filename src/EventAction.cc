#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"

void EventAction::BeginOfEventAction(const G4Event* event)
{
  fEventID = event->GetEventID();
}

void EventAction::AddHit(G4double x, G4double y, G4double z, G4double edep,
                         G4double t)
{
  fRun->FillHit(fEventID, x, y, z, edep, t);
}
