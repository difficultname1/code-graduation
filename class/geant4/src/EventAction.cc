#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"

void EventAction::BeginOfEventAction(const G4Event* event)
{
  fEventID = event->GetEventID();
}

void EventAction::AddHit(const G4ThreeVector& p0, const G4ThreeVector& p1,
                         G4double fracInside, G4double edep,
                         G4double t0, G4double t1)
{
  fRun->FillHit(fEventID, p0, p1, fracInside, edep, t0, t1);
}
