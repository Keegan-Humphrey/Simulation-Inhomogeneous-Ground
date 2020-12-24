#include "StackingAction.hh"
#include "G4Track.hh"
#include "G4String.hh"
#include "G4TrackStatus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  G4ClassificationOfNewTrack     classification = fUrgent;

  // kill all secondaries
  // i added it, try to kill all particle which are not mu-
  G4String ParticleName = aTrack->GetDefinition()->GetParticleName();
  if (ParticleName.compare("mu-"))
  { 
	  classification = fKill;
  }
  //if(aTrack->GetParentID() != 0)  classification = fKill;
/*
  {
	  G4int copy_num = aTrack->GetVolume()->GetCopyNo();
	  if (copy_num<1)
	  	classification = fKill;
  }  */
  return classification;
}
