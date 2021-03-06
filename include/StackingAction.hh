#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"
#include "globals.hh"

class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StackingAction : public G4UserStackingAction
{
  public:

    StackingAction();
   ~StackingAction();

    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* );
    
};

#endif
