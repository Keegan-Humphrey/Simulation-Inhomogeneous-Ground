#ifndef Killer_detector_h
#define Killer_detector_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;
//class G4HCofThisEvent;
class G4TouchableHistory;

class Killer_detector : public G4VSensitiveDetector
{

  public:
      Killer_detector(G4String name);
      virtual ~Killer_detector();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      //Killer_detectorHitsCollection* hitsCollection;
      //G4int HCID;
};




#endif
