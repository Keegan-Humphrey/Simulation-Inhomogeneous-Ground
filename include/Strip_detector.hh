#ifndef Strip_detector_h
#define Strip_detector_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;
//class G4HCofThisEvent;
class G4TouchableHistory;

class Strip_detector : public G4VSensitiveDetector
{

  public:
      Strip_detector(G4String name);
      virtual ~Strip_detector();

      virtual void Initialize(G4HCofThisEvent*HCE);
      virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      virtual void EndOfEvent(G4HCofThisEvent*HCE);

  private:
      //Strip_detectorHitsCollection* hitsCollection;
      //G4int HCID;
};




#endif
