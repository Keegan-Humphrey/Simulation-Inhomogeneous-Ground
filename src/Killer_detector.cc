#include "Killer_detector.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"



Killer_detector::Killer_detector(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="Killer_detectorC");
  //HCID = -1;
}


Killer_detector::~Killer_detector(){;}

void Killer_detector::Initialize(G4HCofThisEvent*HCE)
{
  /*hitsCollection = new Killer_detectorHitsCollection
                   (SensitiveDetectorName,collectionName[0]);
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection); }
  HCE->AddHitsCollection(HCID,hitsCollection);*/
}

G4bool Killer_detector::ProcessHits(G4Step* aStep,G4TouchableHistory* /*ROhist*/)
{
  G4int copy_num = aStep->GetTrack()->GetVolume()->GetCopyNo();
  
  
  if ( copy_num==999){ 
  	aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
  	G4cout << "Killing particle within Killer Detector" <<G4endl;
  }  
  return true;

}

void Killer_detector::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{;}
