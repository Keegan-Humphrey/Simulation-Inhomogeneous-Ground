#include "Strip_detector.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"

//#include "G4SystemOfUnits.hh"


extern std::vector<G4int> DAQs;
extern G4ThreeVector upper_hitting_point;
extern G4ThreeVector lower_hitting_point;
extern G4double nBars;
//extern std::vector<G4double> BarsEventStatus;
extern std::vector <std::vector<G4double> > BarsEventStatus;



Strip_detector::Strip_detector(G4String name)
	:G4VSensitiveDetector(name)
{
	G4String HCname;
	collectionName.insert(HCname = "Strip_detectorC");
	//HCID = -1;
}


Strip_detector::~Strip_detector() { ; }

void Strip_detector::Initialize(G4HCofThisEvent*HCE)
{
	return;
}


G4bool Strip_detector::ProcessHits(G4Step* aStep, G4TouchableHistory* /*ROhist*/)
{
	//G4String ofileName = G4String("Out_")+Strip_detector::GetName();
	//std::ofstream outFile1(ofileName, std::ios::app);

	//static G4int last_evt_ID=-1;

	//G4double R=120*CLHEP::cm;
	//G4double center=0.0*CLHEP::cm;
	G4StepPoint* preStepPoint = aStep->GetPreStepPoint();  //get the step point
	G4StepPoint* postStepPoint = aStep->GetPostStepPoint();	//for my debug, delete it later
	G4ThreeVector position = preStepPoint->GetPosition();
	G4double Time1 = preStepPoint->GetGlobalTime();
	G4int TrkID = aStep->GetTrack()->GetTrackID();
	G4int EvtId = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
	G4int copy_num = aStep->GetTrack()->GetVolume()->GetCopyNo();
	// kill the track if the particle is not a mu-
	G4String ParticleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
	if (ParticleName.compare("mu-"))
	{
		//aStep->GetTrack()->SetTrackStatus(fStopAndKill);
		aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
		return true;
	}
	/*G4String name = aStep->GetTrack()->GetDefinition()->GetParticleName();*/
	//if (  aStep->IsFirstStepInVolume()  ||  aStep->IsLastStepInVolume()  ){ /*THIS CURRENTLY DISABLES THE POSSIBILITY OF CROSSING THE STRIP TWICE.*/
	if (copy_num == 10)
	{
		DAQs[0] = 1;
		upper_hitting_point = preStepPoint->GetPosition(); //from the step point get the (x,y,z) position
	}
	if (copy_num == 20)
	{
		DAQs[1] = 1;
		lower_hitting_point = preStepPoint->GetPosition(); //from the step point get the (x,y,z) position
	}
	if (copy_num >= 100 && copy_num < 500)
	{
		int BarInd = copy_num % 100;
		int LayerInd = floor(copy_num / 100);

		BarsEventStatus[(LayerInd - 1)*nBars + BarInd][1] = 1;
		//BarsEventStatus[(LayerInd - 1)*nBars + BarInd][2] += aStep->GetStepLength();
		BarsEventStatus[(LayerInd - 1)*nBars + BarInd][2] += aStep->GetTotalEnergyDeposit();

	}
	return true;
}

void Strip_detector::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
	;
}
