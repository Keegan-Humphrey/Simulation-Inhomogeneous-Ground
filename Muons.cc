#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RE02RunAction.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "EventAction.hh"
#include "StackingAction.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "G4PhysicalConstants.hh"

using CLHEP::twopi;


unsigned long cnt_total_events;
//G4double StoreE[10000000];
//#define Nstore1 10000001
//G4double StoreT[Nstore1];








G4double UD_SourceLengthScale = 0.5*m;
G4int UD_ProjectionPixelI = 144;
G4int UD_ProjectionPixelJ = 144;


//----- USER (MINE) VARIABLES - STATIC, SHARED BY ALL THE COMPONENTS OG GEANT4 -----
std::vector<G4int> DAQs; /* Temporary Array to Store DAQs. Size: max(Strips) x 3  */
G4ThreeVector upper_hitting_point;
G4ThreeVector lower_hitting_point;
std::vector<std::vector<G4int> > DetectorCounts;
G4String RD2File = "";
G4long LegalEvents = 0;
std::vector<G4double> Iaxis;
std::vector<G4double> Jaxis;


G4int printCounter = 0;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv) {


	cnt_total_events = 0;

	
	DAQs.resize(2);


	/*Initialize counts */

	Iaxis.resize(UD_ProjectionPixelI);
	for (G4int Iind = 0; Iind < Iaxis.size(); Iind++) {
		Iaxis[Iind] = -UD_SourceLengthScale / 2.0 + UD_SourceLengthScale / (Iaxis.size() - 1)*Iind;
	}
	Jaxis.resize(UD_ProjectionPixelJ);
	for (G4int Jind = 0; Jind < Jaxis.size(); Jind++) {
		Jaxis[Jind] = -UD_SourceLengthScale / 2.0 + UD_SourceLengthScale / (Jaxis.size() - 1)*Jind;
	}


	/* Initialize arrays*/
	DetectorCounts.resize(2 * UD_ProjectionPixelI - 1);
	for (G4int Iind = 0; Iind < DetectorCounts.size(); Iind++) {
		DetectorCounts[Iind].resize(2 * UD_ProjectionPixelJ - 1);
		for (G4int Jind = 0; Jind < DetectorCounts[Iind].size(); Jind++) {
			DetectorCounts[Iind][Jind] = 0;
		}
	}
	// Initialize row data file
	std::ofstream RowDataInit("RowData.out");
	RowDataInit << RD2File << G4endl;
	//my Verbose output class
	//G4VSteppingVerbose::SetInstance(new SteppingVerbose);

	// Construct the default run manager
	G4RunManager * runManager = new G4RunManager;

	// set mandatory initialization classes
	DetectorConstruction* det;
	PrimaryGeneratorAction* prim;
	runManager->SetUserInitialization(det = new DetectorConstruction);
	runManager->SetUserInitialization(new PhysicsList);
	runManager->SetUserAction(prim = new PrimaryGeneratorAction(det));

#ifdef G4VIS_USE
	// visualization manager
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();
#endif

	//HistoManager* histo = new HistoManager();

	// set user action classes
	//RunAction* run;
	//runManager->SetUserAction(run = new RunAction(det,prim,histo));
	runManager->SetUserAction(new RE02RunAction);
	runManager->SetUserAction(new EventAction);
	//runManager->SetUserAction(new SteppingAction(histo));
	runManager->SetUserAction(new StackingAction);

	// Start execution
	//
	if (argc > 1) {	// execute an argument macro file if exist
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		G4UImanager::GetUIpointer()->ApplyCommand(command + fileName);

	}
	else
	{
		G4String command = "/control/execute ";
		G4String fileName = "r2.in";
		G4UImanager::GetUIpointer()->ApplyCommand(command + fileName);
	}
	//	else {		// start interactive session
	//		G4UIsession* session = 0;
	//#ifdef G4UI_USE_TCSH
	//		session = new G4UIterminal(new G4UItcsh);
	//#else
	//		session = new G4UIterminal();
	//#endif
	//		session->SessionStart();
	//		delete session;
	//	}

		// job termination
		//

#ifdef G4VIS_USE
	delete visManager;
#endif

	//delete histo;
	delete runManager;
	G4cout << "total number of events= " << cnt_total_events << G4endl;

	//G4double pi=3.1416;

	std::ofstream outCnt2all("Cts_dx_dy.out");
	outCnt2all << cnt_total_events << G4endl;
	for (G4int j = 0; j < (G4int)DetectorCounts.size(); j++) {
		for (G4int i = 0; i < (G4int)DetectorCounts[j].size(); i++) {
			outCnt2all << DetectorCounts[j][i] << "   ";
		}
		outCnt2all << G4endl;
	}
	std::ofstream RowData("RowData.out", std::ios_base::app);
	RowData << RD2File << G4endl;
	RD2File = "";


	/*std::ofstream SST("StoreT.out");
	SST.precision(16);
	for (G4int j=0;j<Nstore1;j++){
	  SST << StoreT[j] <<G4endl;
	}*/


	return 0;
}
