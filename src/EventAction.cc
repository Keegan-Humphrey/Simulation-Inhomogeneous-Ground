#include "EventAction.hh"

#include "EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ThreeVector.hh"
#include "G4VVisManager.hh"

//#include <CLHEP/Units/PhysicalConstants.h>
#include "G4PhysicalConstants.hh"

using CLHEP::twopi;

extern std::vector<G4int> DAQs; /* Temporary Array to Store DAQs. Size: max(Strips) x 3  */
extern G4ThreeVector upper_hitting_point;
extern G4ThreeVector lower_hitting_point;
//extern std::vector<G4double> BarsEventStatus;
extern G4double nBars;
extern G4double UD_Det_planes_seperation;
extern std::vector <std::vector<G4double> > BarsEventStatus;
/*projection related*/
extern std::vector<std::vector<G4int> > DetectorCounts;
extern G4String RD2File;
extern std::vector<G4double> Iaxis;
extern std::vector<G4double> Jaxis;
extern G4double UD_Det_beta;
extern G4int UD_ProjectionPixelI;
extern G4int UD_ProjectionPixelJ;
extern G4int printCounter;
extern G4long LegalEvents;

extern unsigned long cnt_total_events;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
	:printModulo(10000), eventMessenger(0)
{
	eventMessenger = new EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
	delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
	G4int evtNb = evt->GetEventID();

	//printing survey
	if (evtNb%printModulo == 0)
		G4cout << "\n---> Begin of Event: " << evtNb << G4endl;

	// Nulling the hitting points
	upper_hitting_point = G4ThreeVector(0.0, 0.0, 0.0);
	lower_hitting_point = G4ThreeVector(0.0, 0.0, 0.0);
	DAQs[0] = 0;
	DAQs[1] = 0;
	for (int j = 0; j < 4; j++)
	{
		for (int i = 0; i < nBars; i++)
		{
			BarsEventStatus[j*nBars + i][0] = 100 * (j + 1) + i;
			BarsEventStatus[j*nBars + i][1] = 0;
			BarsEventStatus[j*nBars + i][2] = 0;
		}
	}
	return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
	/* VISUALIZATION*/
	if (G4VVisManager::GetConcreteInstance())
	{
		G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
		G4int n_trajectories = 0;
		if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
		for (G4int i = 0; i < n_trajectories; i++) {
			G4Trajectory* trj = (G4Trajectory*)
				((*(evt->GetTrajectoryContainer()))[i]);
			trj->DrawTrajectory();
		}
	}

	/******** PROCESSSING HITS FROM THIS EVENT  *********/

	/* BENCHMARK - THE REAL TRAJECTORY OF THIS EVENT*/
	G4PrimaryVertex *v = evt->GetPrimaryVertex(0);
	G4PrimaryParticle *par = v->GetPrimary(0);
	G4ThreeVector pos = v->GetPosition();
	G4ThreeVector momentum = par->GetMomentum();

	/*G4String ofileName = G4String("HitPositions");
	std::ofstream outFile1(ofileName, std::ios::app);  */

	/******** TRAJECTORY RECONSTRUCTION FROM DAQs ************/
	/*** Currently, Using Triplets only */
	if ((DAQs[0] == 1) && (DAQs[1] == 1)) {
		LegalEvents++;
		G4ThreeVector direction = upper_hitting_point - lower_hitting_point;
		G4ThreeVector center = 0.5*(upper_hitting_point + lower_hitting_point);

		G4int EvtId = evt->GetEventID();
		//G4cout << "Evt: " <<  EvtId << ". DIR " << direction.x() << "  "  << direction.y() << "  "  << direction.z() << "  " << G4endl;

		//Calculating local dx, dy on planes
		G4RotationMatrix Detector_container_rotm = G4RotationMatrix();
		Detector_container_rotm.rotateX(0.0);
		Detector_container_rotm.rotateY(UD_Det_beta);
		Detector_container_rotm.rotateZ(0.0);

		G4ThreeVector up_axis = Detector_container_rotm * G4ThreeVector(0.0, 0.0, UD_Det_planes_seperation/2);
		G4ThreeVector dw_axis = Detector_container_rotm * G4ThreeVector(0.0, 0.0, -UD_Det_planes_seperation/2);
		G4ThreeVector y_axis = Detector_container_rotm * G4ThreeVector(0.0, 1.0, 0.0);
		G4ThreeVector x_axis = Detector_container_rotm * G4ThreeVector(1.0, 0.0, 0.0);

		G4ThreeVector diff_up = upper_hitting_point - up_axis;
		G4double x_proj_up = diff_up * x_axis;
		G4double y_proj_up = diff_up * y_axis;

		G4ThreeVector diff_dw = lower_hitting_point - dw_axis;
		G4double x_proj_dw = diff_dw * x_axis;
		G4double y_proj_dw = diff_dw * y_axis;

		G4double dx = x_proj_up - x_proj_dw;
		G4double dy = y_proj_up - y_proj_dw;

		/*  	G4cout << "Evt: " <<  EvtId << "   " << G4endl;
			G4cout << "UPPER POINT " << upper_hitting_point <<  G4endl;
			G4cout << "up_axis  " << up_axis <<  G4endl;
			G4cout << "diff_up  " << diff_up <<  G4endl;
			G4cout << "x_axis  " << x_axis <<  G4endl;
			G4cout << "y_axis  " << y_axis <<  G4endl;
			G4cout << "x_proj_up  " << x_proj_up <<  G4endl;
			G4cout << "y_proj_up  " << y_proj_up <<  G4endl;

			G4cout << "LOWER POINT " << lower_hitting_point <<  G4endl;
			G4cout << "dw_axis  " << dw_axis <<  G4endl;
			G4cout << "diff_dw  " << diff_dw <<  G4endl;
			G4cout << "x_axis  " << x_axis <<  G4endl;
			G4cout << "y_axis  " << y_axis <<  G4endl;
			G4cout << "x_proj_dw  " << x_proj_dw <<  G4endl;
			G4cout << "y_proj_dw  " << y_proj_dw <<  G4endl;*/

			//std::ofstream outFile1("temp.dat", std::ios::app);
			//outFile1  << (G4int)round(x_proj_up) << "  "  << (G4int)round(y_proj_up) << "  "  << (G4int)round(x_proj_dw) << "    "  << (G4int)round(y_proj_dw) << G4endl; 



			/* Add the counts to the REAL ARRAY*/
		G4int IindUp = (G4int)round((x_proj_up - Iaxis[0]) / (Iaxis[1] - Iaxis[0])); // Binning according to signel plate!
		G4int JindUp = (G4int)round((y_proj_up - Jaxis[0]) / (Jaxis[1] - Jaxis[0]));

		G4int IindDw = (G4int)round((x_proj_dw - Iaxis[0]) / (Iaxis[1] - Iaxis[0]));
		G4int JindDw = (G4int)round((y_proj_dw - Jaxis[0]) / (Jaxis[1] - Jaxis[0]));

		//G4cout  << IindUp << "   " <<  JindUp  << "  "  << IindDw << "   " <<  JindDw << "...."; 

		//G4cout  << Iaxis[0] << "   " <<  dx   << "  "  << (Iaxis[1]-Iaxis[0]) << "   " <<  Iind << "    " << Jaxis[0] << "    "  << dy << (Jaxis[1]-Jaxis[0]) << "   " <<  Jind << G4endl; 
		G4int Iind = IindUp - IindDw + UD_ProjectionPixelI;
		G4int Jind = JindUp - JindDw + UD_ProjectionPixelJ;
		//G4cout << "   " << Iind << "    "   << Jind   <<  "   " <<  G4endl;
		if ((Iind >= 0) && (Jind >= 0)) {
			if ((Iind < DetectorCounts.size()) && (Jind < DetectorCounts[0].size())) {
				DetectorCounts[Iind][Jind]++;
				//G4cout  << "Cts[" <<  Iind << "][" << Jind << "]=" <<  DetectorCounts[Iind][Jind] << G4endl; 
			}
		}


		G4String tmpLine = "";
		time_t _tm = time(NULL);
		struct tm * curtime = localtime(&_tm);
		tmpLine = "*Event*\n" + std::to_string(LegalEvents) + "\n" + (G4String)asctime(curtime);
		
		for (int i = 0; i < 4*nBars; i++)
		{
			if (BarsEventStatus[i][1] == 1)
			{
				tmpLine += std::to_string((int)BarsEventStatus[i][0]) + "," + std::to_string(BarsEventStatus[i][2]) + "\n";
			}
		}
		tmpLine += std::to_string(upper_hitting_point[0]) + "," + std::to_string(upper_hitting_point[1]) + "," + std::to_string(upper_hitting_point[2]) + "\n";
		tmpLine += std::to_string(lower_hitting_point[0]) + "," + std::to_string(lower_hitting_point[1]) + "," + std::to_string(lower_hitting_point[2]) + "\n";
		RD2File += tmpLine;
		




		if ((evt->GetEventID()) > printCounter*printModulo) {
			//G4cout << "Print To Files!!! "  << G4endl;
			printCounter++;

			std::ofstream outCnt2all("Cts_dx_dy.out");
			outCnt2all << cnt_total_events << G4endl;
			for (G4int j = 0; j < DetectorCounts.size(); j++) {
				for (G4int i = 0; i < DetectorCounts[j].size(); i++) {
					outCnt2all << DetectorCounts[j][i] << "   ";
				}
				outCnt2all << G4endl;
			}
			std::ofstream RowData("RowData.out", std::ios_base::app);
			RowData << RD2File << G4endl;
			RD2File = "";
		}
	}
	return;

}
