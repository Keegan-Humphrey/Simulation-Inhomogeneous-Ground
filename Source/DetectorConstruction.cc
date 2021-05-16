#include "G4Tet.hh"
 
#include "Relief.hh"

#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Ellipsoid.hh"
#include "G4Polycone.hh"
#include "G4Sphere.hh"
#include "G4TwistedTubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "Strip_detector.hh"
#include "Killer_detector.hh"

#include "G4PhysicalConstants.hh"

using CLHEP::twopi;

extern std::vector<std::vector<G4int> > DAQs; /* Temporary Array to Store DAQs. Size: max(Strips) x 3  */


G4double UD_World_sizeX = 300 * m;
G4double UD_World_sizeY = 300 * m;
G4double UD_World_sizeZ = 300 * m;
G4double UD_Det_planes_seperation = 25*cm;
G4double BarWidth = 32 * mm;
G4double BarHight = 17 * mm;
G4double BarDepth = 40 * cm;
G4double nBars = 27;
//G4double Trigger_size = (std::max(BarDepth, (nBars - 1)*BarWidth / 2)) + 5 * cm;
G4double Trigger_size = (std::max(BarDepth, (nBars - 1)*BarWidth / 2)) ;
G4double Trigger_width = 1 * cm;
G4double DetectorDepth = 5.0*m;
G4double UD_Det_centerX = 0*cm;
G4double UD_Det_centerY = 0*cm;
G4double UD_Det_centerZ = Trigger_width + 2 * BarHight + UD_Det_planes_seperation / 2;
G4double UD_Det_alpha = 0*deg;
G4double UD_Det_beta = -22.5*deg;
G4double UD_Det_gamma = 0*deg;

extern std::vector<G4double> Iaxis;
extern std::vector<G4double> Jaxis;

//std::vector<G4double> BarsEventStatus;
std::vector <std::vector<G4double> > BarsEventStatus(4*nBars, std::vector<G4double>(3));


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
	/*Initialize counts with the right trigger size... */
	for (G4int Iind = 0; Iind < Iaxis.size(); Iind++) {
		Iaxis[Iind] = -Trigger_size / 2.0 + Trigger_size / (Iaxis.size() - 1)*Iind;
	}
	for (G4int Jind = 0; Jind < Jaxis.size(); Jind++) {
		Jaxis[Jind] = -Trigger_size / 2.0 + Trigger_size / (Jaxis.size() - 1)*Jind;
	}
	return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	//--------- Material definition ---------

	G4double a, z;
	G4double density, abundance, fractionmass;
	G4int nel, ncomponents, i;
	G4int natoms;


	// IMAGING THE USER PARAMETERS (STATIC - EXTERNAL, SHARED BY ALL PARTS OF GEANT4) 
	// WHICH COME FROM THE CARDIN INTO THE INNER PARAMETERS OF DETECTOR
	//G4cout << "DetectorConstruction::World_sizeX = " << UD_World_sizeX << G4endl;
	//G4cout << "DetectorConstruction::World_sizeY = " << UD_World_sizeY << G4endl;
	//G4cout << "DetectorConstruction::World_sizeZ = " << UD_World_sizeZ << G4endl;
	//G4cout << "DetectorConstruction::UD_Det_alpha = " << UD_Det_alpha << G4endl;
	//G4cout << "DetectorConstruction::UD_Det_beta = " << UD_Det_beta << G4endl;
	//G4cout << "DetectorConstruction::UD_Det_gamma = " << UD_Det_gamma << G4endl;
	//G4cout << "DetectorConstruction::UD_Det_planes_width = " << UD_Det_planes_width << G4endl;
	//G4cout << "DetectorConstruction::UD_Det_planes_size = " << UD_Det_planes_size << G4endl;
	//G4cout << "DetectorConstruction::UD_Det_planes_seperation = " << UD_Det_planes_seperation << G4endl;


	G4NistManager* nist = G4NistManager::Instance();

	//Elements
	/////////////////////////////
	G4Element* H = new G4Element("Hydrogen", "H", z = 1., a = 1.00794*g / mole);
	G4Element* C = new G4Element("Carbon", "C", z = 6., a = 12.011*g / mole);
	G4Element* N = new G4Element("Nitrogen", "N", z = 7., a = 14.0067*g / mole);
	G4Element* O = new G4Element("Oxygen", "O", z = 8., a = 16.00*g / mole);
	G4Element* Ca = new G4Element("Calcium", "Ca", z = 20., a = 40.08*g / mole);

	G4Material* U = nist->FindOrBuildMaterial("G4_U");


	//Materials
	////////////////////////////////



	//Delrin 1:0.06713   6:0.40002   8:0.53285
	G4Material* Delrin = new G4Material("Delrin", density = 1.425*g / cm3, ncomponents = 3);
	Delrin->AddElement(H, fractionmass = 0.06713);
	Delrin->AddElement(C, fractionmass = 0.40002);
	Delrin->AddElement(O, fractionmass = 0.53285);


	//Air
	G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");

	//Vacuum
	G4Material* Vacuum =
		new G4Material("Galactic", z = 1., a = 1.01*g / mole, density = universe_mean_density,
			kStateGas, 2.73*kelvin, 3.e-18*pascal);

	//Solid Rock
	G4double rho = 1.5;
	G4Material* GroundMat = new G4Material("SolidRock", rho*g / cm3, ncomponents = 3);
	GroundMat->AddElement(Ca, natoms = 1);
	GroundMat->AddElement(C, natoms = 1);
	GroundMat->AddElement(O, natoms = 3);

	//Bad Mother Fucker 
	G4Material* BMF = new G4Material("BMF", 10000*g / cm3, ncomponents = 1);
	BMF->AddElement(Ca, natoms = 1);
	
	 
	
	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;



	////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //--------- Materials of the principal geometrical components (solids)  ---------
	DetMatter = Delrin;  //amplifier Matter
	////////////////////////////////////////////////////////////////////////////////////////////////////////


	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	G4String SDname;


	//------------------------------
	// World
	//------------------------------
	solidWorld = new G4Box("world", UD_World_sizeX, UD_World_sizeY, UD_World_sizeZ);
	logicWorld = new G4LogicalVolume(solidWorld, Vacuum, "World");

	//  Must place the World Physical volume unrotated at (0,0,0).
	//
	physiWorld = new G4PVPlacement(0,               // no rotation
		G4ThreeVector(), // at (0,0,0)
		logicWorld,      // its logical volume
		"World",         // its name
		0,               // its mother  volume
		false,           // no boolean operations
		0);              // copy number


////    GROUND CONTAINER - Vacuum
//	G4double DetectorDepth = 1.0*m;
//	G4VSolid* solidMotherGround = new G4Box("Ground", UD_World_sizeX, UD_World_sizeY, DetectorDepth / 2.0);
//	G4LogicalVolume* logicMotherGround = new G4LogicalVolume(solidMotherGround, Vacuum, "Ground", 0, 0, 0);
//	G4ThreeVector GroundCenter = G4ThreeVector(0.0, 0.0, DetectorDepth / 2.0);
//	new G4PVPlacement(0,               // no rotation
//		GroundCenter, // at (0,0,0)
//		logicMotherGround,      // its logical volume
//		"GroundMother",         // its name
//		logicWorld,               // its mother  volume
//		false,           // no boolean operations
//		0);              // copy number

	/* HOMOGENEUOUS GROUND */
	G4VSolid* solidGround = new G4Box("HomoGround", UD_World_sizeX, UD_World_sizeY, DetectorDepth / 2);
	G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, GroundMat, "HomoGround");
	//G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, Vacuum, "HomoGround");
	new G4PVPlacement(0,               // no rotation
		G4ThreeVector(0.0, 0.0, DetectorDepth / 2.0), //origin
		logicGround,
		"HomoGround",         // its name
		logicWorld,  // its mother  volume
		false,           // no boolean operations
		0);              // copy number   		                                 

	/* U Box that i added */
	G4double UBoxHight = 2 * m;
	G4VSolid* UBox = new G4Box("UBox", 0.5*m, 50*m, UBoxHight/2);
	//G4VSolid* UBox = new G4Box("UBox", 1*m, 1*m, UBoxHight);
	G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, Air, "UBox");
	//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, GroundMat, "UBox");
	//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, BMF, "UBox");
	//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, Vacuum, "UBox");
	new G4PVPlacement(0,               // no rotation
                      G4ThreeVector(2 * m, 0.0, DetectorDepth / 2 - UBoxHight / 2), //origin
		      //G4ThreeVector(0.0, 2 * m, 3 * UBoxHight - UBoxHight / 2), //origin
		UBoxLogicVolume,
		"UBox",         // its name
		logicGround,  // its mother  volume
		false,           // no boolean operations
		0);

	//    DETECTOR CONTAINER - Air
	//G4double Detector_container_height = 2 * UD_Det_planes_width + UD_Det_planes_seperation;
	G4double Detector_container_height = 4 * BarHight + UD_Det_planes_seperation + 2 * Trigger_width;
	G4double Detector_container_width = Trigger_size;
	G4VSolid* solid_Detector_container = new G4Box("Det_container", Detector_container_width / 2.0, Detector_container_width / 2.0, Detector_container_height / 2.0);
	G4LogicalVolume* logic_Detector_container = new G4LogicalVolume(solid_Detector_container, Vacuum, "Det_container");
	G4ThreeVector Detector_container_center = G4ThreeVector(UD_Det_centerX, UD_Det_centerY, -DetectorDepth / 2.0 + Detector_container_height / 2);
	
    /*Rotate the Detector*/
    G4RotationMatrix* transform = new G4RotationMatrix();
    transform->rotateX(UD_Det_alpha);
    transform->rotateY(UD_Det_beta);
    transform->rotateZ(UD_Det_gamma);
    
    //else 0, for transform
    new G4PVPlacement(transform,
        Detector_container_center,
        logic_Detector_container,      // its logical volume
        "Det_container",         // its name
        logicGround,// its mother  volume
        false,           // no boolean operations
        0);              // copy number

	G4VSensitiveDetector* WC_Det = new Strip_detector(SDname = "Wires_Chamber");
	SDman->AddNewDetector(WC_Det);


	//    UPPER PLANE
		//Upper trigger
	G4VSolid* UpperTrig = new G4Box("UpperTrig", Trigger_size / 2, Trigger_size / 2, Trigger_width / 2);
	G4LogicalVolume* UpperTrigLogicVolume = new G4LogicalVolume(UpperTrig, DetMatter, "UpperTrig");
	new G4PVPlacement(0,               // no rotation
		G4ThreeVector(0.0, 0.0, Detector_container_height / 2 - Trigger_width / 2), //origin
		UpperTrigLogicVolume,
		"UpperTrig",         // its name
		logic_Detector_container,  // its mother  volume
		false,           // no boolean operations
		10);              // copy number
	UpperTrigLogicVolume->SetSensitiveDetector(WC_Det);


		//Bars
		// Define the contours of the bar
	std::vector<G4TwoVector> newBar(3);
	newBar[0] = G4TwoVector(-BarWidth/2, -BarHight/2);
	newBar[1] = G4TwoVector(BarWidth/2, -BarHight / 2);
	newBar[2] = G4TwoVector(0.0, BarHight/2);
	G4ExtrudedSolid* solid_upper_plane = new G4ExtrudedSolid("upper_plane", newBar, BarDepth/2, G4TwoVector(0, 0), 1.0, G4TwoVector(0, 0), 1.0);
	G4LogicalVolume* logic_upper_plane = new G4LogicalVolume(solid_upper_plane, DetMatter, "upper_plane");

	for (int i = 0; i < nBars; i++)
	{
		G4RotationMatrix* rotBarX = new G4RotationMatrix();
		rotBarX->rotateX(pow(-1, i)*90.*deg);
		G4ThreeVector upper_plane_X_center = G4ThreeVector(-(nBars / 4 - 0.25) * BarWidth + BarWidth / 2 * i, 0, UD_Det_planes_seperation / 2.0 + BarHight * 3 / 2.0);
		new G4PVPlacement(rotBarX,
			upper_plane_X_center,
			logic_upper_plane,      // its logical volume
			"upper_plane",         // its name
			logic_Detector_container,               // its mother  volume
			false,           // no boolean operations
			i+100);              // copy number
		G4RotationMatrix* rotBarY = new G4RotationMatrix();
		rotBarY->rotateX(pow(-1, i)*90.*deg);
		rotBarY->rotateY(90.*deg);
		G4ThreeVector upper_plane_Y_center = G4ThreeVector(0, -(nBars / 4 - 0.25) * BarWidth + BarWidth / 2 * i, UD_Det_planes_seperation / 2.0 + BarHight / 2.0);
		new G4PVPlacement(rotBarY,
			upper_plane_Y_center,
			logic_upper_plane,      // its logical volume
			"upper_plane",         // its name
			logic_Detector_container,               // its mother  volume
			false,           // no boolean operations
			i+200);
		/*logic_upper_plane->SetSensitiveDetector(WC_Det);*/
	}
	logic_upper_plane->SetSensitiveDetector(WC_Det);

	//    LOWER PLANE
		//Lower trigger
	G4VSolid* LowerTrig = new G4Box("LowerTrig", Trigger_size / 2, Trigger_size / 2, Trigger_width / 2);
	G4LogicalVolume* LowerTrigLogicVolume = new G4LogicalVolume(LowerTrig, DetMatter, "LowerTrig");
	new G4PVPlacement(0,               // no rotation
		G4ThreeVector(0.0, 0.0, -Detector_container_height / 2 + Trigger_width / 2), //origin
		LowerTrigLogicVolume,
		"LowerTrig",         // its name
		logic_Detector_container,  // its mother  volume
		false,           // no boolean operations
		20);              // copy number
	LowerTrigLogicVolume->SetSensitiveDetector(WC_Det);

	G4ExtrudedSolid* solid_lower_plane = new G4ExtrudedSolid("lower_plane", newBar, BarDepth / 2, G4TwoVector(0, 0), 1.0, G4TwoVector(0, 0), 1.0);
	G4LogicalVolume* logic_lower_plane = new G4LogicalVolume(solid_lower_plane, DetMatter, "lower_plane");


	for (int i = 0; i < nBars; i++)
	{
		G4RotationMatrix* rotBarX = new G4RotationMatrix();
		rotBarX->rotateX(pow(-1, i)*90.*deg);
		G4ThreeVector lower_plane_X_center = G4ThreeVector(-(nBars / 4 - 0.25) * BarWidth + BarWidth / 2 * i, 0, -UD_Det_planes_seperation / 2.0 - BarHight / 2.0);
		new G4PVPlacement(rotBarX,
			lower_plane_X_center,
			logic_lower_plane,      // its logical volume
			"lower_plane",         // its name
			logic_Detector_container,               // its mother  volume
			false,           // no boolean operations
			i+300);              // copy number
		G4RotationMatrix* rotBarY = new G4RotationMatrix();
		rotBarY->rotateX(pow(-1, i)*90.*deg);
		rotBarY->rotateY(90.*deg);
		G4ThreeVector lower_plane_Y_center = G4ThreeVector(0, -(nBars / 4 - 0.25) * BarWidth + BarWidth / 2 * i, -UD_Det_planes_seperation / 2.0 - BarHight * 3 / 2.0);
		new G4PVPlacement(rotBarY,
			lower_plane_Y_center,
			logic_lower_plane,      // its logical volume
			"lower_plane",         // its name
			logic_Detector_container,               // its mother  volume
			false,           // no boolean operations
			i + 400);
		/*logic_lower_plane->SetSensitiveDetector(WC_Det);*/
	}
	logic_lower_plane->SetSensitiveDetector(WC_Det);
	

	//using namespace RData;
    for (int i = 0; i < iaxis - 1; i++)
    {
        //using namespace RData;

        for (int j = 0; j < jaxis - 1; j++)
        {
            
            // Makes boxes between four relief data points (above homogeneous ground),
            // takes height to be the lowest of the four (RMins[i][j])
            if (RMins[i][j] > 0.0*m)
            {
                G4VSolid* RBox = new G4Box("RBox", Pnt_Sep/2, Pnt_Sep/2, RMins[i][j]/2);
                G4LogicalVolume* RBoxLogicVolume = new G4LogicalVolume(RBox, GroundMat, "RBox");
                new G4PVPlacement(0,               // no rotation
                        G4ThreeVector((1-iaxis) * Pnt_Sep /2 + Pnt_Sep * i,(1-jaxis) * Pnt_Sep /2 + Pnt_Sep * j, DetectorDepth + RMins[i][j]/2),
                        RBoxLogicVolume,
                        "RBox",         // its name
                        logicWorld,  // its mother  volume
                        false,           // no boolean operations
                        0);
            }
            
            //////////////
            // Makes surface continuous by placing tetraheadra on top of each box
            auto make_tets = [&](int anchorx, int anchory, double signx, double signy)
            {
                // anchorx, anchory - indices of the anchor point (minima of the adjacent data points)
                // signx, signy - (anchorx, anchory > i,j) ==> +1.0, else -1.0
                //              - relative direction of the other data points from anchor
                
                
                // make a tetrahedron and place it at the anchor point
                auto make_a_tet = [&](G4ThreeVector p1, G4ThreeVector p2, G4ThreeVector p3)
                {
                    // p[i] - points of the tetrahedron relative to anchor
                    G4VSolid* RTet = new G4Tet("RTet", G4ThreeVector {0, 0, 0}, p1, p2, p3, 0);
                    G4LogicalVolume* RTetLogicVolume = new G4LogicalVolume(RTet, GroundMat, "RTet");
                    new G4PVPlacement(0,               // no rotation
                            G4ThreeVector((2 * anchorx - iaxis) * Pnt_Sep /2, (2 * anchory - jaxis) * Pnt_Sep /2, DetectorDepth + RMins[i][j]),
                            RTetLogicVolume,
                            "RTet",         // its name
                            logicWorld,  // its mother  volume
                            false,           // no boolean operations
							100 * i + j);
							// 0);
                };
                
                int n = 0;
                
//                auto check_volume = [&](G4ThreeVector p2, G4ThreeVector p3, G4ThreeVector p4)
//                {
//                    G4double signed_vol=p2.cross(p3).dot(p4);
//
//                    for (int k = 1; k < 4; k++)
//
//                }
                
                // x adjascent
                if  (Relief[i+signx][j] != RMins[i][j])
                {
//                    G4ThreeVector P1 = G4ThreeVector(signx * Pnt_Sep, 0, 0);
//                    G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);
//                    G4ThreeVector P3 = G4ThreeVector(signx * Pnt_Sep, 0, Relief[i+signx][j]-RMins[i][j]);
                    
                    G4ThreeVector P1 = G4ThreeVector(signx * Pnt_Sep, 0, 0);
                    G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[anchorx+signx][anchory+signy]-RMins[i][j]);
                    G4ThreeVector P3 = G4ThreeVector(signx * Pnt_Sep, 0, Relief[anchorx+signx][anchory]-RMins[i][j]);
                    
                    if (P1.cross(P2).dot(P3) != 0.0)
                    {
                        make_a_tet(P1, P2, P3);
                    
                        n++;
                    }
                }
                // y adjascent
                if  (Relief[i][j+signy] != RMins[i][j])
                {
//                    G4ThreeVector P1 = G4ThreeVector(0,signy * Pnt_Sep, 0);
//                    G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);
//                    G4ThreeVector P3 = G4ThreeVector(0, signy * Pnt_Sep, Relief[i][j+signy]-RMins[i][j]);
                    
                    G4ThreeVector P1 = G4ThreeVector(0,signy * Pnt_Sep, 0);
                    G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[anchorx+signx][anchory+signy]-RMins[i][j]);
                    G4ThreeVector P3 = G4ThreeVector(0, signy * Pnt_Sep, Relief[anchorx][anchory+signy]-RMins[i][j]);
                    
                    if (P1.cross(P2).dot(P3) != 0.0)
                    {
                        make_a_tet(P1, P2, P3);
                        
                        n++;
                    }
                }
                // across
                if  (Relief[i+signy][j+signx] != RMins[i][j])
                {
//                    // across and x adjascent
//                    G4ThreeVector P1 = G4ThreeVector(signx * Pnt_Sep, 0, 0);
//                    G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, 0);
//                    G4ThreeVector P3 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);
//
//                    // across and y adjascent
//                    G4ThreeVector P4 = G4ThreeVector(0,signy * Pnt_Sep, 0);
//                    G4ThreeVector P5 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, 0);
//                    G4ThreeVector P6 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);

                    // across and x adjascent
                    G4ThreeVector P1 = G4ThreeVector(signx * Pnt_Sep, 0, 0);
                    G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, 0);
                    G4ThreeVector P3 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[anchorx+signx][anchory+signy]-RMins[i][j]);
                    
                    // across and y adjascent
                    G4ThreeVector P4 = G4ThreeVector(0,signy * Pnt_Sep, 0);
                    G4ThreeVector P5 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, 0);
                    G4ThreeVector P6 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[anchorx+signx][anchory+signy]-RMins[i][j]);
                    
                    if (P1.cross(P2).dot(P3) != 0.0 and P4.cross(P5).dot(P6) != 0.0)
                    {
                        make_a_tet(P1, P2, P3);
                        
                        make_a_tet(P4, P5, P6);
                        
                        n++;
                    }
                }
                // cap
                if (n > 0)
                {
//                    G4ThreeVector P1 = G4ThreeVector(0, signy * Pnt_Sep, Relief[i][j+signy]-RMins[i][j]);
//                    G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, 0, Relief[i+signx][j]-RMins[i][j]);
//                    G4ThreeVector P3 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);
                    
                    G4ThreeVector P1 = G4ThreeVector(0, signy * Pnt_Sep, Relief[anchorx][anchory+signy]-RMins[i][j]);
                    G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, 0, Relief[anchorx+signx][anchory]-RMins[i][j]);
                    G4ThreeVector P3 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[anchorx+signx][anchory+signy]-RMins[i][j]);
                    
                    if (P1.cross(P2).dot(P3) != 0.0)
                    {
                        //make_a_tet(P1, P2, P3);
                    }
                }
            };
            
            // Finds lowest data point to anchor tetrahedra
            if (Relief[i][j] == RMins[i][j])
            {
                    int Anchorx = i; // ----remove anchor varaibles
                    int Anchory = j;
                    int Signx = 1.0;
                    int Signy = 1.0;
                    
                    make_tets(Anchorx, Anchory, Signx, Signy);
            }
            else if (Relief[i+1][j] == RMins[i][j])
            {
                int Anchorx = i + 1;
                int Anchory = j;
                int Signx = - 1.0;
                int Signy = 1.0;
                
                make_tets(Anchorx, Anchory, Signx, Signy);
            }
            else if (Relief[i][j+1] == RMins[i][j])
            {
                int Anchorx = i;
                int Anchory = j + 1;
                int Signx = 1.0;
                int Signy = - 1.0;
                
                make_tets(Anchorx, Anchory, Signx, Signy);
            }
            else
            {
                int Anchorx = i + 1;
                int Anchory = j + 1;
                int Signx = - 1.0;
                int Signy = - 1.0;
                
                make_tets(Anchorx, Anchory, Signx, Signy);
            }/////////////
        }
    }


	//--------- Visualization attributes -------------------------------
	G4VisAttributes* BoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	logicWorld->SetVisAttributes(BoxVisAtt);
	G4Colour col;

	col = G4Colour(0.0, 0.0, 1.0);
	G4VisAttributes* upper_VisAtt = new G4VisAttributes(col);
	logic_upper_plane->SetVisAttributes(upper_VisAtt);

	col = G4Colour(1.0, 0.0, 0.0);
	G4VisAttributes* lower_VisAtt = new G4VisAttributes(col);
	logic_lower_plane->SetVisAttributes(lower_VisAtt);

	col = G4Colour(0.0, 1.0, 0.0);
	G4VisAttributes* Trigger_VisAtt = new G4VisAttributes(col);
	LowerTrigLogicVolume->SetVisAttributes(Trigger_VisAtt);
	UpperTrigLogicVolume->SetVisAttributes(Trigger_VisAtt);
	//RBoxLogicVolume->SetVisAttributes(Trigger_VisAtt);
	return physiWorld;
}
