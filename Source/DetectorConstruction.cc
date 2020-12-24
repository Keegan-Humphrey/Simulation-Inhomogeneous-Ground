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
G4double nBars = 23;
//G4double Trigger_size = (std::max(BarDepth, (nBars - 1)*BarWidth / 2)) + 5 * cm;
G4double Trigger_size = (std::max(BarDepth, (nBars - 1)*BarWidth / 2)) ;
G4double Trigger_width = 1 * cm;
G4double DetectorDepth = 40.0*m;
G4double UD_Det_centerX = 0*m;
G4double UD_Det_centerY = 0*m;
G4double UD_Det_centerZ = Trigger_width + 2 * BarHight + UD_Det_planes_seperation / 2;
G4double UD_Det_alpha = 0*deg;
G4double UD_Det_beta = 0*deg;
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
	
	//POLYSTYRENE
	G4Material* POLYSTYRENE = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
	 
	
	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;



	////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //--------- Materials of the principal geometrical components (solids)  ---------
	//DetMatter = Delrin;  //amplifier Matter
	DetMatter = POLYSTYRENE;
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
	//G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, GroundMat, "HomoGround");
	G4LogicalVolume* logicGround = new G4LogicalVolume(solidGround, GroundMat, "HomoGround");
	new G4PVPlacement(0,               // no rotation
		G4ThreeVector(0.0, 0.0, DetectorDepth / 2.0), //origin
		logicGround,
		"HomoGround",         // its name
		logicWorld,  // its mother  volume
		false,           // no boolean operations
		0);              // copy number   		                                 

	/* U Box */
	G4double UBoxHight = 2 * m;
    G4VSolid* UBox = new G4Box("UBox", 100*m, 0.5*m, UBoxHight/2);
	//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, Air, "UBox");
	//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, GroundMat, "UBox");
	//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, BMF, "UBox");
	//G4LogicalVolume* UBoxLogicVolume = new G4LogicalVolume(UBox, Vacuum, "UBox");
    
    G4double UboxDepth = 0 * m;
	new G4PVPlacement(0,               // no rotation
		//G4ThreeVector(0.0, 0.0, DetectorDepth / 2 - UBoxHight / 2), //origin
		G4ThreeVector(0.0, 0.0, DetectorDepth / 2 - UBoxHight / 2 - UboxDepth), //origin
		UBoxLogicVolume,
		"UBox",         // its name
		logicGround,  // its mother  volume
		false,           // no boolean operations
		0);


	/* myTorus */
	//G4double myTorusHight = 1 * m;
	//G4VSolid* myTorus = new G4Tubs("myTorus", 1 * m, 2 * m, myTorusHight / 2, 0, twopi);
	//G4LogicalVolume* myTorusLogicVolume = new G4LogicalVolume(myTorus, Air, "myTorus");
	//G4LogicalVolume* myTorusLogicVolume = new G4LogicalVolume(myTorus, GroundMat, "myTorus");
	//G4LogicalVolume* myTorusLogicVolume = new G4LogicalVolume(myTorus, BMF, "myTorus");
	//G4LogicalVolume* myTorusLogicVolume = new G4LogicalVolume(myTorus, Vacuum, "myTorus");
	//new G4PVPlacement(0,               // no rotation
	//	G4ThreeVector(0.0, 10 * m, DetectorDepth / 2 - myTorusHight / 2 - 20 * m),
	//	myTorusLogicVolume,
	//	"myTorus",         // its name
	//	logicGround,  // its mother  volume
	//	false,           // no boolean operations
	//	0);
	

	//    DETECTOR CONTAINER - Air
	//G4double Detector_container_height = 2 * UD_Det_planes_width + UD_Det_planes_seperation;
	G4double Detector_container_height = 4 * BarHight + UD_Det_planes_seperation + 2 * Trigger_width;
	G4double Detector_container_width = Trigger_size;
	G4VSolid* solid_Detector_container = new G4Box("Det_container", Detector_container_width / 2.0, Detector_container_width / 2.0, Detector_container_height / 2.0);
	G4LogicalVolume* logic_Detector_container = new G4LogicalVolume(solid_Detector_container, Vacuum, "Det_container");
	G4ThreeVector Detector_container_center = G4ThreeVector(UD_Det_centerX, UD_Det_centerY, -DetectorDepth / 2.0 + Detector_container_height / 2);
	/*G4RotationMatrix Detector_container_rotm = G4RotationMatrix();
	Detector_container_rotm.rotateX(UD_Det_alpha);
	Detector_container_rotm.rotateY(UD_Det_beta);
	Detector_container_rotm.rotateZ(UD_Det_gamma);
	G4Transform3D transform = G4Transform3D(Detector_container_rotm, Detector_container_center);*/
	new G4PVPlacement(0,
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
	
    
    
    
	//--------- Visualization attributes -------------------------------
	G4VisAttributes* BoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	logicWorld->SetVisAttributes(BoxVisAtt);
	G4Colour col;

	col = G4Colour(0.0, 0.0, 1.0);
	G4VisAttributes* upper_VisAtt = new G4VisAttributes(col);
	logic_upper_plane->SetVisAttributes(upper_VisAtt);

	col = G4Colour(0.0, 0.69, 0.94);
	G4VisAttributes* UBox_VisAtt = new G4VisAttributes(col);
	UBoxLogicVolume->SetVisAttributes(UBox_VisAtt);

	col = G4Colour(0.518, 0.235, 0.047);
	G4VisAttributes* GND_VisAtt = new G4VisAttributes(col);
	logicGround->SetVisAttributes(GND_VisAtt);

	col = G4Colour(1.0, 0.0, 0.0);
	G4VisAttributes* lower_VisAtt = new G4VisAttributes(col);
	logic_lower_plane->SetVisAttributes(lower_VisAtt);

	col = G4Colour(0.0, 1.0, 0.0);
	G4VisAttributes* Trigger_VisAtt = new G4VisAttributes(col);
	LowerTrigLogicVolume->SetVisAttributes(Trigger_VisAtt);
	UpperTrigLogicVolume->SetVisAttributes(Trigger_VisAtt);
	return physiWorld;
}
