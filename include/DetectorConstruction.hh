#ifndef RE02DetectorConstruction_h
#define RE02DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MultiFunctionalDetector.hh"
#include "CLHEP/Units/SystemOfUnits.h"



class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

//
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
   // constructor and destructor.
   DetectorConstruction();
   virtual ~DetectorConstruction();

public:
  // virtual method from G4VUserDetectorCOnstruction.
  virtual G4VPhysicalVolume* Construct();

public:
  void GetNumberOfSegmentsInPhantom(G4int& nx, G4int& ny, G4int& nz)
     const{ nx=fNx; ny = fNy; nz = fNz; }

  G4double GetZsize(){
	  /*return -World_sizeZ;*/
	  return -300*CLHEP::m;
  }

private:
  // Data members
  #include "UserParameters.hh"
  G4int         fNx,fNy,fNz;    // Number of segmentation of water phantom.

  G4Box*             solidWorld;    // pointer to the solid envelope
  G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
  G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope

  G4Box*             solidShield;    // pointer to the solid envelope
  G4LogicalVolume*   logicShield;    // pointer to the logical envelope
  G4VPhysicalVolume* physiShield;    // pointer to the physical envelope

  G4Box*             solidR;    // pointer to the solid envelope
  G4LogicalVolume*   logicR;    // pointer to the logical envelope
  G4VPhysicalVolume* physiR;    // pointer to the physical envelope

  G4Box*             solidDet;    // pointer to the solid envelope
  G4LogicalVolume*   logicDet;    // pointer to the logical envelope

  G4LogicalVolume**   logicFujiSens;
  G4LogicalVolume** logicFujiA;
  G4LogicalVolume** logYRep;
  G4LogicalVolume** logXRep;
  G4MultiFunctionalDetector** MFDet;
  
  G4LogicalVolume** logicStereoStrip;
  G4LogicalVolume** logicMiddleStrip;
  G4LogicalVolume** logicOuterStrip;
/*
  G4LogicalVolume*   logicFujiSens[NL];
  G4LogicalVolume* logicFujiA[NL];
  G4LogicalVolume* logYRep[NL];
  G4LogicalVolume* logXRep[NL];
  G4MultiFunctionalDetector* MFDet[NL]; */

  G4Material* FujiMatter;
  G4Material* BackMatter;
  G4Material* DetMatter;

  G4double** R;
  G4double** Z;
  G4int NumRho;
  G4int* Num_points;
  G4double* rho;
  G4int* zones;

};
#endif
