#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:Detector(Det)
{ 
  testemDir = new G4UIdirectory("/testem/");
  testemDir->SetGuidance("commands specific to this example");
  
  detDir = new G4UIdirectory("/testem/det/");
  detDir->SetGuidance("detector construction");
        
  MaterCmd = new G4UIcmdWithAString("/testem/det/setMat",this);
  MaterCmd->SetGuidance("Select material of the box.");
  MaterCmd->SetParameterName("choice",false);
  MaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SizeCmd = new G4UIcmdWithADoubleAndUnit("/testem/det/setSize",this);
  SizeCmd->SetGuidance("Set size of the box");
  SizeCmd->SetParameterName("Size",false);
  SizeCmd->SetRange("Size>0.");
  SizeCmd->SetUnitCategory("Length");
  SizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    
  UpdateCmd = new G4UIcmdWithoutParameter("/testem/det/update",this);
  UpdateCmd->SetGuidance("Update calorimeter geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete MaterCmd;
  delete SizeCmd; 
  delete UpdateCmd;
  delete detDir;  
  delete testemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
/*  if( command == MaterCmd )
   { Detector->SetMaterial(newValue);}
   
  if( command == SizeCmd )
   { Detector->SetSize(SizeCmd->GetNewDoubleValue(newValue));}
     
  if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }*/
  G4cout<<"Not updating values Yet."<<G4endl; 
}
