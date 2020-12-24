#include "globals.hh"
#include "Randomize.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"

#include "RE02RunAction.hh"
#include "F04RunActionMessenger.hh"

F04RunActionMessenger::F04RunActionMessenger(RE02RunAction* RA)
  : runAction (RA)
{
  FileDir2 = new G4UIdirectory("/SetSeed/");
  FileDir2->SetGuidance("set a seed.");

  seedCmd = new G4UIcmdWithAnInteger("/SetSeed/set",this);
  seedCmd->SetGuidance("Set the seed");
  seedCmd->SetParameterName("seed",true);
  seedCmd->SetDefaultValue (1);
  //seedCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

F04RunActionMessenger::~F04RunActionMessenger()
{

}

void F04RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == seedCmd)
  {
	 G4cout << "the seed is: " << newValue << G4endl;

     runAction->SetSeedNumber(seedCmd->GetNewIntValue(newValue));
  }
}
