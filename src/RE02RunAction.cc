#include "RE02RunAction.hh"
//#include "RE02Run.hh"
#include "F04RunActionMessenger.hh"

//-- In order to obtain detector information.
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
//=======================================================================
// RE02RunAction
//
//
//
//=======================================================================
// Constructor
RE02RunAction::RE02RunAction()
{
	runMessenger = new F04RunActionMessenger(this);
}

// Destructor.
RE02RunAction::~RE02RunAction()
{
	delete runMessenger;

}

void RE02RunAction::BeginOfRunAction(const G4Run* aRun)
{

     CLHEP::HepRandom::setTheSeed(seedNum);
     CLHEP::HepRandom::showEngineStatus();
     //read matrix data

}

//
//==
/*G4Run* RE02RunAction::GenerateRun()
{
  // Generate new RUN object, which is specially
  // dedicated for MultiFunctionalDetector scheme.
  //  Detail description can be found in RE02Run.hh/cc.
  return new RE02Run(theSDName);
}

//
//==


//
//==
void RE02RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "entering end of run action" << G4endl;
  //- RE02Run object.
  RE02Run::GetRun()->DumpAllScorer();
  G4cout << "finishing run action" << G4endl;
  //writing an empty file that says "finished"
  std::ofstream outFile_temp("finished");


}*/
