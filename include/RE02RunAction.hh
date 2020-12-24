#ifndef RE02RunAction_h
#define RE02RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <vector>

class F04RunActionMessenger;
class G4Run;

//=======================================================================
// RE02RunAction
//
//  Generate Run object and Dumping Run summary.
//
//  T.Aso Created. 2007.Nov.
//
//=======================================================================
//
class RE02RunAction : public G4UserRunAction
{
public:
  // constructor and destructor
  RE02RunAction();
  virtual ~RE02RunAction();

public:
  // virtual method from G4UserRunAction.
  //virtual G4Run* GenerateRun();
  virtual void BeginOfRunAction(const G4Run*);
  //virtual void EndOfRunAction(const G4Run*);*/

public:
  void SetSeedNumber(G4int num){
	  seedNum=num;
	  G4cout << "the seed is: " << seedNum << G4endl;
  }



private:
  F04RunActionMessenger* runMessenger;
  G4int seedNum;

};

//

#endif
