#ifndef F04RunActionMessenger_h
#define F04RunActionMessenger_h 1

#include "globals.hh"

#include "G4UImessenger.hh"

class RE02RunAction;

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithABool;

class F04RunActionMessenger : public G4UImessenger
{
  public:

    F04RunActionMessenger(RE02RunAction* );
    ~F04RunActionMessenger();

    void SetNewValue(G4UIcommand* ,G4String);

  private:

    RE02RunAction*             runAction;

    G4UIdirectory*             FileDir2;
    G4UIcmdWithAnInteger*      seedCmd;

};

#endif
