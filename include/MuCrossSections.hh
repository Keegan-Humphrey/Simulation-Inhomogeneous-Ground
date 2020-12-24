#ifndef MuCrossSections_h
#define MuCrossSections_h 1

#include "globals.hh"

class G4Material;
class G4Element;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MuCrossSections
{
  public:
    MuCrossSections();
   ~MuCrossSections();

  public:
    G4double CR_Macroscopic (const G4String&, G4Material*, G4double, G4double);   
    G4double CR_PerAtom     (const G4String&, G4Element* , G4double, G4double);
                       
  private:
    G4double CRB_Mephi (G4double, G4double, G4double, G4double);
    G4double CRK_Mephi (G4double, G4double, G4double, G4double);
    G4double CRN_Mephi (G4double, G4double, G4double, G4double);
    G4double CRP_Mephi (G4double, G4double, G4double, G4double);    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
