#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

#define Nangles 10000
#define NanglesInv 1500
#define Nenergies 4990
#define NenergiesInv 30000
#define Emax 1000
#define MuonsRestMass 0.10566

class G4Event;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(DetectorConstruction*);    
   ~PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);
    G4ParticleGun* GetParticleGun() {return particleGun;};

  private:
    G4ParticleGun*        particleGun;
    DetectorConstruction* detector;
    G4double Ftheta[Nangles];
    G4double Phi[Nangles][Nenergies];
    G4double cost[Nangles];
    G4double E[Nenergies];    
    G4double Inverse_Theta[NanglesInv];
    G4double PhiInterMed[NanglesInv][Nenergies];
    G4double Inverse_Phi[NanglesInv][NenergiesInv];
    G4double StoreE[10000000];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
