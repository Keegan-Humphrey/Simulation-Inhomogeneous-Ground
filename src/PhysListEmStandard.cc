#include "PhysListEmStandard.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuMultipleScattering.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4eMultipleScattering.hh"
//#include "G4VUserPhysicsList.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::PhysListEmStandard(const G4String& name)
	: G4VPhysicsConstructor(name)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::~PhysListEmStandard()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructProcess()
{
	// Add standard EM Processes

	GetParticleIterator()->reset();
	while ((*GetParticleIterator())()) {
		G4ParticleDefinition* particle = GetParticleIterator()->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();

		if (particleName == "gamma") {
			// gamma         
			pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
			pmanager->AddDiscreteProcess(new G4ComptonScattering);
			pmanager->AddDiscreteProcess(new G4GammaConversion);

		}
		else if (particleName == "e-") {
			//electron
			pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4eIonisation, -1, 1, 1);
			pmanager->AddProcess(new G4eBremsstrahlung, -1, 2, 2);

		}
		else if (particleName == "e+") {
			//positron
			pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4eIonisation, -1, 1, 1);
			pmanager->AddProcess(new G4eBremsstrahlung, -1, 2, 2);
			pmanager->AddProcess(new G4eplusAnnihilation, 0, -1, 3);

		}
		else if (particleName == "mu+" ||
			particleName == "mu-") {
			//muon  
			pmanager->AddProcess(new G4MuIonisation, -1, 1, 1);
			pmanager->AddProcess(new G4MuBremsstrahlung, -1, 2, 2);
			pmanager->AddProcess(new G4MuPairProduction, -1, 3, 3);
			pmanager->AddProcess(new G4MuMultipleScattering, -1, 1, 1);

		}
		else if (particleName == "alpha" || particleName == "GenericIon") {
			pmanager->AddProcess(new G4ionIonisation, -1, 1, 1);

		}
		else if ((!particle->IsShortLived()) &&
			(particle->GetPDGCharge() != 0.0) &&
			(particle->GetParticleName() != "chargedgeantino")) {
			//all others charged particles except geantino
			pmanager->AddProcess(new G4hIonisation, -1, 1, 1);
		}
	}//TBD
}
