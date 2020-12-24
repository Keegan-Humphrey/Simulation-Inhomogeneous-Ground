#include "PrimaryGeneratorAction.hh"

#include "DetectorConstruction.hh"
#include "Relief.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"

#define M_PI       3.1415926535897931160E0
#define M_PI_2     1.5707963267948965580E0



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern unsigned long cnt_total_events;
extern G4double UD_Det_centerZ;
extern G4double UD_Det_centerY;
extern G4double UD_Det_centerX;
extern G4double UD_Det_planes_seperation;
extern G4double Trigger_size;
extern G4double DetectorDepth;
extern G4double UD_World_sizeZ;


//#define Nstore2 10000001
//extern G4double StoreT[Nstore2];
PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
	:detector(det)
{
	particleGun = new G4ParticleGun(1);
	G4ParticleDefinition* particle
		= G4ParticleTable::GetParticleTable()->FindParticle("mu-");
	particleGun->SetParticleDefinition(particle);
	particleGun->SetParticleEnergy(10 * TeV);
	particleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));

	//G4double Emin = 40;
	G4double Emin = DetectorDepth * 0.4; // Muons lose ~0.6GeV per meter in the ground, so muons with less than DetectorDepth*0.4GeV are not likely to reach the detector

	//prepare density of states:
	G4double dc = 1.0 / ((double)(Nangles - 1));
	G4double Ar = 0.00253;
	G4double a0 = 0.2455;
	G4double a1 = 1.288;
	G4double a2 = -0.2555;
	G4double a3 = 0.0209;

	G4double diffE = ((double)Emax - (double)Emin) / ((double)(Nenergies - 1));

	/* simply the formula of Rayna */
	for (G4int k = 0; k < Nangles; k++) {
		// cos(theta) ranges:     0.0 <   cost    < 1.0  
		cost[k] = ((double)k)*dc;
		//G4cout << cost[k] <<G4endl;
		for (G4int Ei = 0; Ei < Nenergies; Ei++) {
			E[Ei] = Emin + Ei*diffE;
			G4double p1 = std::sqrt(E[Ei] * E[Ei] - MuonsRestMass*MuonsRestMass);
			G4double p = p1*cost[k];
			G4double y = std::log10(p);
			G4double Phi0 = Ar*std::pow(p, -(a3*y*y*y + a2*y*y + a1*y + a0));
			Phi[k][Ei] = cost[k] * cost[k] * cost[k] * Phi0;
			Phi[k][Ei] *= cost[k]; /* ARTIFICIAL CORRECTION FOR VIRTUAL PLANE!!!*/
		}
	}

	/* Generating the functions of \theta etc */
	/* Generating Ftheta - Storing it in 'cumFtheta'*/
	for (G4int k = 0; k < Nangles; k++) {
		Ftheta[k] = 0.0;
		for (G4int Ei = 0; Ei < Nenergies; Ei++) {
			Ftheta[k] += Phi[k][Ei] * diffE;
		}
	}

	/* Generating the Cumulative functions of Theta - and counting total flux */
	Ftheta[0] = 0;
	for (G4int k = 1; k < Nangles; k++) {
		Ftheta[k] += Ftheta[k - 1];
	}

	G4double Total_flux = Ftheta[Nangles - 1] * dc * 2 * M_PI;
	for (G4int k = 0; k < Nangles; k++) {
		Ftheta[k] /= Total_flux / (dc * 2 * M_PI);
	}

	//std::ofstream CTL("Cumulative_Theta.out");
	//CTL.precision(16);
	//for (G4int k=0;k<Nangles;k++){
	//	CTL << cost[k] << "    " << Ftheta[k]   <<"    " << G4endl;
	//}

	/* Building The Inverse function*/
	G4double dU = 1.0 / (NanglesInv - 1);
	for (G4int vXi = 0; vXi < NanglesInv; vXi++) {
		G4double Cur_value = vXi*dU;
		G4int k = 0;
		while (Cur_value > Ftheta[k]) { /* The loop break at cumFtheta[Nangles-1]==1.0, but if cumFtheta[Nangles-1]<1 it will break due to the if()*/
			k++;
			if (k == Nangles) { break; }
		}
		if (k == 0) { Inverse_Theta[vXi] = 0.0; continue; }
		if (k == Nangles) { Inverse_Theta[vXi] = 1.0; continue; }
		Inverse_Theta[vXi] = cost[k - 1] + (cost[k] - cost[k - 1]) / (Ftheta[k] - Ftheta[k - 1])*(Cur_value - Ftheta[k - 1]);
	}

	//std::ofstream CT("Cumulative_ThetaInv.out");
	//CT.precision(16);
	//for (G4int k=0;k<NanglesInv;k++){
	//		CT<< k*dU <<"   "<<  Inverse_Theta[k]   <<G4endl;
	//}

	for (G4int vCi = 0; vCi < NanglesInv; vCi++) {
		G4double cosTheta = Inverse_Theta[vCi];
		for (G4int Ei = 0; Ei < Nenergies; Ei++) {
			G4double p1 = std::sqrt(E[Ei] * E[Ei] - MuonsRestMass*MuonsRestMass);
			G4double p = p1*cosTheta;
			G4double y = std::log10(p);
			G4double Phi0 = Ar*std::pow(p, -(a3*y*y*y + a2*y*y + a1*y + a0));
			PhiInterMed[vCi][Ei] = cosTheta*cosTheta*cosTheta*Phi0;
			PhiInterMed[vCi][Ei] = cosTheta*PhiInterMed[vCi][Ei]; /* ARTIFICIAL CORRECTION FOR VIRTUAL PLANE!!!*/
		}
	}

	/* Generating the Cumulatize functions of Theta Energy - and counting total flux */
	for (G4int vCi = 0; vCi < NanglesInv; vCi++) {
		PhiInterMed[vCi][0] = 0;
		for (G4int Ei = 1; Ei < Nenergies; Ei++) {
			PhiInterMed[vCi][Ei] += PhiInterMed[vCi][Ei - 1];
		}
		for (G4int Ei = 0; Ei < Nenergies; Ei++) {
			PhiInterMed[vCi][Ei] /= PhiInterMed[vCi][Nenergies - 1];
		}
	}

	//std::ofstream CTE("Cumulative_ThetaEnergy.out");
	//for (G4int k=0;k<Nangles;k++){
	//	for (G4int Ei=0; Ei<Nenergies; Ei++){
	//		CTE<<  PhiInterMed[NanglesInv-1][Ei]   <<"    " ;
	//	}  
	//	CTE << G4endl;
	//}

	G4double dUE = 1.0 / (NenergiesInv - 1);
	for (G4int vCi = 0; vCi < NanglesInv; vCi++) {
		for (G4int vEi = 0; vEi < NenergiesInv; vEi++) {
			G4double Cur_value = vEi*dUE;
			G4int k = 0;
			while (Cur_value > PhiInterMed[vCi][k]) {
				k++;
				if (k == Nenergies) { break; };
			}
			if (k == 0) {
				Inverse_Phi[vCi][vEi] = Emin;
				continue;
			}
			if (k == Nenergies) {
				Inverse_Phi[vCi][vEi] = Emax;
				continue;
			}
			Inverse_Phi[vCi][vEi] = E[k - 1] + (E[k] - E[k - 1]) / (PhiInterMed[vCi][k] - PhiInterMed[vCi][k - 1])*(Cur_value - PhiInterMed[vCi][k - 1]);
		}
	}

	//std::ofstream CTR("Cumulative_PhiInv.out");
	//for (G4int k=0;k<NanglesInv;k++){
	//	for (G4int vEi=0;vEi<NenergiesInv;vEi++){
	//		CTR<<  Inverse_Phi[NanglesInv-1][vEi]<<"   ";
	//	}
	//	CTR<<G4endl;
	//}


	std::ofstream inputCurrent("inputCurrent");
	inputCurrent << "total number of mu- above energy " << Emin << " GeV is " << Total_flux << " number/cm2/sec" << G4endl;
	inputCurrent << "per 1 hour " << Total_flux*60.0*60.0 << " number/cm2" << G4endl;
	inputCurrent << "per 1 hour and size (300*m)^2 " << Total_flux*60.0*60.0*30000.0*30000.0 << " number" << G4endl;



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
	//this function is called at the begining of event   
	bool accept_event = 0; //if the event that is drawn is not likely to reach the detector kill it but count it!
	//G4double L = 300.0*m;
	G4double L = Trigger_size*1.2;
	//G4double zInitial = 21.0*m;
	//G4double zInitial = 5.5*m;
	G4double x, y, momx, momy, momz, energy;


	while (accept_event == 0) {
		//position  
		x = -L / 2.0 + G4UniformRand()*L;
		y = -L / 2.0 + G4UniformRand()*L;

		//x = 1.2*cm;
		//y = 1.2*cm;

		//theta
		G4double rnd1 = G4UniformRand();
		G4double rnd1_index = (rnd1)*((G4double)(NanglesInv - 1));
		G4int rnd1_floor = (G4int)floor(rnd1_index);

		/*INTERPOLATION*/
		//G4double costheta = 1.0;
		G4double costheta;
		if (rnd1_floor == NanglesInv - 1) {
			costheta = 1.0;
		}
		else {
			costheta = Inverse_Theta[rnd1_floor] + (Inverse_Theta[rnd1_floor + 1] - Inverse_Theta[rnd1_floor])*(rnd1_index - (G4double)rnd1_floor);
		}

		//G4double Phi = 0.0;
		G4double Phi = 2 * pi*G4UniformRand();
		momz = -costheta;
		momx = std::sqrt(1.0 - costheta*costheta)*cos(Phi);
		momy = std::sqrt(1.0 - costheta*costheta)*sin(Phi);
		//energy- based on k:
		G4double rnd2 = G4UniformRand();
		G4double rnd2_index = (rnd2)*((G4double)(NenergiesInv - 1));
		G4int rnd2_floor = (G4int)floor(rnd2_index);

		/*INTERPOLATION*/
		if (rnd2_floor == NenergiesInv - 1) {
			energy = Emax;
		}
		else {
			energy = Inverse_Phi[rnd1_floor][rnd2_floor] + (Inverse_Phi[rnd1_floor][rnd2_floor + 1] - Inverse_Phi[rnd1_floor][rnd2_floor])*(rnd2_index - (G4double)rnd2_floor);
		}

		//energy = 10;

		G4double r = UD_Det_planes_seperation / costheta;
		G4double DownX = momx*r + x;
		G4double DownY = momy*r + y;
		if (abs(DownX) <= L/2 && abs(DownY) <= L / 2)
		{
			r = (DetectorDepth + max_landscape_height * m - UD_Det_planes_seperation - 1 * cm) / costheta;
            //r = (UD_World_sizeZ * 0.5 - UD_Det_planes_seperation - 1 * cm) / costheta;
			x = -r*momx + x + UD_Det_centerX;
			y = -r*momy + y + UD_Det_centerY;
			accept_event = 1;
		}


		/*
		// Calculating crossing of the muon ray with a blocking sphere of the detecotr region
		G4double RequireR = 0.3*m;
		G4double o1 = x - UD_Det_centerX;
		G4double o2 = y - UD_Det_centerY;
		G4double o3 = zInitial - UD_Det_centerZ;
		G4double A = 1.0;
		G4double B = 2.0*(o1*momx + o2*momy + o3*momz);
		G4double C = (o1*o1 + o2*o2 + o3*o3) - RequireR*RequireR;
		G4double determinant = B*B - 4 * A*C;
		if (determinant > 0) {
			accept_event = 1;
		}


		*/
		cnt_total_events++;
	}

	//StoreT[cnt_total_events]=energy;
	//G4cout << "Particle lanuched. Energy=  " << energy << ".   momz=" << momz << G4endl;

	//SEND THE PARTICLE
	particleGun->SetParticlePosition(G4ThreeVector(x, y, DetectorDepth));
	particleGun->SetParticleEnergy(energy*GeV);
	particleGun->SetParticleMomentumDirection(G4ThreeVector(momx, momy, momz));
	particleGun->GeneratePrimaryVertex(anEvent);
}
