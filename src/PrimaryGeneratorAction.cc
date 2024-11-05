#include "PrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SliceTimer.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include <chrono>
#include <ctime>

typedef std::chrono::duration< double > double_prec_seconds;
//typedef std::chrono::time_point< std::chrono::system_clock,double_prec_seconds > timepoint_t;

extern G4int eventGoal;
G4double primary_vertex_energy; //This is made global, so that it can be used to fill the ion histogram
G4int progress=1;

G4SliceTimer t;
G4SliceTimer* pt = &t;

G4int next = 2;
G4double before = 0;
G4double totaltime = 0;

PrimaryGeneratorAction::PrimaryGeneratorAction(){
	ParticleGenerator	= new G4GeneralParticleSource;
	G4cout << "### AtRIS: Starting the timer " << G4endl;
	std::cout << "### AtRIS: Starting the timer " << G4endl;
	pt->Start();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
	delete ParticleGenerator;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
	G4CsvAnalysisManager* analysis  = G4CsvAnalysisManager::Instance();
	ParticleGenerator->GeneratePrimaryVertex(anEvent);  // We generate the primary particle 

	G4int current  = anEvent->GetEventID();
	if (current % next == 0){
		pt->Stop();
		G4double now = pt->GetRealElapsed();
		G4double thetime = now-before;
		totaltime += thetime;
		G4double speed = progress / thetime;
		G4double tleft = (eventGoal - current)/speed;
		//This code here is only for the purpose of estimating and printing the sim time
		auto now_time = std::chrono::system_clock::now();
		auto eta_time = now_time + std::chrono::duration_cast< std::chrono::seconds >(double_prec_seconds(tleft));
	       	//this is to remove the \n from the ctime print out
		std::time_t eta_c = std::chrono::system_clock::to_time_t(eta_time);
		G4String eta_c2 = ctime(&eta_c);
		eta_c2.pop_back();

		//G4cout << "NOW: " << now_c2 << ":: ETA: " << eta_c2 << " :: " << current << " / " << eventGoal;
		G4cout << "ETA: " << eta_c2 << " :: " << current << " / " << eventGoal;
        	G4cout << ":: Speed:" << speed << "/s. " << G4endl;
		std::cout << "ETA: " << eta_c2 << " :: " << current << " / " << eventGoal;
		std::cout << ":: Speed:" << speed << "/s. " << G4endl;
		G4double timelimit1=100000*s;
		G4double timelimit2=10000*s;
		if ((current > 10000) and (tleft > timelimit1)){
			next = next * 1.2;
			progress = next - next/1.2;
		} else if ((tleft > timelimit2) and (current > 10000)){
			next = next * 1.6;
			progress = next - next/1.6;
		} else {
			next = next * 2;
			progress = next - next/2;
		}
		before = now;
		pt->Start();
		// At this point, we need to apply the histogram modification schemes using g4tools, in order to:
		//    * convert the ionization from:
		//           ***  net energy / complete detector volume (for all particles)
		//           to
		//           ***  eV / ( cm^3 particle)
		analysis->Write();
	}

	// Calculate direction
	G4ThreeVector dir = ParticleGenerator->GetParticleMomentumDirection();
	G4double dirtheta = dir.theta(), dirphi=dir.phi();
	if (dirphi < 0.) dirphi += twopi;
	G4double dircost  = std::cos(dirtheta);

	// Get kinetic energy
	primary_vertex_energy = ParticleGenerator->GetParticleEnergy();

	// Calculate position
	G4ThreeVector pos = ParticleGenerator->GetParticlePosition();
	G4double postheta=pos.theta(), posphi=pos.phi();
	if (posphi < 0.) posphi += twopi;
	G4double poscost = std::cos(postheta);

	// now we fil the histograms
	analysis->FillH1(0,primary_vertex_energy,primary_vertex_energy);
	analysis->FillH2(0,dirphi,dircost);
	analysis->FillH2(1,posphi,poscost);
}

