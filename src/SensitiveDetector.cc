#include "SensitiveDetector.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "DetectorConstruction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4Run.hh"
#include "G4UserRunAction.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4PhysicalConstants.hh"
#include <map>
#include <vector>
extern G4double primary_vertex_energy; //coming from PrimaryGeneratorAction
extern G4int nrofdetectors;
extern FILE * tuples1;
extern FILE * tuples2;
extern G4int intsize;
extern G4int floatsize;
extern G4int ushortsize;
extern G4int ucharsize;
extern G4int mode;
extern G4int nricrubins;
extern G4int nrpdgs;
extern G4double icrustep;
extern G4double mmaxangle;

extern std::vector<G4int> useornot;
extern std::vector<G4double> icru;
extern std::vector<G4double> icru2;
extern std::map<G4int, G4int> heavy;


//std::ofstream G4UserRunAction::outstream;
G4Allocator<SDHit>* SDHit::SDHitAllocator=0;


SDHit::SDHit()
 : G4VHit()
{}

SDHit::~SDHit() {}

void SDHit::Print() {}


SensitiveDetector::SensitiveDetector(const G4String& name, G4bool zero_E)
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL),
   zeroE(zero_E)
{
	collectionName.insert(name);
}

SensitiveDetector::~SensitiveDetector()
{}

void SensitiveDetector::Initialize(G4HCofThisEvent* hce){
	// Create hits collection
	fHitsCollection	= new SDHitsCollection(SensitiveDetectorName, collectionName[0]);
	// Add this collection in hce
	G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	hce->AddHitsCollection( hcID, fHitsCollection );
}

G4bool SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*){
	G4Track *trk   = step->GetTrack();
    	G4double einp = primary_vertex_energy;
	G4double edep  = step->GetTotalEnergyDeposit();
	G4int    detid = stoi(this->GetName());
	if (detid == 0){
		trk->SetTrackStatus(fStopAndKill);
		return false;
	}
	if (detid == (nrofdetectors-1)){
		trk->SetTrackStatus(fStopAndKill);
		return false;
	}
	G4double nioni = step->GetNonIonizingEnergyDeposit();
	G4double ioni  = edep - nioni;
	G4CsvAnalysisManager* analysis = G4CsvAnalysisManager::Instance();
	analysis->FillH2(2,einp,detid,ioni);
	if (mode >= -1){
		// mode -2:      default: only ionization
		// mode -1:      dose rate determination
		// mode 0:   	 output binary data for all atmospheric sheets, at sheet interface
		//         	 purpose: altitude particle profiles
		// mode 1...665: outbut binary data for the bottom interface of the sheet with SDID=mode
		//               purpose: dose rate at a given altitude
		// mode 666: 	 output asci data at each step 
		//           	 purpose: shower visualization
		if (mode > 665) {
			if (mode == 666){
				// It outputs a lot of data in ascii format and can be used to visualize particle showers
				// take care to properly configure the cut in range value for the relevant particles
				const G4Event *evt = G4RunManager::GetRunManager()->GetCurrentEvent();
				const G4Run   *run = G4RunManager::GetRunManager()->GetCurrentRun();
				G4StepPoint* pre = step->GetPreStepPoint();
				const G4StepPoint *post = step->GetPostStepPoint();
				//IDs: Integers
				G4UserRunAction::outstream << run->GetRunID()         << " ";		// RunID
				G4UserRunAction::outstream << evt->GetEventID()       << " ";		// EventID
				G4UserRunAction::outstream << trk->GetParentID()      << " ";		// ParentID
				G4UserRunAction::outstream << trk->GetTrackID()       << " ";		// TrackID
				G4UserRunAction::outstream << detid                   << " ";		// Sensitive Detector ID
				G4UserRunAction::outstream << trk->GetParticleDefinition()->GetPDGEncoding() << " "; // PDG Code
				//Position, momentum and energ
				G4ThreeVector prevec = pre->GetPosition();
				G4ThreeVector postvec = post->GetPosition();
				G4UserRunAction::outstream << prevec.x()/cm << " ";
				G4UserRunAction::outstream << prevec.y()/cm << " ";
				G4UserRunAction::outstream << prevec.z()/cm << " ";
				G4UserRunAction::outstream << postvec.x()/cm << " ";
				G4UserRunAction::outstream << postvec.y()/cm << " ";
				G4UserRunAction::outstream << postvec.z()/cm << " ";
				G4UserRunAction::outstream << pre->GetKineticEnergy()/MeV << " ";
				G4UserRunAction::outstream << step->GetTotalEnergyDeposit()/MeV << " ";
				G4UserRunAction::outstream << step->GetNonIonizingEnergyDeposit()/MeV << " ";
				G4UserRunAction::outstream << post->GetKineticEnergy()/MeV << " ";
				//Information: Strings
				G4UserRunAction::outstream << trk->GetParticleDefinition()->GetParticleName() << " "; // Particle name 
				G4UserRunAction::outstream << post->GetProcessDefinedStep()->GetProcessName() << " "; // Process name
				G4UserRunAction::outstream << G4endl;
			} else {
			G4cout << "### AtRIS: ERROR!! mode can not be above 666! " << G4endl;
			std::cout << "### AtRIS: ERROR!! mode can not be above 666! " << G4endl;
			G4RunManager *runManager = new G4RunManager;
			delete runManager;
			return false;
			}
		}

		//before the step, the particle should have been at the boundary. 
		G4StepPoint* pre = step->GetPreStepPoint();
		if (pre->GetStepStatus() == fGeomBoundary){
			//Let's check which sdid we have:
			G4double ekin = pre->GetKineticEnergy();
			G4ThreeVector position     = pre->GetPosition().unit();
			G4ThreeVector direction    = pre->GetMomentumDirection();
			G4double inclination       = std::acos(position*direction); 
			bool upward                = inclination < halfpi;
			bool withinangle           = ((pi-inclination)<=mmaxangle) || (inclination<=mmaxangle);
			//G4cout << "debuging max=" << mmaxangle << ". inclination=" << inclination/deg << " upward " << upward <<  " result " << withinangle << G4endl;
			G4int diid                 = detid - (1*upward);
			G4int    partPDG           = trk->GetParticleDefinition()->GetPDGEncoding();
			G4int    apartPDG          = abs(partPDG);
			G4int    k = (G4int) (log10(ekin/eV)/icrustep);
			if (k < 0) k=0;
			//G4cout << "PDG=" << partPDG << ", apartPDG=" << apartPDG << ", k=" << k;
			if (withinangle && (apartPDG<2214)) {
				G4int    dpartPDG          = partPDG + 2213;
				bool use = ((useornot[dpartPDG]) or (partPDG==22));
				//G4cout << " is light with use=" << use << ",dpartPDG="<<dpartPDG << ", and ekin=" << ekin/MeV << "MeV";
				if (use) { // light relevant particles
					//G4cout << " and should be used ";
					G4int    j = useornot[dpartPDG]*nricrubins+k;
					// retrieve the corresponding average dose
					G4double f = icru[j]*ekin;
					G4double f2= icru2[j]*ekin; // For the second phantom
					analysis->FillH2(3,einp,diid,f);
					analysis->FillH2(5,einp,diid,f2); // For the second phantom
					// ##########################
					// Extension: equivalent dose
					//            factor = 1 for all particles except:
					//                     2 for protons and pions
					//                    20 for heavy nuclei and fission fragments
					//                       for neutrons, we have a continious function wneut
					G4double wfactor = 1;
					if ((apartPDG == 2212) or (apartPDG == 111)) wfactor=2;
					if (apartPDG == 2112){
						if (ekin <= 1/MeV) {
							wfactor = 3.5 + 18.2*exp(-log(pow(ekin,2.0))/6.0) ;
						} else if (ekin <= 50/MeV) {
							wfactor = 5.0 + 17.0*exp(-log(pow(2*ekin,2.0))/6.0) ;
						} else if (ekin > 50/MeV) {
						        wfactor = 2.5 + 3.25*exp(-log(pow(0.04*ekin,2.0))/6.0) ;
						}
					}
					//G4cout << ", j=" << j << ", f=" << f << ", f2=" << f2 << ", wfactor=" << wfactor;
					analysis->FillH2(4,einp,diid,f*wfactor);
					analysis->FillH2(6,einp,diid,f2*wfactor); // For the second phantom
				}
			} else if (withinangle && (1000010019 < partPDG)) {
				//G4cout << " is heavy ";
				// We first check if we are dealing with one of the smaller, relevant isotopes
				// We then check if we are dealing with one of: H2, H3, He3, He4
				// If so, we add the contribution to the dose-rate histogram
				if (heavy.find(partPDG) != heavy.end()){
					G4int    j = heavy[partPDG]*nricrubins+k;
					G4double f = icru[j]*ekin;
					G4double f2= icru2[j]*ekin;
					//G4cout << j << std::flush << G4endl;
					analysis->FillH2(3,einp,diid,f);
					analysis->FillH2(5,einp,diid,f2); // For the second phantom
					// ##########################
					// Extension: equivalent dose
					//            for heavy nuclei, we have a factor of 20!
					analysis->FillH2(4,einp,diid,f*20);
					analysis->FillH2(6,einp,diid,f2*20); // For the second phantom
					//G4cout << "is heavy with  j=" << j << ", f=" << f << ", f2=" << f2 << ", wfactor=20, and ekin=" << ekin/MeV << "MeV";
				}
			}
			//G4cout << G4endl;
			if (mode >= 0) {
				if (mode == 0) {
				G4float       esource      = (G4float) einp;
				G4float       ekinf         = (G4float) pre->GetKineticEnergy();
				unsigned short entering    = (unsigned short) diid;
				unsigned short sangle        = (unsigned short)round(inclination/deg);
				//G4int is a 4 byte signed integer and this size is necessary to be able to hold ions
				// Here casting is essential, since otherwise fwrite makes wierd stuff.
				//G4cout << "angle before casting " << inclination/deg << G4endl;
				//G4cout << "angle before writing " << angle << " and diid " << diid << G4endl;
				// Now we have everything: PDG code, kinetic energy, inclination, altitude
				// For rescaling we will also include input energy. Because of the variability of input
				// energy, and the sequential nature of the simulation, this is enough info to discriminate 
				// between daughters of different parent particles. 
				fwrite(&partPDG,intsize,1,tuples1);
				fwrite(&entering,ushortsize,1,tuples1);
				fwrite(&sangle,ushortsize,1,tuples1);
				//fwrite(&inclination,floatsize,1,tuples1);
				fwrite(&ekinf,floatsize,1,tuples1);
				fwrite(&esource,floatsize,1,tuples1);
				return true;
			 	} else if (mode == diid) {
				G4float       esource      = (G4float) einp;
				G4float       ekinf         = (G4float) pre->GetKineticEnergy();
				unsigned short entering    = (unsigned short) diid;
				unsigned short sangle        = (unsigned short)round(inclination/deg);
				 //unsigned short entering    = (unsigned short) diid;
				 //unsigned short sangle        = (unsigned short)round(inclination/deg);
				 fwrite(&partPDG,intsize,1,tuples2);
				 fwrite(&entering,ushortsize,1,tuples2);
				 fwrite(&sangle,ushortsize,1,tuples2);
				 fwrite(&ekinf,floatsize,1,tuples2);
				 fwrite(&esource,floatsize,1,tuples2);
				 return true;
				}
			}
		}
	}
	return true;
}


