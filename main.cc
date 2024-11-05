#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4GDMLParser.hh"
#include "atrisSession.hh"

#include "DetectorConstruction.hh"
#include "UserActionInitialization.hh"
#include "G4UserRunAction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "SensitiveDetector.hh"
#include "Randomize.hh"
#include "G4PhysListFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include <vector>
#include <string>

G4int year = 1987;
G4int mode = -2;
G4int eventGoal = 0;
G4int nrofdetectors = 0;

int main(int argc,char** argv){
	// This part is for forcing the loading of a GDML file
	G4String prefix = G4String (argv[2]);
	if (argc == 4){
		mode = std::stoi(argv[3]);
		G4cout << "### AtRIS: custom mode set: " << mode << G4endl;
	}
	G4cout << "### AtRIS: ./AtRIS QGSP_BERT_HP " << prefix << " " << mode << G4endl;
	G4cout << "### AtRIS: make sure that the following files exist:" << G4endl;
	G4cout << "\t\t "<< prefix <<".gdml\t\t:: geometry specification" << G4endl;
	G4cout << "\t\t "<< prefix <<".mac\t\t:: macro file" << G4endl;
	G4cout << "\t\t "<< prefix <<".bins\t\t:: input spectrum energy bins" << G4endl;
	G4cout << "\t\t "<< "icru_response.txt\t\t:: first default phantom (or alternative)" << G4endl;
	G4cout << "\t\t "<< "cell_response.txt\t\t:: second default phantom (or alternative)" << G4endl;
	G4cout << "### AtRIS: the following files will be overwritten!" << G4endl;
	G4cout << "\t\t "<< prefix <<".ioni\t\t::        mode -2 custom ionization density output" << G4endl;
	G4cout << "\t\t "<< prefix <<".dose\t\t::        mode -1 custom absorbed dose rate output using first phantom" << G4endl;
	G4cout << "\t\t "<< prefix <<".dose2\t\t::       mode -1 custom absorbed dose rate output using second phantom" << G4endl;
	G4cout << "\t\t "<< prefix <<".equi\t\t::        mode -1 custom equivalent dose rate output using first phantom" << G4endl;
	G4cout << "\t\t "<< prefix <<".equi2\t\t::       mode -1 custom equivalent dose rate outputi using second phantom" << G4endl;
	G4cout << "\t\t "<< prefix <<"0.binary\t\t::     mode 0 binary output (all interfaces)" << G4endl;
	G4cout << "\t\t "<< prefix <<"1.binary\t\t::     mode 1-665 binary output (selected interface)" << G4endl;
	G4cout << "\t\t "<< prefix <<"_h1_h1.1.csv\t::   source energy distribution" << G4endl;
	G4cout << "\t\t "<< prefix <<"_h2_h2.dir.csv\t:: source direction distribution" << G4endl;
	G4cout << "\t\t "<< prefix <<"_h2_h2.pos.csv\t:: source position distribution" << G4endl;
	G4cout << "\t\t "<< prefix <<"_h2_h2.ion.csv\t:: resulting ionization histrogram" << G4endl;
	G4cout << "\t\t "<< prefix <<".asci\t\t::        mode 666 asci output" << G4endl;
	G4cout << "\t\t "<< prefix <<".log\t\t::         detailed log file" << G4endl;
	G4cout << "\t\t "<< prefix <<".error\t\t::       error stream file" << G4endl;
	if (argc <3){
		G4cout << "Error! Mandatory input is not specified!" << G4endl;
		return -1;
	}
	// The log files could be necessary for the data analysis.
	// Therefore, we redirect them to a file. 
	G4String erfile = G4String (argv[2]) + ".error";
	G4String logfile = G4String (argv[2]) + ".log";
	G4cout << "### AtRIS: redirecting ouput stream to " << logfile << " for documentation" << G4endl;
	atrisSession *session   = new atrisSession(logfile,erfile);
	G4UImanager * UImanager	= G4UImanager::GetUIpointer();
	UImanager->SetCoutDestination(session);

	G4Random::setTheEngine(new CLHEP::MTwistEngine); // Choose the Random engine
	G4RunManager *runManager = new G4RunManager;     // Instantiate the run manager
	G4GDMLParser parser;                             //instantiate a GDML parser
	parser.Read(G4String (argv[2]) +".gdml");	 //parse the 
	runManager->SetUserInitialization(new DetectorConstruction(parser));
	G4String physics = argv[1]; 			 //Setting up the physics list
	G4PhysListFactory factory;                       // initialize physics list
	if (factory.IsReferencePhysList(physics) or (factory.IsReferencePhysList(physics.substr(0,physics.size()-4)))) {
		G4VModularPhysicsList* physlist = factory.GetReferencePhysList(physics);
		physlist->RegisterPhysics(new G4StepLimiterPhysics());
		runManager->SetUserInitialization(physlist);
	} else {
		G4cout << "### AtRIS: " << physics << " is not valid! use: "<<G4endl;
		std::cout << "### AtRIS: " << physics << " is not valid! use: "<<G4endl;
		std::vector<G4String> available = factory.AvailablePhysLists();
		for (uint i=0; i<available.size(); ++i){
			G4cout << available[i] << G4endl;
			std::cout << available[i] << endl;
		}
		return -1;
	}

	PrimaryGeneratorAction *pga = new PrimaryGeneratorAction;
	runManager->SetUserAction(pga);
	G4UserRunAction *beg = new G4UserRunAction;
	runManager->SetUserAction(beg);
	
	// configure the bin specification for h1.ene and h2.ion
	HistoManager::binedgesFile = G4String (argv[2])+".bins";
	G4cout << "### AtRIS: using the file " << HistoManager::binedgesFile << " to construct the bin edges " << G4endl;
	std::cout << "### AtRIS: using the file " << HistoManager::binedgesFile << " to construct the bin edges " << G4endl;
	HistoManager::loadBins = true;

	// Get the pointer to the user interface manager and start the batch mode
	G4UserRunAction::outputFileName = argv[2];
	G4String tempString             ="/control/execute "+G4String (argv[2])+".mac";
	G4cout << "### AtRIS: Starting the simulation " << G4endl;
	std::cout << "### AtRIS: Starting the simulation " << G4endl;
	UImanager->ApplyCommand(tempString);
	delete runManager;
	delete session;
	G4cout << "### AtRIS: The logfile has been sucessfully closed... " << G4endl;
	return 0;
}
