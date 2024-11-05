#include "atrisSession.hh"
#include "G4UImanager.hh"
atrisSession::atrisSession(G4String outputFileName,G4String errFileName){
	//outFile = ofstream(outputFileName);
	//errFile = ofstream(inputFileName);
	outFile.open(outputFileName);
	errFile.open(errFileName);
	G4UImanager* UI = G4UImanager::GetUIpointer();
	UI->SetCoutDestination(this);
}
atrisSession::~atrisSession(){
	G4UImanager* UI = G4UImanager::GetUIpointer();
	UI->SetCoutDestination(NULL);
	outFile.close();
	errFile.close();
}
G4UIsession* atrisSession::SessionStart(){
	return NULL;
}
G4int atrisSession::ReceiveG4cout(const G4String& output){
	outFile<<output;
	outFile.flush();
	return 0;
}
G4int atrisSession::ReceiveG4cerr(const G4String& err){
	errFile<<err;
	errFile.flush();
	return 0;
}
