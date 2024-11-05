#ifndef LoggedSession_hh
#define LoggedSession_hh

#include "G4UIsession.hh"
#include <iostream>
#include <fstream>
using namespace std;

class atrisSession : public G4UIsession {
public:
	atrisSession(G4String,G4String);
	~atrisSession();
	G4UIsession* SessionStart();
	G4int ReceiveG4cout(const G4String&);
	G4int ReceiveG4cerr(const G4String&);
private:
	ofstream outFile;
	ofstream errFile;
};
#endif


