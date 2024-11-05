#ifndef G4UserRunAction_h
#define G4UserRunAction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4GenericMessenger.hh"
#include "G4Cache.hh"
#include "G4SystemOfUnits.hh"
#include "G4Types.hh"
#include "G4String.hh"
#include "G4Run.hh"
class HistoManager;

class G4UserRunAction
{
  public:
    G4UserRunAction();
    virtual ~G4UserRunAction();

  public:
    virtual G4Run* GenerateRun();
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);
    void atrisSave();
    G4double		icrumass;
    G4double		geomfact;
    G4double 		getPhantomMass();
    G4String 		phantomFileName;
    static G4String outputFileName;
    static G4String outputFileName2;
    static std::ofstream outstream;
    G4double		icrumass2;
    G4double		geomfact2;
    G4double 		getPhantomMass2();
    G4String 		phantomFileName2;

  private:
    	void DefineCommands();
	HistoManager* fHistoManager;
	G4GenericMessenger*	fMessenger;
	void SetPhantomName(const G4String& newFileName);
	void SetPhantomName2(const G4String& newFileName2);

  protected:
    G4bool isMaster;


  public:
    inline virtual void SetMaster(G4bool val=true)
    { isMaster = val; }
    inline G4bool IsMaster() const
    { return isMaster; }
};

#endif


