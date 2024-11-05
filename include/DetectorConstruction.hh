#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "SensitiveDetector.hh"

class G4Material;
class G4GDMLParser;
//some of the above can probably be deleted


class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    	DetectorConstruction(const G4GDMLParser& parser);
    	~DetectorConstruction();

    	virtual void ConstructSDandField();
    	virtual G4VPhysicalVolume* Construct();
	G4double		icrusurface;
	G4double		icrusurface2;

	G4double 		getPhantomArea();
	G4double 		getPhantomArea2();


private:
    	void DefineCommands();
	const G4GDMLParser& fParser;

	G4GenericMessenger*	fMessenger;
};

#endif
