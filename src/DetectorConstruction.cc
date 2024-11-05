#include "DetectorConstruction.hh"
#include "G4UserRunAction.hh"
#include "SensitiveDetector.hh"
#include "G4GDMLParser.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include <map>
#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4ApplicationState.hh"


// These vectors hold data necessary to perform proper scaling of ionization and dose rate histograms
std::vector<G4double> sdid_volume; // used for normalization in G4UserRunAction
std::vector<G4double> sdid_mass;
std::vector<G4double> sdid_rmin;
std::vector<G4double> diid_icru_factor;
std::vector<G4double> diid_icru_factor2;
G4double mmaxangle;

DetectorConstruction::DetectorConstruction(const G4GDMLParser& parser):
       	G4VUserDetectorConstruction(),
	fParser(parser),
       	fMessenger(nullptr){
	DefineCommands();
	icrusurface  = 4*15*15*pi*cm2;
	icrusurface2 = 4*0.05*0.05*pi*mm2;
	mmaxangle = 90*deg;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	  return fParser.GetWorldVolume();
}


/*DetectorConstruction::~DetectorConstruction() {
	;
}*/

void DetectorConstruction::ConstructSDandField() {
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	// Dictionaries to hold basic volume information, before restructuring in to vector
	// This is necessary since the gdml parsing is not in-order. Sometimes detectors with 
	// higher SDIDs are parsed before others with lower. If we were to construct a vector
	// to hold this data immidiately, we would not have the required order for later scale.
	std::map<G4int, G4double> SDIDmass;
	std::map<G4int, G4double> SDIDvolume;
	std::map<G4int, G4double> SDIDsurface;
	std::map<G4int, G4double> SDIDdensity;
	std::map<G4int, G4double> DIIDsurface;
	std::map<G4int, G4double> ICRUfactor;
	std::map<G4int, G4double> ICRUfactor2;

	const G4GDMLAuxMapType* auxmap = fParser.GetAuxMap();
	for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin(); iter!=auxmap->end(); iter++) {
		for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin(); vit!=(*iter).second.end(); vit++) {
			// This iteration is not in order! 
			G4cout << "--> Type: " << (*vit).type << " Value: " << (*vit).value << ", ";
		}
	}
	G4cout << G4endl;
	// We now generate the sensitive detectors and assign them to the proper volumes
	// We also count the number of generated sensitive detectors
	extern G4int nrofdetectors;
	for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin(); iter!=auxmap->end(); iter++) {
		for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin(); vit!=(*iter).second.end();vit++) {
			if ((*vit).type=="SensDet") {
				 nrofdetectors += 1;
				 G4cout << "SD " << (*vit).value << " to " << ((*iter).first)->GetName() << ", ";
				 G4String detname            = (*vit).value;
				 G4int SDID                  = std::stoi(detname);
				 SensitiveDetector* detSD    = new SensitiveDetector(detname);
				 SDman->AddNewDetector( detSD);
				 G4VSensitiveDetector* mydet = SDman->FindSensitiveDetector((*vit).value);
				 if(mydet) {
					 G4LogicalVolume* myvol = (*iter).first;
					 myvol->SetSensitiveDetector(mydet);
					 // At this point we implement the calculation of volume, mass, surface
					 G4double mass     = myvol->GetMass();
					 G4VSolid* mysol   = myvol->GetSolid();
					 G4double surface  = mysol->GetSurfaceArea();
					 G4double volume   = mysol->GetCubicVolume();
					 G4double density  = mass / volume;
					 // Some debugging output:
					 G4cout << "m= " << mass/kg << " kg. A= " << surface/cm2 << " cm^2. V= " << volume/cm3 << " cm^3. density=" << density/(g/cm3) << G4endl;
					 // populating the dictionaries
					 SDIDmass[SDID]    = mass;
					 SDIDvolume[SDID]  = volume;
					 SDIDsurface[SDID] = surface;
					 SDIDdensity[SDID] = density;
				 }
				 else {
					 G4cout << (*vit).value << " detector not found" << G4endl;
				 }
			}
		}
	}
	//At this point we need to calculate and populate the DIID surfaces. 
	G4cout << "calculating DIID surfaces using a phantom surface of " << icrusurface/cm2 << "cm2, and a second phantom surface of "<< icrusurface2/cm2 << "cm2 while using a maxangle of " << mmaxangle/deg << "degrees " << G4endl;
	DIIDsurface[0] = SDIDsurface[0];
	G4cout << " set initial DIIDsurface with index 0 to: " << DIIDsurface[0]/cm2 << "cm2"  << G4endl;
	ICRUfactor[0]  = icrusurface/SDIDsurface[0];
	ICRUfactor2[0] = icrusurface2/SDIDsurface[0];
	for (int i=1; i<nrofdetectors; i++){
		G4cout << "calculating with SDIDsurface" << SDIDsurface[i] << " DIIDsurface:" << DIIDsurface[i-1] << " ";
		G4cout << "DIID: " << i << " has surface: " ;
		G4double surfdiff = (SDIDsurface[i] - DIIDsurface[i-1]);
		DIIDsurface[i]    = surfdiff;
		ICRUfactor[i]     = icrusurface  / (surfdiff*2);
		ICRUfactor2[i]    = icrusurface2 / (surfdiff*2);
		G4cout << surfdiff/cm2 << " and an ICRUfactor of " << ICRUfactor[i] << " and a second ICRUfactor2 of " << ICRUfactor2[i] << G4endl;
		// these so defined icrufactors make sense only for planets. Ignore when simulating the sphere itself.
	}
	// It makes sence to convert now the dictionaries to vectors or arrays.
	// We had to use dictionaries because the GDML parsing is not in order. Now that we know the SDIDs,
	// we can construct an orderly data type containing the relevant data.
	sdid_rmin.push_back(0);
	for (int i=0; i<nrofdetectors; i++){
		diid_icru_factor.push_back(ICRUfactor[i]);
		diid_icru_factor2.push_back(ICRUfactor2[i]);
		sdid_rmin.push_back(sqrt(DIIDsurface[i] / (4.*pi)));
		sdid_mass.push_back(SDIDmass[i]);
		sdid_volume.push_back(SDIDvolume[i]);
	}
	// Debugging: printout the new vector data:
	G4cout << "### AtRIS: restructuring of dictionaries into vectors has happened " << G4endl;
	for (int i=0; i<nrofdetectors; i++){
		G4cout << "SDID=" << i << " V=" << sdid_volume[i]/cm3 << "cm^3 m=" <<sdid_mass[i]/kg << "kg f=" << diid_icru_factor[i] << G4endl;
	}
	G4cout << "### AtRIS: printing detailed information about all materials" << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
	G4cout << "### AtRIS: gdml construction complete. In total we have " << nrofdetectors << " sensitive detectors" << G4endl;
	std::cout << "### AtRIS: gdml construction complete. In total we have " << nrofdetectors << " sensitive detectors" << G4endl;
}



void DetectorConstruction::DefineCommands() {
	fMessenger  = new G4GenericMessenger(this,"/atris/","atris custom interface for configuration");
	G4GenericMessenger::Command& phantomArea = fMessenger->DeclarePropertyWithUnit("setPhantomSurface","centimeter2",icrusurface,"surface of the radiation phantom");
	phantomArea.SetParameterName("A", true);
	phantomArea.SetRange("A>0.");
	phantomArea.SetDefaultValue("2827.433388");
	// same for the second phantom.
	G4GenericMessenger::Command& phantomArea2 = fMessenger->DeclarePropertyWithUnit("setPhantom2Surface","millimeter2",icrusurface2,"surface of the second radiation phantom");
	phantomArea2.SetParameterName("A", true);
	phantomArea2.SetRange("A>0.");
	phantomArea2.SetDefaultValue("0.031416");
	G4GenericMessenger::Command& maxa = fMessenger->DeclarePropertyWithUnit("setMaxAngle","radian",mmaxangle,"max particle angle for the dose rate calculation");
	maxa.SetParameterName("g", true);
	maxa.SetRange("g>0.");
	maxa.SetDefaultValue("1.5707963268");
}

DetectorConstruction::~DetectorConstruction(){
	delete fMessenger;
}

G4double DetectorConstruction::getPhantomArea(){
	return DetectorConstruction::icrusurface;
}
G4double DetectorConstruction::getPhantomArea2(){
	return DetectorConstruction::icrusurface2;
}
