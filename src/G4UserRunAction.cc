#include "G4UserRunAction.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DetectorConstruction.hh"
#include "G4GenericMessenger.hh"
#include "G4UnitsTable.hh"
#include "HistoManager.hh"
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

extern G4int eventGoal;
extern G4int nrofdetectors;
extern G4int mode;
extern std::vector<G4double> enbins;
G4int intsize = sizeof(G4int);
G4int floatsize = sizeof(G4float);
G4int ushortsize = sizeof(unsigned short);
G4int ucharsize  = sizeof(unsigned char);
G4String G4UserRunAction::outputFileName;
G4String G4UserRunAction::outputFileName2;
std::ofstream G4UserRunAction::outstream;
FILE * tuples1;
FILE * tuples2;

extern std::vector<G4double> sdid_volume;
extern std::vector<G4double> sdid_mass;
extern std::vector<G4double> sdid_rmin;
extern std::vector<G4double> diid_icru_factor;
extern std::vector<G4double> diid_icru_factor2;
extern G4double mmaxangle;

G4int nricrubins = 650;
G4int nrpdgs;
G4double icrustep = 0.02;
std::vector<G4int> useornot(4427,0);
std::vector<G4double> icru;
std::vector<G4double> icru2;
std::vector<G4int> pdgs;
std::map<G4int, G4int> heavy;

G4UserRunAction::G4UserRunAction()
:fHistoManager(0), fMessenger(nullptr),isMaster(true) 
{
 DefineCommands();
 icrumass = 14137.1669*g;               // default phantom mass
 icrumass2= 0.00052208*mg;              // default second phantom mass
 phantomFileName  = "icru_response.txt"; // default phantom file name
 phantomFileName2 = "cell_response.txt"; // default second phantom file name
 if(!(G4ParticleTable::GetParticleTable()->GetReadiness()))
 {
   G4String msg =  " You are instantiating G4UserRunAction BEFORE your G4VUserPhysicsList \n";
   G4Exception("G4UserRunAction::G4UserRunAction()","Run0041",FatalException,msg);
 }
}

G4UserRunAction::~G4UserRunAction(){
   delete fMessenger;
   delete G4CsvAnalysisManager::Instance();
}

G4Run* G4UserRunAction::GenerateRun()
{ return 0; }
void G4UserRunAction::BeginOfRunAction(const G4Run* aRun){
	G4cout << "### AtRIS: executing BeginOfRunAction" << G4endl;
	std::cout << "### AtRIS: executing BeginOfRunAction" << G4endl;
    	eventGoal = ((G4Run *)(aRun))->GetNumberOfEventToBeProcessed();
	if (mode >= -1){
		// We process icru_response.txt
		// Both phantoms must include same pdgs!
		G4cout << "### AtRIS: loading ICRU response matrix from filename: " << phantomFileName <<  "- Found pdgs: " << G4endl << std::flush;
		std::ifstream is(phantomFileName);
		G4int pdgfound;
		if (is.is_open()){
			G4String bLine;         // instantiate the header string
			G4String temp;
			std::getline(is,bLine); // read first line
			std::stringstream ss;    // create a string stream
			ss << bLine;            // put the header in to the string stream
			G4int    j = 0;             // initialize the column index
			while (!ss.eof()) {
				ss >> temp;
				if (std::stringstream(temp) >> pdgfound){
					pdgs.push_back(pdgfound);
					if (abs(pdgfound)<2214){
						G4cout << j << " adding sensitivity for light particle with PDG= " << pdgfound << " "<<  "  at index " << pdgfound + 2213 << G4endl;
					       	useornot[pdgfound+2213] = j; // modify the corresponding useornot location to j
						// this offset is for antiparticles up to neutron and proton!
					}else {
						G4cout << "Added sensitivity for heavy particle with PDG= " << pdgfound << " under the index j=" << j << G4endl;
						heavy[pdgfound] = j;
					}
					j++;                         // increment j
				}
			}

			G4cout << "total" << j << G4endl;
			nrpdgs = j;
			// Now we need to process all the data and put it inside the icru vector
			G4int expected = nricrubins * (j);
			G4double f;
			G4int count = 0;
		        while ( std::getline(is,bLine)){
				std::stringstream dd;    // create a string stream
				dd << bLine;
				while (!dd.eof()){
					dd >> temp;
					if (std::stringstream(temp) >> f){
					       icru.push_back(f);
					       count++;
					}
				}
			}	
			G4cout << "### AtRIS: found " << count << " out of expected "<< expected  << G4endl;
			std::cout << "### AtRIS: found " << count << " out of expected "<< expected  << std::endl;
			if (count != expected) G4cout << "check phantoms file for empty ending charachters!!" << G4endl;
		} else {
			G4cout << "### AtRIS: couldn't find ICRU response matrix for filename" << phantomFileName << G4endl;
		}
		G4cout << "### AtRIS: loading second ICRU response matrix from filename: " << phantomFileName2 <<  "-for previously determined number of pdgs: " << nrpdgs << G4endl << std::flush;
		std::ifstream is2(phantomFileName2);
		if (is2.is_open()){
			// We skip the first line, and process the rest
			G4String bLine;         // instantiate the header string
			G4String temp;
			G4String dummyLine;
			std::getline(is2, dummyLine);
			G4int expected = nricrubins * (nrpdgs);
			G4double f;
			G4int count = 0;
		        while ( std::getline(is2,bLine)){
				std::stringstream dd;    // create a string stream
				dd << bLine;
				while (!dd.eof()){
					dd >> temp;
					if (std::stringstream(temp) >> f){
					       icru2.push_back(f);
					       count++;
					}
				}
			}	
			G4cout << "### AtRIS: found " << count << " out of expected "<< expected  << G4endl;
			std::cout << "### AtRIS: found " << count << " out of expected "<< expected  << std::endl;
		} else {
			G4cout << "### AtRIS: couldn't find second response matrix for filename " << phantomFileName2  << G4endl;
		}

		if (mode >= 0){
			G4String bFN = G4UserRunAction::outputFileName + (G4String) "1.binary";
			G4cout << "### AtRIS: opening file:" << bFN;
			G4cout << " for binary output with floatsize = " << floatsize;
			G4cout << " and intsize = " << intsize << G4endl;
			G4cout << " and unsigned shortsize= " << ushortsize << G4endl;
			tuples1 = std::fopen(bFN,"wb");
			if (mode > 0){
				G4String bFN2 = G4UserRunAction::outputFileName + (G4String) "2.binary";
				tuples2       = std::fopen(bFN2,"wb");
			}
			if (mode == 666){
				G4String bFN3 = G4UserRunAction::outputFileName + (G4String) ".asci";
				G4UserRunAction::outstream.close();
				G4UserRunAction::outstream.open(bFN3, std::ofstream::out);
				G4String header = "#runID eventID parentID trackID detectorID codePDG startX startY startZ endX endY endZ startE eDep ioniDep nonioniDep stopE partName processName";
				G4UserRunAction::outstream << header << G4endl;
			}
		}
	}
	fHistoManager = new HistoManager();
	G4CsvAnalysisManager * analysis = G4CsvAnalysisManager::Instance();
	if ( analysis->IsActive() ) {
		analysis->SetFileName(G4UserRunAction::outputFileName);
		analysis->OpenFile();
	}
}

void G4UserRunAction::EndOfRunAction(const G4Run*){
  G4cout << "### AtRIS: now writing out histograms... " << G4endl;
  std::cout << "### AtRIS: now writing out histograms... " << G4endl;
  // At this point we implement the scaling according to mass, volume, surface area
  atrisSave();
  auto analysis = G4CsvAnalysisManager::Instance();
  if ( analysis->IsActive() ) {
	  analysis->SetVerboseLevel(1);
	  analysis->Write();
	  analysis->CloseFile();
  }
  if (G4UserRunAction::outstream.is_open()){
	  G4cout << "### AtRIS: simulation completed sucessfully!" << G4endl;
	  std::cout << "### AtRIS: simulation completed sucessfully!" << G4endl;
	  G4UserRunAction::outstream.close();
  }
  if (mode == 0) {std::fclose(tuples1);}
  if (mode > 0)	{std::fclose(tuples2);}
}

void G4UserRunAction::atrisSave(){
  // IDs:
  //    * H1 ID 0 is the primary vertex energy
  //    * H2 ID 0 is the source direction
  //    * H2 ID 1 is the source position
  //    * H2 ID 2 is the ionization
  //    * H2 ID 3 is the dose_rate
  // STRATEGY:
  //    * acccess the histograms
  //    * remove the over/underbin bins
  //    * apply the per particle scaling using the h1 ID 0 entries 
  //    * apply for the ionization the per cm3 scaling using the sdid_volume data from DetectorConstruction
  //    * apply for the dose rate the per icru sphere scaling using the diid_icru_factor
  //    * saving information about statistical quality makes no sence since we almost always have more then 1e6 entries
  //    * include bin edges in the header
  //    * include units in the header
  //    * include short description in the header
  G4cout << "### AtRIS: now starting normalization of matrices" << G4endl;
  G4CsvAnalysisManager* analysis = G4CsvAnalysisManager::Instance();
  G4H1* ener = analysis->GetH1(0);
  G4H2* ioni = analysis->GetH2(2);
  // These are the entries from the energy histogram INCLUDING the overflow and underflow bins
  std::vector<unsigned int> entries = ener->bins_entries();
  std::vector<unsigned int> entries2 = ioni->bins_entries();
  std::vector<double>        values = ener->bins_sum_w();
  std::vector<double>        values2 = ioni->bins_sum_w();
  // Task: clean up, normalize with particle number, normalize with sensitive detector volume, convert to eV/cm^3
  // We define:
  //  nrofdetectors = number of sensitive detectors (including core, crust and loss volumes)
  //  nrbins        = number of 
  G4int rowsize = enbins.size() + 1;
  G4int colsize = nrofdetectors + 2;
  G4cout << "### AtRIS: h2_h2.ion should have " << rowsize*colsize << " entries. It has: " << entries2.size() << G4endl;
  // first h2_h2.ion entry corresponds to SDID=0, energy udnerflowbin, then we have SDID=0, first energy bin...
  // We remove all underflow and overflow bins
  // We reformat the data into a 2d array
  std::vector<double> norm_ioni;
  G4int end2 = values2.size();
  for (auto i=0; i<end2; i++){
	  bool cond1 = i > 2*rowsize;                // not planet underflow and not in core volume
	  bool cond2 = i < (end2 - 2*rowsize);        // not planet overflow and not in loss volume 
	  bool cond3 = (i%rowsize) != 0;             // not energy underflow
	  bool cond4 = (i%rowsize) != (rowsize - 1); // not energy overflow
	  bool cond = cond1 and cond2 and cond3 and cond4; // all together to keep
	  G4cout << "### AtRIS: Conditions for " << i << " " << cond1 << " " << cond2 << " " << cond3 << " " << cond4 << " keep point: "<< cond << " v=" << values2[i]/eV << " eV " << G4endl;
	  if (cond){
		  // We keep the point and normalize
		  G4int    nparticles = entries[i%rowsize];
		  G4double   sdvolume = sdid_volume[int (i/rowsize)-1]; //the -1 is because of the underflow bin
		  if (nparticles == 0){
			G4cout << "replacing number of particles 0 with 1 to aviod nan. The result will remain zero!" << G4endl;
			G4cout << "if you are getting this message, you are probbably doing something wront!" << G4endl;
			nparticles = 1;
		  }
		  G4double value = values2[i] /( sdvolume * G4double (nparticles));
		  G4cout << "npart = " << nparticles << " sdvol= " << sdvolume/cm3 <<"cm3" << " v=" << value/(eV/cm3) << "eV" << G4endl;
		  norm_ioni.push_back(value);
	  }
  }

  G4cout << "### AtRIS: cleanup and normalization are complete. Now savin custom output." << G4endl;
  // We:
  //   * open a file
  //   * generate a header with:
  //     - file description
  //     - unit specification
  //     - column specification
  //   * create a row (atmospheric or crust sheet)
  //   * write down 
  G4String header1 = "# AtRIS average ionization density\n";
  G4String header2 = "# SDID alt_min/cm I_0(eV/cm^3) I_1(eV/cm^3) ... I_N(eV/cm^3)\n";
  G4String header3 = "# I_j corresponds to the average ionization density in the j-th energy bin.\n";
  G4String header4 = "# alt_min is measured from the planetary center. We assume no gaps in thet planet\n";
  G4String header5 = "# Find the exact energy binning in the .bins file which you provided for the simulation.\n";
  G4String header6 = "# If you stack multiple files, you need to divide columns 3 and on with the number of files.";
  // Now saving
  G4String fname_ioni = G4UserRunAction::outputFileName + ".ioni"; 
  G4cout << "### AtRIS: using " << fname_ioni << " to store custom ionization response matrices" << G4endl;
  std::ofstream ionistream;
  ionistream.close();
  ionistream.open(fname_ioni, std::ofstream::out);
  ionistream << header1 << header2 << header3 << header4 << header5 << header6;
  for (unsigned int k=0; k<norm_ioni.size(); k++){
	// If at the beginning of a line, add a new line, the SDID and altitude. 
	if ((k%(rowsize-2))==0){
	       	ionistream << std::endl;
		G4int sdid = k/(rowsize-2)+1;
		ionistream << sdid << " " << sdid_rmin[sdid]/cm << " ";
	}
	ionistream << norm_ioni[k]/(eV/cm3) << " ";
  }
  ionistream << std::endl;
  ionistream.close();
  // Now we do the dose stuff

  //I'm completely aware that the following could have been implemented together with the above.
  //Since it does not cost a significant amount of cpu time, I've decided to implement it separaterly
  if (mode >= -1){
	  std::vector<double> norm_dose;
	  std::vector<double> norm_equi_dose;
	  std::vector<double> norm_dose2;
	  std::vector<double> norm_equi_dose2;
	  G4H2* dose = analysis->GetH2(3);
	  G4H2* equi = analysis->GetH2(4);
	  G4H2* dose2= analysis->GetH2(5);
	  G4H2* equi2= analysis->GetH2(6);
	  std::vector<unsigned int> entries3 = dose->bins_entries();
	  std::vector<double>        values3 = dose->bins_sum_w();
	  std::vector<unsigned int> entries4 = equi->bins_entries();
	  std::vector<double>        values4 = equi->bins_sum_w();
	  std::vector<unsigned int> entries5 = dose2->bins_entries();
	  std::vector<double>        values5 = dose2->bins_sum_w();
	  std::vector<unsigned int> entries6 = equi2->bins_entries();
	  std::vector<double>        values6 = equi2->bins_sum_w();
	  G4int end3 = values3.size();

	 G4cout << "using phantom mass of " << icrumass/kg << " kg" << G4endl;
	 G4cout << "using second phantom mass of " << icrumass2/kg << " kg" << G4endl;
	 G4double diidgeoffac = 2*pi*0.5*(1-std::cos(mmaxangle)*std::cos(mmaxangle));
	 G4cout << "The geometric factor of a detector (consult Sullivan paper)  using maxangle=" << mmaxangle/deg << "deg, is " << diidgeoffac <<" times the surface. For dose rate, this is not used for normalization." << G4endl;

	  // data formating and filtering
  	  for (auto i=0; i<end3; i++){
		  bool cond1 = i > 2*rowsize;                // not planet underflow and not in core volume
		  bool cond2 = i < (end3 - 3*rowsize);       // not planet overflow and not in loss volume 
		  bool cond3 = (i%rowsize) != 0;             // not energy underflow
		  bool cond4 = (i%rowsize) != (rowsize - 1); // not energy overflow
		  bool cond = cond1 and cond2 and cond3 and cond4; // all together to keep
		  G4cout << "### AtRIS: Conditions for " << i << " " << cond1 << " " << cond2 << " " << cond3 << " " << cond4 << " keep point: "<< cond << " en for dose1=" << values3[i]/eV << " eV, en for dose2=" << values5[i]/eV << "eV" << G4endl;
		  if (cond){
			  // We keep the point and normalize
			  G4double      value = values3[i];
			  G4double     evalue = values4[i];
			  G4double      value2= values5[i];
			  G4double     evalue2= values6[i];
			  G4int    nparticles = entries[i%rowsize];

			  G4double        fac = diid_icru_factor[int (i/rowsize) - 1];
			  G4double        fac2= diid_icru_factor2[int (i/rowsize) - 1];
			  if (nparticles == 0){
	  			  G4cout << "replacing number of particles 0 with 1 to aviod nan. The result will remain zero!" << G4endl;
	  			  G4cout << "if you are getting this message, you are probbably doing something wront!" << G4endl;
	  			  nparticles = 1;
			  }
			  value  =  value * fac /( G4double (nparticles) * icrumass);
			  evalue = evalue * fac /( G4double (nparticles) * icrumass);
			  value2 =  value2* fac2/( G4double (nparticles) * icrumass2);
			  evalue2= evalue2* fac2/( G4double (nparticles) * icrumass2);
			  G4cout << "npart = " << nparticles << ", fac= " << fac <<", mass=" << icrumass/kg << "kg, V=" << value/(eV/g) << " eV/g= "<< "eV, " << value/(joule/kg) << " J/kg, which is in gray: " << value/gray << ",  and also an equivalent " << evalue/gray << G4endl;
			  G4cout << "2ND npart = " << nparticles << ", fac= " << fac2<<", mass=" << icrumass2/kg << "kg, V=" << value2/(eV/g) << " eV/g= "<< "eV, " << value2/(joule/kg) << " J/kg, which is in gray: " << value2/gray << ",  and also an equivalent " << evalue2/gray << G4endl;
			  norm_dose.push_back(value);
			  norm_equi_dose.push_back(evalue);
			  norm_dose2.push_back(value2);
			  norm_equi_dose2.push_back(evalue2);
	  }
  }
	  
	  // File output
	  G4String bheader1 = "# AtRIS absorbed dose \n";
	  G4String eheader1 = "# AtRIS equivalent absorbed dose \n";
	  G4String bheader2 = "# DIID alt/cm D_1(Gy) ... D_N(Gy)Ni\n";
	  G4String eheader2 = "# DIID alt/cm H_1(Gy) ... H_N(Gy)Ni\n";
	  G4String bheader3 = "# D_j corresponds to the average absorbed energy dose rate for a particle in the j-th energy bin.\n";
	  G4String eheader3 = "# H_j corresponds to the average equivalent absorbed energy dose rate for a particle in the j-th energy bin.\n";
	  G4String bheader4 = "# alt_min is measured from the planetary center. We assume no gaps in thet planet\n";
	  G4String bheader5 = "# Find the exact energy binning in the .bins file which you provided for the simulation.\n";
	  G4String bheader6 = "# If you stack multiple files, you need to divide columns 3 and on with the number of files.\n";
	  G4String bheader71= "# using phantom file: " + phantomFileName;
	  G4String bheader72= "# using phantom file: " + phantomFileName2;
	  G4String fname_dose = G4UserRunAction::outputFileName + ".dose"; 
	  G4String fname_equi = G4UserRunAction::outputFileName + ".equi"; 
	  G4String fname_dose2= G4UserRunAction::outputFileName + ".dose2"; 
	  G4String fname_equi2= G4UserRunAction::outputFileName + ".equi2"; 
	  G4cout << "### AtRIS: using " << fname_dose << " and " << fname_equi << " to store custom response matrices" << G4endl;
	  G4cout << "### AtRIS: using " << fname_dose2<< " and " << fname_equi2<< " to store custom response matrices for the second phantom" << G4endl;
	  std::ofstream dosestream;
	  std::ofstream equistream;
	  std::ofstream dosestream2;
	  std::ofstream equistream2;
	  dosestream.close();
	  dosestream.open(fname_dose, std::ofstream::out);
	  dosestream << bheader1 << bheader2 << bheader3 << bheader4 << bheader5 << bheader6 <<bheader71 ;
	  equistream.close();
	  equistream.open(fname_equi, std::ofstream::out);
	  equistream << eheader1 << eheader2 << eheader3 << bheader4 << bheader5 << bheader6 <<bheader71 ;
	  dosestream2.close();
	  dosestream2.open(fname_dose2, std::ofstream::out);
	  dosestream2 << bheader1 << bheader2 << bheader3 << bheader4 << bheader5 << bheader6 <<bheader72 ;
	  equistream2.close();
	  equistream2.open(fname_equi2, std::ofstream::out);
	  equistream2 << eheader1 << eheader2 << eheader3 << bheader4 << bheader5 << bheader6 <<bheader72 ;
  	  for (unsigned int k=0; k<norm_dose.size(); k++){
	  	  // If at the beginning of a line, add a new line, the SDID and altitude. 
		  if ((k%(rowsize-2))==0){
	  		  dosestream << std::endl;
	  		  equistream << std::endl;
	  		  dosestream2<< std::endl;
	  		  equistream2<< std::endl;
	  		  G4int diid = k/(rowsize-2)+1;
	  		  dosestream << diid << " " << sdid_rmin[diid+1]/cm << " ";
	  		  equistream << diid << " " << sdid_rmin[diid+1]/cm << " ";
	  		  dosestream2<< diid << " " << sdid_rmin[diid+1]/cm << " ";
	  		  equistream2<< diid << " " << sdid_rmin[diid+1]/cm << " ";
	  	  }
	  	  dosestream << norm_dose[k]/gray << " ";
	  	  equistream << norm_equi_dose[k]/gray << " ";
	  	  dosestream2<< norm_dose2[k]/gray << " ";
	  	  equistream2<< norm_equi_dose2[k]/gray << " ";
	  }
	  dosestream << std::endl;
	  dosestream.close();
	  equistream << std::endl;
	  equistream.close();
	  dosestream2 << std::endl;
	  dosestream2.close();
	  equistream2 << std::endl;
	  equistream2.close();
	  }
}

void G4UserRunAction::SetPhantomName(const G4String& newFileName){
 	phantomFileName = newFileName; 
}

void G4UserRunAction::SetPhantomName2(const G4String& newFileName2){
 	phantomFileName2 = newFileName2; 
}


void G4UserRunAction::DefineCommands() {
	fMessenger  = new G4GenericMessenger(this,"/atris/","atris custom interface for configuration of processes happening after the simulation");
	// Set mass command
	G4GenericMessenger::Command& phantomMass = fMessenger->DeclarePropertyWithUnit("setPhantomMass","g",icrumass,"mass of the radiation phantom");
	phantomMass.SetParameterName("m", true);
	phantomMass.SetRange("m>0.");
	phantomMass.SetDefaultValue("14137.1669");
	// Set mass command for second phantom
	G4GenericMessenger::Command& phantomMass2 = fMessenger->DeclarePropertyWithUnit("setPhantom2Mass","mg",icrumass2,"mass of the second radiation phantom");
	phantomMass2.SetParameterName("m", true);
	phantomMass2.SetRange("m>0.");
	phantomMass2.SetDefaultValue("0.00052208");

	// set phantom file name
	auto& setPhantomNameCMD = fMessenger->DeclareMethod("setPhantomFilename", &G4UserRunAction::SetPhantomName, "set file name with phantom data");
	setPhantomNameCMD.SetParameterName("phantomFileName", false);
	setPhantomNameCMD.SetDefaultValue("icru_response.txt");
	// set phantom file name for second phantom
	auto& setPhantomNameCMD2 = fMessenger->DeclareMethod("setPhantom2Filename", &G4UserRunAction::SetPhantomName2, "set file name with second phantom data");
	setPhantomNameCMD2.SetParameterName("phantomFileName2", false);
	setPhantomNameCMD2.SetDefaultValue("cell_response.txt");
}

G4double G4UserRunAction::getPhantomMass(){
	return G4UserRunAction::icrumass;
}

G4double G4UserRunAction::getPhantomMass2(){
	return G4UserRunAction::icrumass2;
}
