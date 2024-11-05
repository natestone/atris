#include "HistoManager.hh"
#include "G4CsvAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "DetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserRunAction.hh"
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>

extern G4int year;
extern G4int nrofdetectors;
extern G4int eventGoal;

std::vector<G4double> enbins;

G4String HistoManager::binedgesFile;
G4bool HistoManager::loadBins = false;

HistoManager::HistoManager()
  : fFileName("exgps")
{
  Book();
}

HistoManager::~HistoManager()
{
  G4cout << "HistoManager deleted." << G4endl;
  delete G4CsvAnalysisManager::Instance();
}

void HistoManager::Book()
{
  G4CsvAnalysisManager* analysis = G4CsvAnalysisManager::Instance();
  
  analysis->SetFileName(fFileName);
  analysis->SetVerboseLevel(0);
  analysis->SetActivation(true);     //enable inactivation of histos, nTuples

  // Create all histograms as inactivated 
  G4double newbinedge;
  // Here we load the data file containing binedges in MeV
  std::ifstream is(HistoManager::binedgesFile);
  G4int j = 0; 
  while (is >> newbinedge){
	  enbins.push_back(newbinedge);
	  j++;
  }
  G4cout << "### AtRIS: read " << j << " binedges " <<G4endl;
  std::cout << "### AtRIS: read " << j << " binedges " <<G4endl;

  G4int    id   = analysis->CreateH1("h1.1","kinetic energy", enbins, "MeV");
  analysis->SetH1Activation(id, true);
    
  // This is the direction histogram.
  G4int    dir_nbinx = 120;
  G4int    dir_nbiny = 60 ;
  G4double dir_vminx = 0.;
  G4double dir_vmaxx = twopi;
  G4double dir_vminy = 0.;
  G4double dir_vmaxy = 1.;
  id = analysis->CreateH2("h2.dir","source direction: phi,cos(theta)",
		  dir_nbinx,dir_vminx,dir_vmaxx,dir_nbiny,dir_vminy,dir_vmaxy);
  analysis->SetH2Activation(id,true);
	
  // This is the position histogram.
  G4int    pos_nbinx = 120;
  G4int    pos_nbiny = 60 ;
  G4double pos_vminx = 0.;
  G4double pos_vmaxx = twopi;
  G4double pos_vminy = 0.;
  G4double pos_vmaxy = 1.;
  id = analysis->CreateH2("h2.pos","source position: phi,cos(theta)",
		  pos_nbinx,pos_vminx,pos_vmaxx,pos_nbiny,pos_vminy,pos_vmaxy);
  analysis->SetH2Activation(id,true);

  // This is the new ionization energy histogram
  std::vector<G4double> detbins;
  G4cout << "Generating detid bins. Starting with an empty vector... adding bin edges: " << G4endl;
  for(int i=0;i<=nrofdetectors;i++){
	detbins.push_back(i-0.5);
	G4String detname = std::to_string(i);
  }
  id = analysis->CreateH2("h2.ion","ionization: input energy,SID",enbins,detbins,"MeV","none","none","none");
  analysis->SetH2Activation(id,true);
  // For the DIID, we need less detbins, but we will keep this here and remove it latter when writing out the custom output
  id = analysis->CreateH2("h2.dos","absorbed dose rate in a spherical water phatnom: input energy,SID",enbins,detbins,"MeV","none","none","none");
  analysis->SetH2Activation(id,true);
  // same as above, but only for the equivalent dose rate  
  id = analysis->CreateH2("h2.equi","equivalent absorbed dose rate in a spherical water phatnom: input energy,SID",enbins,detbins,"MeV","none","none","none");
  analysis->SetH2Activation(id,true);
  // second phantom: For the DIID, we need less detbins, but we will keep this here and remove it latter when writing out the custom output
  id = analysis->CreateH2("h2.dos2","absorbed dose rate in a second spherical phatnom: input energy,SID",enbins,detbins,"MeV","none","none","none");
  analysis->SetH2Activation(id,true);
  // second phantom same as above, but only for the equivalent dose rate  
  id = analysis->CreateH2("h2.equi2","equivalent absorbed dose rate in a second spherical phatnom: input energy,SID",enbins,detbins,"MeV","none","none","none");
  analysis->SetH2Activation(id,true);
  //analysis->SetH1Plotting(id,true);
  //analysis->SetH2Plotting(id,true);
}

