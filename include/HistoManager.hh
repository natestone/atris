#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include "G4Types.hh"
#include "G4String.hh"

#include "g4root.hh"

class HistoManager
{
  public:
   HistoManager();
  ~HistoManager();
  static G4String binedgesFile;
  static G4bool loadBins;

  private:
    void Book();
    G4String fFileName;
};


#endif

