#ifndef SensitiveDetector_hh
#define SensitiveDetector_hh

#include "G4VSensitiveDetector.hh"
#include "G4UImessenger.hh"

#include <vector>
#include <fstream>

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"
#include "G4UserRunAction.hh"
#include "DetectorConstruction.hh"
#include <functional>


class SDHit;


class SDHit : public G4VHit
{
public:
	SDHit();
	virtual ~SDHit();

	inline void* operator new(size_t);
	inline void  operator delete(void*);

	virtual void Print();

	// variables
	//G4ThreeVector pos, dir, pv_pos, pv_dir, prepos;
	//G4double ekin, pv_ekin, len, ioni, nioni;
	G4double ioni, nioni;
	//G4String part, det, vol, proc;
	G4String det, part;
	//G4int evid, trid, parid, runid, detid;
	G4int detid, stepid;
	

	static G4Allocator<SDHit>* SDHitAllocator;
};

typedef G4THitsCollection<SDHit> SDHitsCollection;

inline void* SDHit::operator new(size_t)
{
	if(!SDHitAllocator)
		SDHitAllocator = new G4Allocator<SDHit>;
	return (void *) SDHitAllocator->MallocSingle();
}

inline void SDHit::operator delete(void *hit)
{
	SDHitAllocator->FreeSingle((SDHit*) hit);
}


class G4Step;
class G4HCofThisEvent;


class SensitiveDetector : public G4VSensitiveDetector, public G4UImessenger
{
public:
	SensitiveDetector(const G4String& name, G4bool zero_E=false);
	virtual ~SensitiveDetector();

	// methods from base class
	virtual void   Initialize(G4HCofThisEvent* hitCollection);
	virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
	//virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

private:

	SDHitsCollection* fHitsCollection;
	//G4VProcess* fpCreatorProcess;

protected:
	static std::string formatString;
	G4bool zeroE;

	static std::ofstream outstream;
};

#endif
