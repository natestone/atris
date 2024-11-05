#include "UserActionInitialization.hh"
#include "G4UserRunAction.hh"
#include "PrimaryGeneratorAction.hh"

UserActionInitialization::UserActionInitialization() {
  ;
}
UserActionInitialization::~UserActionInitialization() {
  ;
}

void UserActionInitialization::Build() const {
  SetUserAction(new G4UserRunAction);
  SetUserAction(new PrimaryGeneratorAction);
}
  
void UserActionInitialization::BuildForMaster() const {
  SetUserAction(new G4UserRunAction);
}

