/**********************************************************************
 * @file   : tarcVega.cc
 * @author : Abhijit Bhattacharyya
 * @brief  : This ia a code for TARC application using G4
 *           The GDML geometry file has been taken from Alexander Howard
 *           This code is in concurrence with GV version
 *********************************************************************/
// $Id: tarcVega.cc Nov-17-2017 11:03 vega $

#include <vector>
#include <iostream>
#include <string>

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4TARCDetectorConstruction.hh"
#include "G4TARCParallelWorld.hh"
#include "G4TARCPhysicsList.hh"
//#include "G4TARCPhysicsListMessenger.hh"
#include "G4TARCPrimaryGeneratorAction.hh"
#include "G4TARCActionInitialization.hh"
#include "G4TARCRunAction.hh"
#include "G4TARCEventAction.hh"
#include "G4TARCStackingAction.hh"

// I am adding all these and later we may be selective to pick effective ones
#include "G4ParallelWorldPhysics.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "G4NeutronKiller.hh"
#include "G4NeutronKillerMessenger.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4RadioactiveDecayMode.hh"
#include "G4RadioactiveDecaymessenger.hh"
#include "G4ScoringManager.hh"
#include "G4NuclearLevelData.hh"
#include "G4DecayPhysics.hh"
#include "G4ParticleTable.hh"
#include "G4VTrajectoryModel.hh"
#include "G4TrajectoryDrawByParticleID.hh"

#include "FTFP_BERT.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4GDMLParser.hh"
#include "Randomize.hh"

void PrintUsage() {
  G4cerr << G4endl;
  G4cerr << "Usage: tarcVega <input_gdml_file:mandatory> [-m macro] [-u UIsession] [-t nThreads] "
         << G4endl;
  G4cerr << "\t \t Here -t option is available for multithreaded mode only." << G4endl;
  G4cerr << "The geometry that can be used from GeomData directory is :" << G4endl;
  G4cerr << system("ls -l GeomData/");
  G4cerr << G4endl;
}

int main(int argc, char** argv) {
  G4String macro;
  G4String session;

  if (argc < 2 || argc > 8) {
    PrintUsage();
    return -1;
  }


  // Choose the Random Engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

#ifdef G4MULTITHREADED
  G4int nThreads = G4Threading::G4GetNumberOfCores();   //   / 2;
#endif
  for (G4int i = 2; i < argc; i += 2) {
    if      (G4String(argv[i]) == "-m") macro = argv[i+1];
    else if (G4String(argv[i]) == "-u") session = argv[i+1];
#ifdef G4MULTITHREADED
    else if (G4String(argv[i]) == "-t") {
      if (G4String(argv[i+1]) == "NMAX") {
        nThreads = G4Threading::G4GetNumberOfCores();
      } else {
        nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
      }
    }
#endif
    else {
      PrintUsage();
      return -1;
    }
  }

  // construction of default run manager
  #ifdef G4MULTITHREADED
  	auto runManager = new G4MTRunManager;
    if ( nThreads > 0 ) {
      runManager->SetNumberOfThreads(nThreads);
    }
    G4cout << "Code is started with " << runManager->GetNumberOfThreads() << " threads. " << G4endl;
  #else
  	auto runManager = new G4RunManager;
  #endif

  auto scoreManager = G4ScoringManager::GetScoringManager();
  scoreManager->SetVerboseLevel(0);

  G4String fileName = argv[1];

  //  construction of detector geometry setup
  G4String parallelWorldName = "tarcParallelWorld";                // trying to use it
  auto geomConstruct = new G4TARCDetectorConstruction(fileName);
  auto parallelWorld =  new G4TARCParallelWorld(parallelWorldName);
  geomConstruct->RegisterParallelWorld(parallelWorld);
  runManager->SetUserInitialization(geomConstruct);    // RUNMANAGER for Geometry

  // PhysicsList
  G4PhysListFactory          physFactory;   //= new G4PhysListFactory();
  //G4TARCPhysicsList*          phys        = new G4TARCPhysicsList();
  G4VModularPhysicsList*      phys          = 0;
  //G4TARCPhysicsListMessenger* physMess      = 0;
  G4NeutronKiller*            neutKiller    = 0;
  G4NeutronKillerMessenger*   neutKillMess  = 0;
  char*                       physnameInput = getenv("PHYSLIST");
  G4String                    physicsName   = "";
  if (physnameInput)          physicsName   = G4String(physnameInput);

  if (physicsName.size() && physFactory.IsReferencePhysList(physicsName)){
      phys = physFactory.GetReferencePhysList(physicsName);
      phys->RegisterPhysics(new G4RadioactiveDecayPhysics);
      //physMess = new G4TARCPhysicsListMessenger(phys);
      neutKiller = new G4NeutronKiller();
      neutKillMess = new G4NeutronKillerMessenger(neutKiller);
  }
  if (!phys) { phys = new G4TARCPhysicsList(); }
  phys->RegisterPhysics(new G4ParallelWorldPhysics(parallelWorldName, true));
  //physMess = new G4TARCPhysicsListMessenger(phys);
  runManager->SetUserInitialization(phys);      // RUNMANAGER for Physics List

  // Action Initialization
  #ifdef G4MULTITHREADED
    runManager->SetUserInitialization(new G4TARCActionInitialization());
  #else
      runManager->SetUserAction( new G4TARCPrimaryGeneratorAction() );

      runManager->SetUserAction( new G4TARCRunAction() );

      runManager->SetUserAction( new G4TARCEventAction() );

      runManager->SetUserAction( new G4TARCStackingAction() );

  #endif

  // Create new drawByParticleID model
  G4TrajectoryDrawByParticleID* model = new G4TrajectoryDrawByParticleID;
  // Configure model
  model->SetDefault("white");
  model->Set("gamma", "green");
  model->Set("e+", "cyan");
  model->Set("e-", "blue");
  model->Set("proton", "red");
  model->Set("neutron", "yellow"); //G4Color(0.7, 0.1, 0.4));
  //Register model with visualization manager


  // Visualization
  auto visManager = new G4VisExecutive;
  visManager->Initialize();
  visManager->RegisterModel(model);

  // get pointer to user interface manager
  auto UImanager =  G4UImanager::GetUIpointer();

  if (macro.size()){
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  } else {
    #ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv, session);
      ui->SessionStart();
      delete ui;
    #endif
  }

  //termination of job
  #ifdef G4VIS_USE
    delete visManager;
  #endif
  delete runManager;
  //delete physMess;

  return 0;
}
