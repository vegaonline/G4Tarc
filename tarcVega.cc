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

#include "G4Types.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4TARCDetectorConstruction.hh"
#include "G4TARCParallelWorld.hh"
#include "G4TARCPhysicsList.hh"
//#include "G4TARCPhysicsListMessenger.hh"
// #include "G4TARCPrimaryGeneratorAction.hh"
#include "G4TARCActionInitialization.hh"
//  #include "G4TARCRunAction.hh"
#include "G4TARCEventAction.hh"
#include "G4TARCSteppingAction.hh"
#include "G4TARCStackingAction.hh"

// I am adding all these and later we may be selective to pick effective ones  :: for biasing and scoring
#include "G4GeometrySampler.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ImportanceBiasing.hh"
#include "G4IStore.hh"
#include "G4VWeightWindowStore.hh"
#include "G4WeightWindowAlgorithm.hh"


// I am adding all these and later we may be selective to pick effective ones
#include "G4GeneralParticleSource.hh"
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
#include "G4Material.hh"
#include "G4TrajectoryDrawByParticleID.hh"

#include "FTFP_BERT.hh"
#include "FTFP_BERT_HP.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BIC_HP.hh"
#include "QGSP_BIC_AllHP.hh"
#include "QGSP_BIC.hh"
#include "QGS_BIC.hh"
#include "QGSP_INCLXX.hh"
#include "QGSP_INCLXX_HP.hh"
#include "QBBC.hh"
#include "LBE.hh"


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
  G4cerr << system("ls -l Data/GeomData/*.gdml");
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
  G4int nThreads = G4Threading::G4GetNumberOfCores();
#endif


  for (G4int i = 2; i < argc; i += 2) {
    G4String inArg = G4String(argv[i]);
    if (inArg == "-m") macro = argv[i + 1];
    if (inArg == "-u") session = argv[i + 1];
#ifdef G4MULTITHREADED
    if (inArg == "-t")
      nThreads = (G4String(argv[i+1]) == "NMAX")
               ? G4Threading::G4GetNumberOfCores()
               : G4UIcommand::ConvertToInt(argv[i + 1]);
#endif
    if (inArg != "-m" && inArg != "-u" && inArg != "-t") {
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
  scoreManager->SetVerboseLevel(1);

  G4String fileName = argv[1];

  //  construction of detector geometry setup
  G4String parallelWorldName = "tarcParallelWorld";
  auto geomConstruct = new G4TARCDetectorConstruction(fileName);
  runManager->SetUserInitialization(geomConstruct);    // RUNMANAGER for Geometry
  auto parallelWorld =  new G4TARCParallelWorld(parallelWorldName);
  geomConstruct->RegisterParallelWorld(parallelWorld);
  G4GeometrySampler pgsN(parallelWorld->GetWorldVolume(), "neutron");
  //G4GeometrySampler pgsP(parallelWorld->GetWorldVolume(), "proton");
  pgsN.SetParallel(true);
  //pgsP.SetParallel(true);

  // PhysicsList
  G4PhysListFactory           physFactory;
  G4VModularPhysicsList*      phys          = nullptr;
  //G4NeutronKiller*            neutKiller    = nullptr;
  //G4NeutronKillerMessenger*   neutKillMess  = nullptr;
  char*                       physnameInput = getenv("PHYSLIST");
  G4String                    physicsName   = "";
  if (physnameInput)          physicsName   = G4String(physnameInput);
  if (!physicsName.size())    physicsName   = "QGSP_BERT_HP";

  if (physicsName.size() && physFactory.IsReferencePhysList(physicsName)){
      phys = physFactory.GetReferencePhysList(physicsName);
      phys->RegisterPhysics(new G4RadioactiveDecayPhysics);
    //  neutKiller = new G4NeutronKiller();
    //  neutKillMess = new G4NeutronKillerMessenger(neutKiller);
  }
  if (!phys) { phys = new G4TARCPhysicsList(); }
  phys->RegisterPhysics(new G4ImportanceBiasing(&pgsN, parallelWorldName));
  phys->RegisterPhysics(new G4ParallelWorldPhysics(parallelWorldName));

  runManager->SetUserInitialization(phys);      // RUNMANAGER for Physics List

  // Action Initialization
  runManager->SetUserInitialization(new G4TARCActionInitialization());

  // Initialize Runmanager
  runManager->Initialize();
  parallelWorld->CreateImportanceStore();
  // Create new drawByParticleID model
  auto model = new G4TrajectoryDrawByParticleID;
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
    UImanager->ApplyCommand("/control/execute vis.mac");
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv, session);
      ui->SessionStart();
      delete ui;
#endif
  }

  G4GeometryManager::GetInstance()->OpenGeometry();
  // pgsN.ClearSampling();

  //termination of job
#ifdef G4VIS_USE
  delete visManager;
  G4cout << "Vis manager deleted" << G4endl;
#endif
  //delete runManager;
  //G4cout << "run manager deleted" << G4endl;

  return 0;
}
