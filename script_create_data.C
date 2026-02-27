#include "SimulationTools.h"

#include <iostream>
#include <fstream>
#include <time.h>

TString particleName[5] = { "proton", "deuteron", "triton", "3He", "alpha" };

//void deblur_and_update( int fileN=1 )
int main( int argc, char** argv )
{
  clock_t time_start, time_end;
  time_start = clock();

  // -------------------------------------------------------
  // Save parameters from the input text file
  int nEvents = 0;
  int proj_Z = 0; int proj_A = 0;
  int targ_Z = 0; int targ_A = 0;
  double ke_beam = 0;
  bool CallModel = 1;
  bool IsCentral = 0;
  bool fUseTotalQ= 0;
  bool fOptimal  = 0;
  bool fOptAzi   = 0;

  TString dummy;
  TString inputFileName = argv[1];
  ifstream iFile( inputFileName );
  while( !iFile.eof() )
  {
    iFile >> dummy >> nEvents;
    if( iFile.eof() ) break;
    iFile >> dummy >> proj_Z;
    if( iFile.eof() ) break;
    iFile >> dummy >> proj_A;
    if( iFile.eof() ) break;
    iFile >> dummy >> targ_Z;
    if( iFile.eof() ) break;
    iFile >> dummy >> targ_A;
    if( iFile.eof() ) break;
    iFile >> dummy >> ke_beam;
    if( iFile.eof() ) break;
    iFile >> dummy >> CallModel;
    if( iFile.eof() ) break;
    iFile >> dummy >> IsCentral;
    if( iFile.eof() ) break;
    iFile >> dummy >> fUseTotalQ;
    if( iFile.eof() ) break;
    iFile >> dummy >> fOptimal;
    if( iFile.eof() ) break;
    iFile >> dummy >> fOptAzi;
    if( iFile.eof() ) break;

    cout << "nEvents : " << nEvents << endl;
    cout << "proj_Z : " << proj_Z << " and proj_A : " << proj_A << endl;
    cout << "targ_Z : " << targ_Z << " and targ_A : " << targ_A << endl;
    cout << "Beam Energy : " << ke_beam << endl;
    if( IsCentral ) cout << "Central collision events" << endl;
    else cout << "Mid-central collision events" << endl;
    if( fOptimal ) cout << "Switch to Optimal weight" << endl;
    else cout << "Standard weight" << endl;
  }
  TString collisionEvent = IsCentral ? "central" : "midcentral";
  // Save parameters from the input text file
  // -------------------------------------------------------

  SimulationTools doSimulation;
  doSimulation.SetOptimizedAzimuth( fOptAzi );
  doSimulation.SetCollisionProperties( proj_Z, proj_A, targ_Z, targ_A, ke_beam );
  if( CallModel ) 
  {
    double avg_charge = 50;
    LocalEquilibriumModel* thisModel = new LocalEquilibriumModel();
    thisModel->PrepareParameters( ke_beam, avg_charge );
    thisModel->CalcAvgMulti();
    doSimulation.CallParticleFrom( thisModel );


    TString effFileName = Form( "/home/jhpark/work/deblurring_spirit/merging/efficiency/eff_Sn%d_%s.root", proj_A, collisionEvent.Data() );
    auto efficiencyFile = new TFile( effFileName, "READ" ); // Like Exp
    doSimulation.EnableRandomReactionPlane();
    doSimulation.EnableEfficiency();
    doSimulation.SetEfficiencyFile( efficiencyFile );

    nEvents *= 5;
  }
  else 
  {
    TFile* dataFile = new TFile( Form("expdata/tree_Sn%d.root", proj_A), "READ" );
    TTree* dataTree = (TTree*) dataFile->Get( "tree" );
    doSimulation.CallParticleFrom( dataTree );
    nEvents = dataTree->GetEntries();
  }


  TString fileName = Form( "output_dist_Sn%d.root", proj_A );
  doSimulation.SetOutputfileName( fileName );
  int initialized = doSimulation.Init();
  if( initialized==-1 )
  {
    cout << "Initialized failure" << endl;
    return -1;
  }

  doSimulation.SetNumberOfEvents( nEvents );

  // Create data first
  doSimulation.CreateTheFirstData();
  // RL deblurring w/ central limit blurring function
  doSimulation.RLdeblurring();

  // If optimal weights you are going to use,
  if( fOptimal )
  {
    // Setup optimal weight flag
    doSimulation.SetUseOptimal( fOptimal );
    // Calculate optimal weight
    doSimulation.SetupOptimal();
    // Set flag not to update restored distribution in saving data
    doSimulation.DoNotUpdateRestoredDistribution();
    // Then produce the data again
    doSimulation.CreateTheFirstData();
  }

  doSimulation.SaveOutputFiles();

  time_end = clock();
  double res = (double) (time_end - time_start) / CLOCKS_PER_SEC;
  cout << res << endl;

  return 0;

}
