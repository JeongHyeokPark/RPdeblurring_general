#include "SimulationTools.h"

#include <iostream>
#include <fstream>
#include <time.h>

TString particleName[5] = { "proton", "deuteron", "triton", "3He", "alpha" };

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



  // -------------------------------------------------------
  // Call ROOT file from the input parameter
  TString fileROOTinputName = argv[2];
  auto file = new TFile( fileROOTinputName, "READ" );
  if( file==nullptr )
  {
    cout << "Please make sure that you are feeding input ROOT file containing distribution (TH3)" << endl;
    return -1;
  }
  // And call TH3D from the file
  std::vector< TH3D* > hist_pt_rap_phi_rest;
  for( int ipid=0; ipid<5; ipid++ ) hist_pt_rap_phi_rest.push_back( (TH3D*) file->Get( Form("hist_pt_rap_phi_rest_%s", particleName[ipid].Data()) ) );
  std::vector< TH3D* > hist_pt_rap_phi_meas;
  for( int ipid=0; ipid<5; ipid++ ) hist_pt_rap_phi_meas.push_back( (TH3D*) file->Get( Form("hist_pt_rap_phi_meas_%s", particleName[ipid].Data()) ) );
  // Call ROOT file from the input parameter
  // -------------------------------------------------------


  // -------------------------------------------------------
  // Setting random seed from the input number
  TString input_value_string  = argv[3];
  int input_value = input_value_string.Atoi();

  int fileN = input_value;

#if defined(__linux__)
  gRandom->SetSeed( time(0) * long(fileN) * gethostid() ); // To set random seed depending on time * computer id
#elif defined(__APPLE__)
  gRandom->SetSeed( time(0) * long(fileN) ); // For MAC os
#endif
  // Setting random seed from the input number
  // -------------------------------------------------------




  SimulationTools doSimulation;
  doSimulation.SetCollisionProperties( proj_Z, proj_A, targ_Z, targ_A, ke_beam );
  doSimulation.CallParticleFrom( hist_pt_rap_phi_rest );
  doSimulation.SetNumberOfEvents( nEvents );
  doSimulation.EnableRandomReactionPlane();

  // Efficiency file needs to be assigned here to construct TM
  TString effFileName = Form( "/home/jhpark/work/deblurring_spirit/merging/efficiency/eff_Sn%d_%s.root", proj_A, collisionEvent.Data() );
  //TString effFileName = Form( "./efficiency/eff_Sn%d_%s.root", proj_A, collisionEvent.Data() );
  auto efficiencyFile = new TFile( effFileName, "READ" ); // Like Exp
  cout << effFileName << endl;
  doSimulation.EnableEfficiency();
  doSimulation.SetEfficiencyFile( efficiencyFile );

  doSimulation.SetUseTotalQ( fUseTotalQ );
  doSimulation.SetUseOptimal( fOptimal );
  if( fOptimal ) doSimulation.SetMeasuredHist( hist_pt_rap_phi_meas );

  TString fileName = Form( "output_%d.root", fileN );
  doSimulation.SetOutputfileName( fileName );
  doSimulation.Init();


  // Set flag not to update restored distribution in saving data
  doSimulation.DoNotUpdateRestoredDistribution();
  doSimulation.SampleAndMakeTM();

  doSimulation.SaveOutputFiles();

  time_end = clock();
  double res = (double) (time_end - time_start) / CLOCKS_PER_SEC;
  cout << res << endl;

  return 0;

}
