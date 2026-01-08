#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TEfficiency.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <unordered_map>

#include "LocalEquilibriumModel.h"

using namespace std;

class SimulationTools {

  public:
    // Main functions
    void CreateTheFirstData();
    void SampleAndMakeTM();
    void RLdeblurring();
    // Main functions

    int Init();
    void PrepareMass();
    void SetSeed( long val ) { fSeed = val; };

    void SetNumberOfEvents( int val ) { fNumberOfEvents = val; };
    void SetCollisionProperties( int proj_Z, int proj_A, int targ_Z, int targ_A, double collisionE );
    void EnableRandomReactionPlane() { fRandomPhi = true; };

    void EnableEfficiency() { fEfficiencyMode = true; };

    void SetEfficiency( TEfficiency* eff ) { fEfficiency = eff; };
    void SetEfficiencyFile( TFile* file ) { fEfficiencyFile = file; };
    double GetEfficiency( int ipid, double pt, double y0, double phiLab );

    void CallParticleFrom( TTree* tree ) { fDataTree = tree; }; // Sample particles from distribution
    void CallParticleFrom( std::vector<TH3D*> &histRestored ) { fHistogramRestored = histRestored; }; // Sample particles from distribution
    void CallParticleFrom( LocalEquilibriumModel* model ) { fModel = model; }; // Sample particles from flow model
    void CallParticle(); // Used in "SampleAndMakeTM" method
    void DoNotUpdateRestoredDistribution() { fUpdateRestored = false; };

    void SetOutputfileName( TString fileName ) { fOutputName = fileName; };
    void SaveOutputFiles();

    void SetUseTotalQ( bool val ) { fUseTotalQ = val; };
    void ConstructTMfromQ();

    void SetMeasuredHist( std::vector<TH3D*> &histMeasured ) { fHistogramMeasured = histMeasured; };
    void UpdateInternalDistribution() { fUpdateInternal = true; };
    void SetUseOptimal( bool val ) { bUseOptimal = val; };
    void SetupOptimal();

    void PrepareRLDeblurring( TFile* thisFile ) { fRLfile = thisFile; };
    void RandomizeRegularization();

    void SetTotalQMatrix( TH2D* hist ) { hist_totalq = hist; };
    void NormalizeEstimatedDistribution();
    void EstimaedError( int nSample );


  private:
    static const int fMultiMax = 200;
    static const int nParticleType = 5;
    TString particleName[6] = { "proton", "deuteron", "triton", "3He", "alpha", "neutron" };

    double fDelta = 0.17; // Cuts for weighting factor
    double fEffCoarse = 0.7;

    long fSeed = 1;

    int nEvents = 0;

    static const int nBinPtMax = 25;
    int nBinPt = -1;
    double binWidthPt[nBinPtMax] = {0};
    double ptMin = -1.;
    double ptMax = -1.;

    int nBinRap = -1;
    double binWidthRap = -1;
    double yyMin = -1.;
    double yyMax = -1.;

    int nBinPhi = -1; // Total number of bins of histogram to be saved
    int nBinPhy = -1 ; // Number of bins in half way
    double binWidthPhy = -1;
    double phiMin = -99;
    double phiMax = -99;

    int nBinQphy = -1;
    double binWidthQphy = -1;
    double QphiMin = -1;
    double QphiMax = -1;
    int nBinQphi = -1;


    bool   fUpdateRestored = 1;
    int    fNumberOfEvents = -1;
    bool   fRandomPhi = false;
    bool   fEfficiencyMode = false;
    TEfficiency* fEfficiency = nullptr;
    TFile* fEfficiencyFile = nullptr;
    std::vector< TEfficiency* > fEfficiency_vector;
    std::vector< TH3D* > fEfficiency_hist;

    bool fUseTotalQ = false;

    int fProjZ = -1;
    int fProjA = -1;
    int fTargZ = -1;
    int fTargA = -1;
    double fEn = -1;
    double p_beam = -1;
    double y_cm = -99;

    int tree_nTracks = -1;
    int tree_pid[fMultiMax]     = {0};
    double tree_px[fMultiMax]   = {0.};
    double tree_py[fMultiMax]   = {0.};
    double tree_pz[fMultiMax]   = {0.};
    double tree_mass[fMultiMax] = {0.};
    TTree* fDataTree = nullptr;


    LocalEquilibriumModel* fModel = nullptr;
    std::vector<TH3D*> fHistogramRestored;
    std::vector<TH3D*> fHistogramMeasured;
    std::vector< std::vector< int > > fCheckedBin;

    std::vector<TLorentzVector> *fVector;
    std::vector<int> *fZ;
    std::vector<int> *fA;
    std::vector<double> *fMass;

    std::vector<std::vector<std::vector<TH2D*>>> hist_phi_diff; // in terms of pID, pT, rapidity
    std::vector<TH3D*> hist_pt_rap_phi_real; // in terms of pID
    std::vector<TH3D*> hist_pt_rap_phi_meas; // in terms of pID
    std::vector<TH3D*> hist_pt_rap_phi_rest; // in terms of pID
    std::vector<std::vector<std::vector<TH2D*>>> hist_prev_phi_diff; // in terms of pID, pT, rapidity

    TFile* fRLfile = nullptr;
    std::vector<TH3D*> hist_pt_rap_phi_blur; // For RL Deblurring
    std::vector<TH3D*> hist_pt_rap_phi_temp; // For RL Deblurring


    double particle_mass[6]; // p, d, t, 3He, a, n
    int particle_a[6]       = { 1, 2, 3, 3, 4, 1 }; // p, d, t, 3He, a, n
    int particle_z[6]       = { 1, 1, 1, 2, 2, 0 }; // p, d, t, 3He, a, n

    TH2D* hist_totalq;
    TH1D* hist_central;
    TH2D* hist_rp_phi_diff;

    TString fOutputName = "output.root";

    bool fUpdateInternal = false;

    double dnCut = 0.15;
    std::vector< std::vector<double> > optWeight;
    bool bUseOptimal = false;

    double fLambda = 0.007;
    double fLambda_org = -1;
    int modulo( int a, int b ) { return ((a % b) + b) % b; }

};

void SimulationTools::PrepareMass()
{
  double amu  = .931502; // amu in GeV
  double ampu = 1.00727647;
  double amnu = 1.00866501;
  double amdu = 2.014102;
  double amtu = 3.016049;
  double amhu = 3.016029;
  double amau = 4.002603;

  double mass_p = amu*ampu;
  double mass_n = amu*amnu;
  double mass_d = amu*amdu;
  double mass_t = amu*amtu;
  double mass_h = amu*amhu;
  double mass_a = amu*amau;


  particle_mass[0] = mass_p;
  particle_mass[1] = mass_d;
  particle_mass[2] = mass_t;
  particle_mass[3] = mass_h;
  particle_mass[4] = mass_a;
  particle_mass[5] = mass_n;
  // ---------------- contents in flowsym_mod.for -------------------------
}

void SimulationTools::SetCollisionProperties( int proj_Z, int proj_A, int targ_Z, int targ_A, double collisionE )
{
  fProjZ = proj_Z;
  fProjA = proj_A;
  fTargZ = targ_Z;
  fTargA = targ_A;
  fEn = collisionE;
}

void SimulationTools::RandomizeRegularization()
{
  double y = gRandom->Uniform( -1., 1. );
  double z = TMath::ACos( y );
  double sinz = TMath::Sqrt( 1. - y*y );

  while( sinz < gRandom->Uniform() )
  {
    y = gRandom->Uniform( -1., 1. );
    z = TMath::ACos( y );
    sinz = TMath::Sqrt( 1. - y*y );
  }
  double sin2ran = z/TMath::Pi();

  fLambda = fLambda_org * (0.5 + sin2ran);
}

int SimulationTools::Init()
{
  nEvents = 0;

  if( fDataTree==nullptr && fHistogramRestored.size() == 0 && fModel == nullptr ) cout << "Please set model or TH3D class as an input to emit the particles" << endl;
  if( fDataTree==nullptr && fHistogramRestored.size() == 0 && fModel == nullptr ) return -1;

  PrepareMass();

  fVector = new std::vector<TLorentzVector>;
  fZ      = new std::vector<int>;
  fA      = new std::vector<int>;
  fMass   = new std::vector<double>;

  // ---------------------- Calculate boost vector ----------------------
  double mass_proj = .931502 * 131.9178239/132.;
  //double mass_proj = .931502 * 107.9119653/108.;
  double e_beam = mass_proj + fEn;
  p_beam = std::sqrt( (e_beam+mass_proj) * (e_beam-mass_proj) );
  double y_beam = std::log( (e_beam+p_beam) / mass_proj );
  y_cm = 0.5 * y_beam;
  cout << "y_cm: " << y_cm << endl;
  // ---------------------- Calculate boost vector ----------------------



  // ---------------------------------- Prepare histograms ---------------------------------- 
  if( fDataTree || fModel )
  {
    nBinPt = 25;
    for( int ipid=0; ipid<nParticleType; ipid++ ) binWidthPt[ipid] = .10;
    binWidthPt[0] = .10; binWidthPt[1] = .15; binWidthPt[2] = .23; binWidthPt[3] = .29; binWidthPt[4] = .23;

    // Rapidity in CM / reduced rapidity
    binWidthRap = 0.20;
    yyMin = -2. - binWidthRap/2.;
    yyMax =  2. + binWidthRap/2.;
    nBinRap = (yyMax - yyMin) / binWidthRap + 0.1;

    nBinPhy = 12;
    binWidthPhy = TMath::Pi() / nBinPhy;
    phiMin = -nBinPhy*binWidthPhy - binWidthPhy/2.;
    phiMax =  nBinPhy*binWidthPhy + binWidthPhy/2.;
    nBinPhi = (phiMax - phiMin) / binWidthPhy + 0.1;
  }
  if( fHistogramRestored.size()!=0 )
  {
    nBinPt = fHistogramRestored[0]->GetNbinsX();
    for( int ipid=0; ipid<nParticleType; ipid++ ) binWidthPt[ipid] = fHistogramRestored[ipid]->GetXaxis()->GetBinWidth( 1 );
    ptMin = fHistogramRestored[0]->GetXaxis()->GetXmin();
    ptMax = fHistogramRestored[0]->GetXaxis()->GetXmax();

    nBinRap = fHistogramRestored[0]->GetNbinsY();
    yyMin = fHistogramRestored[0]->GetYaxis()->GetXmin();
    yyMax = fHistogramRestored[0]->GetYaxis()->GetXmax();
    binWidthRap = fHistogramRestored[0]->GetYaxis()->GetBinWidth( 1 );

    nBinPhi = fHistogramRestored[0]->GetNbinsZ();
    phiMin = fHistogramRestored[0]->GetZaxis()->GetXmin();
    phiMax = fHistogramRestored[0]->GetZaxis()->GetXmax();
    binWidthPhy = fHistogramRestored[0]->GetZaxis()->GetBinWidth( 1 );
    nBinPhy = (nBinPhi-1)/2;
  }

  nBinQphy = nBinPhy;
  binWidthQphy = TMath::Pi() / nBinQphy;
  QphiMin = -nBinQphy*binWidthQphy - binWidthQphy/2.;
  QphiMax =  nBinQphy*binWidthQphy + binWidthQphy/2.;
  nBinQphi = (QphiMax - QphiMin) / binWidthQphy + 0.1;

  hist_central = new TH1D( "hist_central",  "", nBinQphi, QphiMin, QphiMax );
  hist_totalq  = new TH2D( "hist_totalq",  "", nBinQphi, QphiMin, QphiMax, nBinQphi, QphiMin, QphiMax );


  if( bUseOptimal ) SetupOptimal();


  // Resize first
  hist_pt_rap_phi_real.resize( nParticleType ); // in terms of pID
  hist_pt_rap_phi_meas.resize( nParticleType ); // in terms of pID
  hist_pt_rap_phi_rest.resize( nParticleType ); // in terms of pID

  hist_pt_rap_phi_blur.resize( nParticleType ); // For Deblurring
  hist_pt_rap_phi_temp.resize( nParticleType ); // For Deblurring


  hist_phi_diff.resize( nParticleType ); // in terms of pID, pT, rapidity
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_phi_diff.at(ipid).resize( nBinPt ); // Slice pT bins
  for( int ipid=0; ipid<nParticleType; ipid++ ) for( int ipt=0; ipt<nBinPt; ipt++ ) hist_phi_diff.at(ipid).at(ipt).resize( nBinRap ); // Slice rapidity bins

  hist_prev_phi_diff.resize( nParticleType ); // in terms of pID, pT, rapidity
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_prev_phi_diff.at(ipid).resize( nBinPt ); // Slice pT bins
  for( int ipid=0; ipid<nParticleType; ipid++ ) for( int ipt=0; ipt<nBinPt; ipt++ ) hist_prev_phi_diff.at(ipid).at(ipt).resize( nBinRap ); // Slice rapidity bins
  // Resize first


  for( int ipid=0; ipid<nParticleType; ipid++ ) 
  {
    ptMin = binWidthPt[ipid]*0.;
    ptMax = binWidthPt[ipid]*nBinPt;
    hist_pt_rap_phi_real[ipid] = new TH3D( Form("hist_pt_rap_phi_real_%s", particleName[ipid].Data()), "real dist; pt (GeV/c); Rapidity (y_{0}); #phi",      nBinPt, ptMin, ptMax, nBinRap, yyMin, yyMax, nBinPhi, phiMin, phiMax );
    hist_pt_rap_phi_meas[ipid] = new TH3D( Form("hist_pt_rap_phi_meas_%s", particleName[ipid].Data()), "measured dist; pt (GeV/c); Rapidity (y_{0}); #phi'", nBinPt, ptMin, ptMax, nBinRap, yyMin, yyMax, nBinPhi, phiMin, phiMax );
    hist_pt_rap_phi_rest[ipid] = new TH3D( Form("hist_pt_rap_phi_rest_%s", particleName[ipid].Data()), "restored dist; pt (GeV/c); Rapidity (y_{0}); #phi'", nBinPt, ptMin, ptMax, nBinRap, yyMin, yyMax, nBinPhi, phiMin, phiMax );

    // For Deblurring
    hist_pt_rap_phi_blur[ipid] = new TH3D( Form("hist_pt_rap_phi_blur_%s", particleName[ipid].Data()), "", nBinPt, ptMin, ptMax, nBinRap, yyMin, yyMax, nBinPhi, phiMin, phiMax );
    hist_pt_rap_phi_temp[ipid] = new TH3D( Form("hist_pt_rap_phi_temp_%s", particleName[ipid].Data()), "", nBinPt, ptMin, ptMax, nBinRap, yyMin, yyMax, nBinPhi, phiMin, phiMax );

    for( int ipt=0; ipt<nBinPt; ipt++ )
      for( int iy=0; iy<nBinRap; iy++ ) 
        hist_phi_diff[ipid][ipt][iy] = new TH2D( Form("hist_phi_diff_%s_p%d_y%d",particleName[ipid].Data(),ipt,iy) , "estimated ptcles phi; #phi; #phi'", nBinPhi, phiMin, phiMax, nBinPhi, phiMin, phiMax );
  }



  int nBinRPphi = 20;
  double binWidthRPphi = TMath::Pi() / nBinRPphi;
  double RPphiMin1 = -nBinRPphi*binWidthRPphi - binWidthRPphi/2.;
  double RPphiMax1 =  nBinRPphi*binWidthRPphi + binWidthRPphi/2.;
  int nBinRPphi1 = (RPphiMax1 - RPphiMin1) / binWidthRPphi + 0.1;
  hist_rp_phi_diff = new TH2D( "hist_rp_phi_diff", "", nBinRPphi1, RPphiMin1, RPphiMax1, nBinRPphi1, RPphiMin1, RPphiMax1 );
  // ---------------------------------- Prepare histograms ---------------------------------- 


  // ----------------------------------------------------------------------------------------------------
  if( fDataTree ) 
  {
    fDataTree->SetBranchAddress( "nTracks", &tree_nTracks );
    fDataTree->SetBranchAddress( "pid",     tree_pid     );
    fDataTree->SetBranchAddress( "px",      tree_px      );
    fDataTree->SetBranchAddress( "py",      tree_py      );
    fDataTree->SetBranchAddress( "pz",      tree_pz      );
    fDataTree->SetBranchAddress( "mass",    tree_mass    );
  }
  // ----------------------------------------------------------------------------------------------------

  if( fEfficiencyMode )
  {
    fEfficiency_vector.resize( 0 );
    fEfficiency_hist.resize( 0 );
    if( fEfficiencyFile != nullptr )
      for( int i=0; i<nParticleType; i++ )
      {
        auto theEfficiency = (TEfficiency*) fEfficiencyFile->Get( Form("e3_%s_HighRes", particleName[i].Data()) );
        fEfficiency_vector.push_back( theEfficiency );

        auto thisHist = (TH3D*) theEfficiency->GetPassedHistogram();
        fEfficiency_hist.push_back( thisHist );
      }
    else if( fEfficiency != nullptr )
      for( int i=0; i<nParticleType; i++ ) fEfficiency_vector.push_back( fEfficiency );
    else cout << "Please indicate the efficiency file first." << endl;
  }


  if( fRLfile )
  {
    cout << "Deblurring with: " << nBinPt << " " << nBinRap << " " << nBinPhi << " binning" << endl;
    for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_real[ipid] = (TH3D*) fRLfile->Get( Form("hist_pt_rap_phi_corr_%s", particleName[ipid].Data()) );
    //for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_real[ipid] = (TH3D*) fRLfile->Get( Form("hist_pt_rap_phi_real_%s", particleName[ipid].Data()) );
    for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_meas[ipid] = (TH3D*) fRLfile->Get( Form("hist_pt_rap_phi_meas_%s", particleName[ipid].Data()) );
    for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_rest[ipid] = (TH3D*) fRLfile->Get( Form("hist_pt_rap_phi_rest_%s", particleName[ipid].Data()) );
    for( int ipid=0; ipid<nParticleType; ipid++ )
      for( int ipt=0; ipt<nBinPt; ipt++ )
        for( int iy=0; iy<nBinRap; iy++ ) hist_phi_diff[ipid][ipt][iy] = (TH2D*) fRLfile->Get( Form("hist_phi_diff_%s_p%d_y%d",particleName[ipid].Data(),ipt,iy) );

    // We don't need to initialize below if we are only going to do RL deblurring
    return 0;
  }

  // For RL deblurring
  fLambda_org = fLambda;;


  // ----------------------------------------------------------------------------------------------------
  if( fHistogramRestored.size() != 0 )
  {
    cout << "Scanning non zero bins of group including " << fHistogramRestored[0]->GetName() << endl;

    // Scanning non-zero bins
    fCheckedBin.resize( nParticleType );
    for( int ipid=0; ipid<nParticleType; ipid++ )
    {
      for( int i=0; i<nBinPt; i++ )
      {
        for( int j=0; j<nBinRap; j++ )
        {
          for( int k=0; k<nBinPhi-1; k++ ) // The both end bins are the same bin
          {
            int thisBin = fHistogramRestored[ipid]->GetBin( i+1, j+1, k+1 );
            double binContent = fHistogramRestored[ipid]->GetBinContent( i+1, j+1, k+1 );
            if( 0. < binContent ) hist_pt_rap_phi_rest[ipid]->SetBinContent( i+1, j+1, k+1, binContent );
            if( 0. < binContent ) fCheckedBin[ipid].push_back( thisBin );

            binContent = fHistogramMeasured[ipid]->GetBinContent( i+1, j+1, k+1 );
            if( 0. < binContent ) hist_pt_rap_phi_meas[ipid]->SetBinContent( i+1, j+1, k+1, binContent );
          }
        }
      }
      cout << particleName[ipid] << " histogram has " << fCheckedBin[ipid].size() << " bins of non zero value." << endl;
    }
  }
  // ----------------------------------------------------------------------------------------------------



  return 0;
}


void SimulationTools::CreateTheFirstData()
{
  for( int ipid=0; ipid<nParticleType; ipid++ )
  {
    hist_pt_rap_phi_real[ipid]->Reset();
    hist_pt_rap_phi_meas[ipid]->Reset();
    if( fUpdateRestored ) hist_pt_rap_phi_rest[ipid]->Reset();
  }

  double qaqa_sum = 0.;
  double qka_sum  = 0.;
  double wga_sum  = 0.;
  double wbar_sum = 0.;

  for( int iEvent=0; iEvent<fNumberOfEvents; iEvent++ ) // loop over events
  {
    if( iEvent%50000==0 )  gRandom->SetSeed( fSeed * time(0) );

    int particleZ[fMultiMax]           = {0};
    int particleA[fMultiMax]           = {0};
    bool particleDet[fMultiMax]        = {0};
    double particlePx[fMultiMax]       = {0};
    double particlePy[fMultiMax]       = {0};
    double particlePz[fMultiMax]       = {0};
    double particleE[fMultiMax]        = {0};
    double particleMass[fMultiMax]     = {0};
    double particlePhi[fMultiMax]      = {0};

    double particleY[fMultiMax]  = {0}; double particleID[fMultiMax] = {0};
    double particleQx[fMultiMax] = {0}; double particleQy[fMultiMax] = {0};

    double phi_rp = 0.;
    if( fRandomPhi ) phi_rp = gRandom->Uniform( -TMath::Pi(), TMath::Pi() );

    if( fDataTree ) fDataTree->GetEntry( iEvent );
    CallParticle();
    int nTracks = fVector->size();

    // ---------------------------------------------------------------
    // Central limit assumption
    double qx_event    = 0.;
    double qy_event    = 0.;
    double qx2         = 0.;
    double qy2         = 0.;
    double qx2_event   = 0.;
    double qy2_event   = 0.;
    double wbar        = 0.;
    double wbar_event  = 0.;
    double wbar2_event = 0.;
    // Central limit assumption
    // ---------------------------------------------------------------


    if( iEvent%100000==0 ) cout << iEvent << " calling " << nTracks << " particles" << endl;
    for( int iTrack=0; iTrack<nTracks; iTrack++ )
    {
      particleDet[iTrack] = 0; // Initialize that this track is not detected

      double mass        = fMass->at( iTrack );
      double mass_square = mass*mass;

      int particle_type = -1;
      int IZI   = fZ->at( iTrack );
      int IAI   = fA->at( iTrack );
      if( IZI==1 && IAI==1 ) particle_type = 0;
      if( IZI==1 && IAI==2 ) particle_type = 1;
      if( IZI==1 && IAI==3 ) particle_type = 2;
      if( IZI==2 && IAI==3 ) particle_type = 3;
      if( IZI==2 && IAI==4 ) particle_type = 4;
      if( IZI==0 && IAI==1 ) particle_type = 5; // In case of neutron analysis

      particleZ[iTrack] = IZI;
      particleA[iTrack] = IAI;
      particleMass[iTrack] = mass;
      particleID[iTrack] = particle_type;


      auto L = fVector->at( iTrack );
      double normalized_rapidity = L.Rapidity();

      double pt = L.Pt();
      double mt = TMath::Sqrt( pt*pt + mass*mass );

      // L.Rapidity should be y0, yCM/yBeam
      double rapidity = normalized_rapidity*y_cm;
      double ylab = rapidity + y_cm;
      double particle_phi = L.Phi();

      double elab   = mt * TMath::CosH( ylab );
      double plab_z = mt * TMath::SinH( ylab );
      double plab_x = pt * TMath::Cos( particle_phi );
      double plab_y = pt * TMath::Sin( particle_phi );

      TVector3 mom( plab_x, plab_y, plab_z );

      // This should be consistent as in producing Exp Data.
      particleY[iTrack] = normalized_rapidity;

      int ipid = particleID[iTrack];
      particlePhi[iTrack] = particle_phi;

      // ----------------------- Rotation with reaction plane angle ----------------------- 
      // No rotation if reading exp data
      mom.RotateZ( phi_rp );
      particlePx[iTrack] = mom.X();
      particlePy[iTrack] = mom.Y();
      particlePz[iTrack] = mom.Z();
      particleE[iTrack]  = elab;
      // ----------------------- Rotation with reaction plane angle ----------------------- 


      // ----------------------- Efficiency Filter ----------------------- 
      double eff = 1.;
      if( fModel && fEfficiencyMode )
      {
        // ----------------------- Efficiency( Pt, Rapidity, Phi ) ----------------------- 
        double phiInDeg = mom.Phi() * TMath::RadToDeg();
        int thisBin = fEfficiency_vector[ipid]->FindFixBin( pt*1e+3, normalized_rapidity, phiInDeg );
        //eff = fEfficiency_vector[ipid]->GetEfficiency( thisBin );
        eff = GetEfficiency( ipid, pt*1e+3, normalized_rapidity, phiInDeg );
        // ----------------------- Efficiency( Pt, Rapidity, Phi ) ----------------------- 

        if( gRandom->Uniform() < eff )
          particleDet[iTrack] = 1; // particle detected
        else 
          particleDet[iTrack] = 0; // particle NOT detected
      }
      else particleDet[iTrack] = 1; // Don't consider efficiency filter
      // ----------------------- Efficiency Filter ----------------------- 


      // ---------------------------------------------------------------------------------
      // -------------------------- Considering over/under-flow --------------------------
      int iBinPt = hist_pt_rap_phi_real[ipid]->GetXaxis()->FindBin( pt );
      int iBinRap = hist_pt_rap_phi_real[ipid]->GetYaxis()->FindBin( particleY[iTrack] );
      //if( !(0<iBinPt && iBinPt<nBinPt && 0<iBinRap && iBinRap<nBinRap) ) continue;
      if( nBinPt<iBinPt || nBinRap<iBinRap ) continue;
      // -------------------------- Considering over/under-flow --------------------------
      // ---------------------------------------------------------------------------------

      if( fModel ) hist_pt_rap_phi_real[ipid]->Fill( pt, particleY[iTrack], particle_phi ); 

      // ----------------------- Calculation for Q vector construction ----------------------- 
      if( particleDet[iTrack]==1 )
      {
        double weight_factor = 0.;
        if( normalized_rapidity > +fDelta ) weight_factor = +1.; // forward
        else if( normalized_rapidity < -fDelta ) weight_factor = -1.; // backward
        else weight_factor = 0.;

        int iBinRapWeight = iBinRap - 1; // To compensate the index in array
        // Don't use optimal weight when initially produce the data
        if( bUseOptimal )
          weight_factor = optWeight[ipid][iBinRapWeight];

        double qx_i = weight_factor * particlePx[iTrack];
        double qy_i = weight_factor * particlePy[iTrack];


        particleQx[iTrack] = qx_i;
        particleQy[iTrack] = qy_i;

        qx_event = qx_event + qx_i;
        qy_event = qy_event + qy_i;

        // --------------------------------------------------------------------------------------------------
        // Central Limit function parametrization
        qx2 = qx_i*qx_i;
        qy2 = qy_i*qy_i;
        qx2_event += qx2;
        qy2_event += qy2;

        wbar = std::abs(weight_factor) * particle_a[ipid];
        wbar_event  += wbar;
        wbar2_event += wbar*wbar;
        // --------------------------------------------------------------------------------------------------

      }
      // ----------------------- Calculation for Q vector construction ----------------------- 


    } // loop end ptcle type & filling of event completed


    // ----------------------------------------------------------------------------------------
    // Central limit assumption
    double qaqa     = 0.;
    double qk_r     = 0.;
    double qka      = 0.;

    double qk = qx2_event + qy2_event;
    double dqq = qx_event*qx_event + qy_event*qy_event - qk;
    double wbwb = wbar_event*wbar_event - wbar2_event;
    if( 1.1 < wbwb )
    {
      qaqa     = dqq  / wbwb;
      qk_r     = qk   - wbar_event*qaqa;
      qka      = qk_r / wbar_event;
      qaqa_sum = qaqa_sum + wbar_event*qaqa;
      qka_sum  = qka_sum  + wbar_event*qka;
      wga_sum  = wga_sum  + wbar_event;
      wbar_sum = wbar_sum + wbar_event*wbar_event;
    }
    // Central limit assumption
    // ----------------------------------------------------------------------------------------


    for( int iTrack=0; iTrack<nTracks; iTrack++ )
    {
      if( particleDet[iTrack]==0 ) continue;

      double phiq = 0.;
      double QX1 = qx_event - particleQx[iTrack];
      double QY1 = qy_event - particleQy[iTrack];
      if( QX1!=0. || QY1!=0. ) phiq = TMath::ATan2( QY1, QX1 );
      else phiq = gRandom->Uniform( -TMath::Pi(), TMath::Pi() );

      double phi_diff = phiq - phi_rp;
      if( phi_diff < -TMath::Pi() ) phi_diff = phi_diff + TMath::TwoPi();
      if( phi_diff >  TMath::Pi() ) phi_diff = phi_diff - TMath::TwoPi();

      TVector3 mom( particlePx[iTrack], particlePy[iTrack], particlePz[iTrack] );
      double pt = mom.Pt();

      double normalized_rapidity = particleY[iTrack];

      int ipid = particleID[iTrack];
      double particlePhi_lab = mom.Phi();
      double particlePhi_rel = particlePhi_lab - phiq;
      if( particlePhi_rel < -TMath::Pi() ) particlePhi_rel = particlePhi_rel + TMath::TwoPi();
      if( particlePhi_rel >  TMath::Pi() ) particlePhi_rel = particlePhi_rel - TMath::TwoPi();

      hist_pt_rap_phi_meas[ipid]->Fill( pt, normalized_rapidity, particlePhi_rel );

      if( hist_rp_phi_diff->GetXaxis()->FindBin( phi_rp ) == hist_rp_phi_diff->GetNbinsX() ) phi_rp = phi_rp - TMath::TwoPi();
      if( hist_rp_phi_diff->GetXaxis()->FindBin( phiq ) == hist_rp_phi_diff->GetNbinsX() ) phiq = phiq - TMath::TwoPi();
      hist_rp_phi_diff->Fill( phi_rp, phiq );
    }


  } // loop end events



  // ------------------------------------------------------------------------------------------------------
  // Normalize the data
  for( int ipid=0; ipid<nParticleType; ipid++ )
    hist_pt_rap_phi_meas[ipid]->Scale( 1./fNumberOfEvents / binWidthPhy / binWidthRap / binWidthPt[ipid] );

  for( int ipid=0; ipid<nParticleType; ipid++ )
    for( int ipt=0; ipt<nBinPt; ipt++ )
      for( int irap=0; irap<nBinRap; irap++ )
        for( int iphi=0; iphi<nBinPhi; iphi++ )
        {
          double val_meas = hist_pt_rap_phi_meas[ipid]->GetBinContent( ipt+1, irap+1, iphi+1 );
          double err_meas = hist_pt_rap_phi_meas[ipid]->GetBinError( ipt+1, irap+1, iphi+1 );
          double cen_meas = hist_pt_rap_phi_meas[ipid]->GetXaxis()->GetBinCenter( ipt+1 );

          hist_pt_rap_phi_meas[ipid]->SetBinContent( ipt+1, irap+1, iphi+1, val_meas/cen_meas );
          hist_pt_rap_phi_meas[ipid]->SetBinError( ipt+1, irap+1, iphi+1, err_meas/cen_meas );
        }
  // Normalize the data
  // ------------------------------------------------------------------------------------------------------


  // Compensating the end of the bins
  for( int ipid=0; ipid<nParticleType; ipid++ )
  {
    for( int irap=0; irap<nBinRap; irap++ )
    {
      for( int ix=0; ix<nBinPt; ix++ )
      {
        double avgReal = ( hist_pt_rap_phi_real[ipid]->GetBinContent( ix+1, irap+1, 1 ) + hist_pt_rap_phi_real[ipid]->GetBinContent( ix+1, irap+1, nBinPhi ) );
        hist_pt_rap_phi_real[ipid]->SetBinContent( ix+1, irap+1, 1,       avgReal );
        hist_pt_rap_phi_real[ipid]->SetBinContent( ix+1, irap+1, nBinPhi, avgReal );

        double avgMeas = ( hist_pt_rap_phi_meas[ipid]->GetBinContent( ix+1, irap+1, 1 ) + hist_pt_rap_phi_meas[ipid]->GetBinContent( ix+1, irap+1, nBinPhi ) );
        hist_pt_rap_phi_meas[ipid]->SetBinContent( ix+1, irap+1, 1,       avgMeas );
        hist_pt_rap_phi_meas[ipid]->SetBinContent( ix+1, irap+1, nBinPhi, avgMeas );

        double errMeas = ( hist_pt_rap_phi_meas[ipid]->GetBinError( ix+1, irap+1, 1 ) + hist_pt_rap_phi_meas[ipid]->GetBinError( ix+1, irap+1, nBinPhi ) );
        hist_pt_rap_phi_meas[ipid]->SetBinError( ix+1, irap+1, 1,       errMeas );
        hist_pt_rap_phi_meas[ipid]->SetBinError( ix+1, irap+1, nBinPhi, errMeas );
      }
    }
  }
  // Compensating the end of the bins

  if( fUpdateRestored )
  {
    for( int ipid=0; ipid<nParticleType; ipid++ )
    {
      hist_pt_rap_phi_rest[ipid] = (TH3D*) hist_pt_rap_phi_meas[ipid]->Clone( hist_pt_rap_phi_rest[ipid]->GetName() );

      for( int ipt=0; ipt<nBinPt; ipt++ )
        for( int irap=0; irap<nBinRap; irap++ )
          for( int iphi=0; iphi<nBinPhi; iphi++ )
            hist_pt_rap_phi_rest[ipid]->SetBinError( ipt+1, irap+1, iphi+1, 0 );

      hist_pt_rap_phi_rest[ipid]->Scale( 1./fEffCoarse );
    }


  // ----------------------------------------------------------------------------------------
  // Central limit assumption
    wga_sum        = std::max( wga_sum, 1e0 );
    double wbar    = wbar_sum / wga_sum;
    double qaqa    = qaqa_sum / wga_sum;
    double qka     = qka_sum  / wga_sum;
    double qkato   = wbar * qka;
    double qato    = wbar * std::sqrt( std::max( qaqa, 0e0 ) );
    double facq    = qato / std::sqrt( qkato );
    double aph     = std::exp( -facq*facq ) / TMath::TwoPi();
    cout.precision( 12 );
    cout << setw(12) << wga_sum << " " << setw(12) << wbar << " " << setw(12) << qaqa << " " << setw(12) << qka << " " << setw(12) << qkato << " " << setw(12) << qato << endl;
    cout << "facq/aph " << setw(12) << facq << " " << setw(12) << aph << endl;

    double summ = 0.;
    for( int iphi=0; iphi<nBinQphi; iphi++ )
    {
      double phi = hist_central->GetBinCenter( iphi+1 );

      double bph = facq * TMath::Cos( phi );
      double val = aph * ( 1. + TMath::Sqrt(TMath::Pi()) * bph * TMath::Exp( bph*bph ) * ( 1. + TMath::Erf(bph) ) );
      val *= hist_central->GetBinWidth( 1 );

      hist_central->SetBinContent( iphi+1, val );
      if( iphi<nBinQphi-1 ) summ += val; // To avoid the last bin
    }
    cout << "summ: " << summ << endl;
    hist_central->Scale( 1./summ );
  // Central Limit assumption
  // ----------------------------------------------------------------------------------------


  // ----------------------------------------------------------------------------------------
  // Register blurring function from central limit
  for( int ipid=0; ipid<nParticleType; ipid++ )
    for( int ipt=0; ipt<nBinPt; ipt++ )
      for( int iy=0; iy<nBinRap; iy++ )
        for( int ireal=0; ireal<nBinPhi; ireal++ )
          for( int imeas=0; imeas<nBinPhi; imeas++ )
          {
            double phi_real = hist_central->GetBinCenter( ireal+1 );
            double phi_meas = hist_central->GetBinCenter( imeas+1 );
            //double phi_real = ( ireal+nBinPhy )*binWidthPhy + gRandom->Uniform( -binWidthPhy/2., binWidthPhy/2. );
            //double phi_meas = ( imeas+nBinPhy )*binWidthPhy + gRandom->Uniform( -binWidthPhy/2., binWidthPhy/2. );
            double phi_diff = phi_meas - phi_real;
            if( phi_diff < -TMath::Pi() ) phi_diff = phi_diff + TMath::TwoPi();
            if( phi_diff >  TMath::Pi() ) phi_diff = phi_diff - TMath::TwoPi();
            double val = hist_central->GetBinContent( hist_central->FindBin(phi_diff) );
            val *= fEffCoarse; // Coarse efficiency

            hist_phi_diff[ipid][ipt][iy]->SetBinContent( ireal+1, imeas+1, val );
          }
  }
  // Register blurring function from central limit
  // ----------------------------------------------------------------------------------------
}

void SimulationTools::SampleAndMakeTM()
{
  for( int ipid=0; ipid<nParticleType; ipid++ )
  {
    hist_pt_rap_phi_real[ipid]->Reset();
    //hist_pt_rap_phi_meas[ipid]->Reset();
    if( fUpdateRestored ) hist_pt_rap_phi_rest[ipid]->Reset();
  }

  for( int iEvent=0; iEvent<fNumberOfEvents; iEvent++ ) // loop over events
  {
    if( iEvent%100000==0 ) cout << iEvent << endl;
    if( iEvent%50000==0 )  gRandom->SetSeed( fSeed * time(0) );

    int particleZ[fMultiMax]           = {0};
    int particleA[fMultiMax]           = {0};
    bool particleDet[fMultiMax]        = {0};
    double particlePx[fMultiMax]       = {0};
    double particlePy[fMultiMax]       = {0};
    double particlePz[fMultiMax]       = {0};
    double particleE[fMultiMax]        = {0};
    double particleMass[fMultiMax]     = {0};
    double particlePhi[fMultiMax]      = {0};

    double particleY[fMultiMax]  = {0}; double particleID[fMultiMax] = {0};
    double particleQx[fMultiMax] = {0}; double particleQy[fMultiMax] = {0};

    double phi_rp = 0.;
    if( fRandomPhi ) phi_rp = gRandom->Uniform( -TMath::Pi(), TMath::Pi() );

    CallParticle();
    int nTracks = fVector->size();

    double qx_event    = 0.;
    double qy_event    = 0.;

    if( iEvent%100000==0 ) cout << iEvent << " calling " << nTracks << " particles" << endl;
    for( int iTrack=0; iTrack<nTracks; iTrack++ )
    {
      particleDet[iTrack] = 0; // Initialize that this track is not detected

      double mass        = fMass->at( iTrack );
      double mass_square = mass*mass;

      int particle_type = -1;
      int IZI   = fZ->at( iTrack );
      int IAI   = fA->at( iTrack );
      if( IZI==1 && IAI==1 ) particle_type = 0;
      if( IZI==1 && IAI==2 ) particle_type = 1;
      if( IZI==1 && IAI==3 ) particle_type = 2;
      if( IZI==2 && IAI==3 ) particle_type = 3;
      if( IZI==2 && IAI==4 ) particle_type = 4;
      if( IZI==0 && IAI==1 ) particle_type = 5; // In case of neutron analysis

      particleZ[iTrack] = IZI;
      particleA[iTrack] = IAI;
      particleMass[iTrack] = mass;
      particleID[iTrack] = particle_type;


      auto L = fVector->at( iTrack );
      double pt = L.Pt();
      double mt = TMath::Sqrt( pt*pt + mass*mass );

      // L.Rapidity should be y0, yCM/yBeam
      double rapidity = L.Rapidity()*y_cm;
      double ylab = rapidity + y_cm;
      double particle_phi = L.Phi();

      double elab   = mt * TMath::CosH( ylab );
      double plab_z = mt * TMath::SinH( ylab );
      double plab_x = pt * TMath::Cos( particle_phi );
      double plab_y = pt * TMath::Sin( particle_phi );


      TVector3 mom( plab_x, plab_y, plab_z );

      // This should be consistent as in producing Exp Data.
      double normalized_rapidity = ylab/y_cm - 1.;
      particleY[iTrack] = normalized_rapidity;

      int ipid = particleID[iTrack];
      particlePhi[iTrack] = particle_phi;


      // ----------------------- Rotation with reaction plane angle ----------------------- 
      mom.RotateZ( phi_rp );
      particlePx[iTrack] = mom.X();
      particlePy[iTrack] = mom.Y();
      particlePz[iTrack] = mom.Z();
      particleE[iTrack]  = elab;
      // ----------------------- Rotation with reaction plane angle ----------------------- 


      // ----------------------- Efficiency Filter ----------------------- 
      double eff = 1.;
      if( fEfficiencyMode )
      {
        // ----------------------- Efficiency( Pt, Rapidity, Phi ) ----------------------- 
        double phiInDeg = mom.Phi() * TMath::RadToDeg();
        int thisBin = fEfficiency_vector[ipid]->FindFixBin( pt*1e+3, normalized_rapidity, phiInDeg );
        //eff = fEfficiency_vector[ipid]->GetEfficiency( thisBin );
        eff = GetEfficiency( ipid, pt*1e+3, normalized_rapidity, phiInDeg );
        // ----------------------- Efficiency( Pt, Rapidity, Phi ) ----------------------- 

        if( gRandom->Uniform() < eff )
          particleDet[iTrack] = 1; // particle detected
        else 
          particleDet[iTrack] = 0; // particle NOT detected
      }
      else particleDet[iTrack] = 1; // Don't consider efficiency filter
      // ----------------------- Efficiency Filter ----------------------- 


      // ---------------------------------------------------------------------------------
      // -------------------------- Considering over/under-flow --------------------------
      int iBinPt = hist_pt_rap_phi_real[ipid]->GetXaxis()->FindBin( pt );
      int iBinRap = hist_pt_rap_phi_real[ipid]->GetYaxis()->FindBin( particleY[iTrack] );
      //if( !(0<iBinPt && iBinPt<nBinPt && 0<iBinRap && iBinRap<nBinRap) ) continue;
      if( nBinPt<iBinPt || nBinRap<iBinRap ) continue;
      // -------------------------- Considering over/under-flow --------------------------
      // ---------------------------------------------------------------------------------

      hist_pt_rap_phi_real[ipid]->Fill( pt, particleY[iTrack], particle_phi ); 


      // ----------------------- Calculation for Q vector construction ----------------------- 
      if( particleDet[iTrack]==1 )
      {
        double theta = mom.Theta();
        double phi = mom.Phi();

        double weight_factor = 0.;
        if( normalized_rapidity > +fDelta ) weight_factor = +1.; // forward
        else if( normalized_rapidity < -fDelta ) weight_factor = -1.; // backward
        else weight_factor = 0.;

        int iBinRapWeight = iBinRap - 1; // To compensate the index in array
        if( bUseOptimal )
          weight_factor = optWeight[ipid][iBinRapWeight];

        double qx_i = weight_factor * particlePx[iTrack];
        double qy_i = weight_factor * particlePy[iTrack];

        particleQx[iTrack] = qx_i;
        particleQy[iTrack] = qy_i;

        qx_event = qx_event + qx_i;
        qy_event = qy_event + qy_i;
      }
      // ----------------------- Calculation for Q vector construction ----------------------- 


    } // loop end ptcle type & filling of event completed


    nEvents++;


    if( fUseTotalQ )
    {
      // We accumulate the reaction plane direction with total Q vector, not omitting particles to be analyzed
      double phiq = 0.;
      if( qx_event!=0. && qy_event!=0. ) phiq = TMath::ATan2( qy_event, qx_event );
      else phiq = gRandom->Uniform( -TMath::Pi(), TMath::Pi() );

      if( phiq < -TMath::Pi() ) phiq = phiq + TMath::TwoPi();
      if( phiq >  TMath::Pi() ) phiq = phiq - TMath::TwoPi();


      double phiqq = phiq;
      double phi_save_rp = phi_rp;
      if( hist_totalq->GetXaxis()->FindBin( phiqq ) == nBinPhi ) phiqq = phiqq - TMath::TwoPi();
      if( hist_totalq->GetXaxis()->FindBin( phi_save_rp ) == nBinPhi ) phi_save_rp = phi_save_rp - TMath::TwoPi();

      // phi_save_rp = REAL RP angle
      // phiq        = ESTIMATED RP angle
      hist_totalq->Fill( phi_save_rp, phiqq );
      hist_rp_phi_diff->Fill( phi_rp, phiq );
    }
    else
    {
      // This section calculates directly the relations of the reaction plane direction with omitting particles to be analyzed
      for( int iTrack=0; iTrack<nTracks; iTrack++ )
      {
        if( particleDet[iTrack]==0 ) continue;

        double phiq = 0.;
        double QX1 = qx_event - particleQx[iTrack];
        double QY1 = qy_event - particleQy[iTrack];
        if( QX1!=0. || QY1!=0. ) phiq = TMath::ATan2( QY1, QX1 );
        else phiq = gRandom->Uniform( -TMath::Pi(), TMath::Pi() );

        double phi_diff = phiq - phi_rp;
        if( phi_diff < -TMath::Pi() ) phi_diff = phi_diff + TMath::TwoPi();
        if( phi_diff >  TMath::Pi() ) phi_diff = phi_diff - TMath::TwoPi();

        TVector3 mom( particlePx[iTrack], particlePy[iTrack], particlePz[iTrack] );
        double pt = mom.Pt();

        double normalized_rapidity = particleY[iTrack];

        int ipid = particleID[iTrack];
        double particlePhi_lab = mom.Phi();
        double particlePhi_rel = particlePhi_lab - phiq;
        if( particlePhi_rel < -TMath::Pi() ) particlePhi_rel = particlePhi_rel + TMath::TwoPi();
        if( particlePhi_rel >  TMath::Pi() ) particlePhi_rel = particlePhi_rel - TMath::TwoPi();

        hist_pt_rap_phi_meas[ipid]->Fill( pt, normalized_rapidity, particlePhi_rel );

        double particlePhi_rel_real = particlePhi[iTrack];
        if( hist_pt_rap_phi_real[ipid]->GetZaxis()->FindBin( particlePhi_rel_real ) == nBinPhi ) particlePhi_rel_real = particlePhi_rel_real - TMath::TwoPi();
        if( hist_pt_rap_phi_meas[ipid]->GetZaxis()->FindBin( particlePhi_rel ) == nBinPhi ) particlePhi_rel = particlePhi_rel - TMath::TwoPi();

        int ptBin = hist_pt_rap_phi_meas[ipid]->GetXaxis()->FindBin( pt ) - 1;
        int rapBin = hist_pt_rap_phi_meas[ipid]->GetYaxis()->FindBin( normalized_rapidity ) - 1;
        if( !fUseTotalQ ) if( 0<=ptBin && ptBin<nBinPt && 0<=rapBin && rapBin<nBinRap ) hist_phi_diff[ipid][ptBin][rapBin]->Fill( particlePhi_rel_real, particlePhi_rel );


        if( hist_rp_phi_diff->GetXaxis()->FindBin( phi_rp ) == hist_rp_phi_diff->GetNbinsX() ) phi_rp = phi_rp - TMath::TwoPi();
        if( hist_rp_phi_diff->GetXaxis()->FindBin( phiq ) == hist_rp_phi_diff->GetNbinsX() ) phiq = phiq - TMath::TwoPi();

        hist_rp_phi_diff->Fill( phi_rp, phiq );
      }
    }


  } // loop end events



  // Compensating the end of the bins
  for( int ipid=0; ipid<nParticleType; ipid++ )
  {
    for( int irap=0; irap<nBinRap; irap++ )
    {
      for( int ix=0; ix<nBinPt; ix++ )
      {
        double avgReal = ( hist_pt_rap_phi_real[ipid]->GetBinContent( ix+1, irap+1, 1 ) + hist_pt_rap_phi_real[ipid]->GetBinContent( ix+1, irap+1, nBinPhi ) );
        hist_pt_rap_phi_real[ipid]->SetBinContent( ix+1, irap+1, 1,       avgReal );
        hist_pt_rap_phi_real[ipid]->SetBinContent( ix+1, irap+1, nBinPhi, avgReal );

        double avgMeas = ( hist_pt_rap_phi_meas[ipid]->GetBinContent( ix+1, irap+1, 1 ) + hist_pt_rap_phi_meas[ipid]->GetBinContent( ix+1, irap+1, nBinPhi ) );
        hist_pt_rap_phi_meas[ipid]->SetBinContent( ix+1, irap+1, 1,       avgMeas );
        hist_pt_rap_phi_meas[ipid]->SetBinContent( ix+1, irap+1, nBinPhi, avgMeas );
      }
    }
  }
  // Compensating the end of the bins



  if( fUseTotalQ )
    ConstructTMfromQ();
  else
  {
    for( int ipid=0; ipid<nParticleType; ipid++ )
    {
      for( int ipt=0; ipt<nBinPt; ipt++ )
      {
        for( int irap=0; irap<nBinRap; irap++ )
        {
          auto hist_proj = hist_pt_rap_phi_real[ipid]->ProjectionZ("",ipt+1,ipt+1, irap+1,irap+1);

          for( int ireal=0; ireal<nBinPhi; ireal++ )
          {
            double normFactor = hist_proj->GetBinContent( ireal+1 );
            for( int imeas=0; imeas<nBinPhi; imeas++ )
            {
              // Normalize the blurring with the real reaction plane direction counts
              double binContent = hist_phi_diff[ipid][ipt][irap]->GetBinContent( ireal+1, imeas+1 );
              if( 0.<normFactor ) hist_phi_diff[ipid][ipt][irap]->SetBinContent( ireal+1, imeas+1, binContent / normFactor );
            }
          }
        }
      }
    }
  }


  return;
}

void SimulationTools::CallParticle()
{
  fVector->resize(0);
  fZ->resize(0);
  fA->resize(0);
  fMass->resize(0);

  if( fModel != nullptr )
  {
    fModel->EmitParticlesInCM( fZ, fA, fMass, fVector );
    return;
  }
  else if( fDataTree!=nullptr )
  {
    for( int iTrack=0; iTrack<tree_nTracks; iTrack++ )
    {
      int ipid = -1;
      if( tree_pid[iTrack]==2212 ) ipid = 0; // proton
      if( tree_pid[iTrack]==1000010020 ) ipid = 1; // deuteron
      if( tree_pid[iTrack]==1000010030 ) ipid = 2; // triton
      if( tree_pid[iTrack]==1000020030 ) ipid = 3; // 3He
      if( tree_pid[iTrack]==1000020040 ) ipid = 4; // alpha

      //double mass = tree_mass[iTrack]; // GeV/c2
      double mass = particle_mass[ipid]; // GeV/c2
      double px = tree_px[iTrack]; // GeV/c
      double py = tree_py[iTrack]; // GeV/c
      double pz = tree_pz[iTrack]; // GeV/c
      double scale_factor = std::pow( 10.0, 7 );
      px   = std::round(px * scale_factor) / scale_factor;
      py   = std::round(py * scale_factor) / scale_factor;
      pz   = std::round(pz * scale_factor) / scale_factor;

      // -----------------------------------------------------------------------
      // Re-calculate total energy and assign into LorentzVector
      double totalE_lab = TMath::Sqrt( px*px + py*py + pz*pz + mass*mass );
      TLorentzVector L_lab( px, py, pz, totalE_lab );
      double ylab = L_lab.Rapidity(); // LAB rapidity
      double yCM  = ylab - y_cm;  // NN Center of Mass rapidity
      double y0   = yCM / y_cm;   // normalized rapidity

      double mt      = TMath::Sqrt( mass*mass + px*px+py*py );
      double totalE  = mt * TMath::CosH( y0 );
             pz      = mt * TMath::SinH( y0 );
      TLorentzVector L( px, py, pz, totalE );
      // Re-calculate total energy and assign into LorentzVector
      // -----------------------------------------------------------------------

      if( ipid==-1 ) continue; // To avoid particles not identified

      fZ->push_back( particle_z[ipid] );
      fA->push_back( particle_a[ipid] );
      fMass->push_back( mass );
      fVector->push_back( L );
    }

    return;
  }
  else if( fHistogramRestored.size() != 0 )
  {
    for( int ipid=0; ipid<nParticleType; ipid++ )
    {
      double mass  = particle_mass[ipid];

      double pt = 0.;
      double pt_width = fHistogramRestored[ipid]->GetXaxis()->GetBinWidth( 1 );
      double rapidity_width = fHistogramRestored[ipid]->GetYaxis()->GetBinWidth( 1 );
      double phi_width = fHistogramRestored[ipid]->GetZaxis()->GetBinWidth( 1 );

      for( auto iBin : fCheckedBin[ipid] )
      {
        int binx, biny, binz;
        fHistogramRestored[ipid]->GetBinXYZ( iBin, binx, biny, binz );

        double pt_center = fHistogramRestored[ipid]->GetXaxis()->GetBinCenter( binx );
        double rapidity_center = fHistogramRestored[ipid]->GetYaxis()->GetBinCenter( biny );
        double phi_center = fHistogramRestored[ipid]->GetZaxis()->GetBinCenter( binz );

        double binContent = fHistogramRestored[ipid]->GetBinContent( iBin );
        binContent = binContent * pt_center; // Since 3D yield is dN/ptdptdydphi
        binContent = binContent * pt_width * rapidity_width * phi_width;
        int nCount = 0;
        if( !binContent ) nCount = 0;
        else if( binContent < 1e-5 ) { if( gRandom->Uniform() < binContent ) nCount = 1; }
        else nCount = gRandom->Poisson( binContent );


        for( int i=0; i<nCount; i++ )
        {
          double pt = pt_center + pt_width*(gRandom->Uniform() - 0.5);
          double rapidity = rapidity_center + rapidity_width*(gRandom->Uniform() - 0.5); // This is reduced rapiditiy

          // Pawel's implement
          int IPE1 = binx-1;
          int IPE21 = binx+IPE1;
          int IPE1K = IPE1*IPE1;
          double PTK = pt_width*pt_width * ( IPE1K + IPE21 * gRandom->Uniform() );
          pt = std::sqrt( PTK );

          double phi = phi_center + phi_width*(gRandom->Uniform() - 0.5);
          if( phi < -TMath::Pi() ) phi += TMath::TwoPi();
          if( phi >  TMath::Pi() ) phi -= TMath::TwoPi();

          double px   = pt * TMath::Cos( phi );
          double py   = pt * TMath::Sin( phi );

          double mt    = TMath::Sqrt( pt*pt + mass*mass );
          double totE  = mt * TMath::CosH( rapidity );
          double pz    = mt * TMath::SinH( rapidity );
          TLorentzVector L( px, py, pz, totE );

          fZ->push_back( particle_z[ipid] );
          fA->push_back( particle_a[ipid] );
          fMass->push_back( mass );
          fVector->push_back( L );
        }
      }
    }

    return;
  }
}



double SimulationTools::GetEfficiency( int ipid, double pt, double y0, double phiLab )
{
  TH3D* thisEff = fEfficiency_hist[ipid];
  if( thisEff==nullptr ) return 0.;


  int nBinPtEff = thisEff->GetNbinsX();
  double binWidthPtEff = thisEff->GetXaxis()->GetBinWidth( 1 );
  if( nBinPtEff < thisEff->GetXaxis()->FindBin( pt ) ) return 0.;
  if( thisEff->GetXaxis()->FindBin( pt ) < 1 ) return 0.;

  int nBinRapEff = thisEff->GetNbinsY();
  double binWidthRapEff = thisEff->GetYaxis()->GetBinWidth( 1 );
  double rapMinEff = thisEff->GetYaxis()->GetXmin();
  if( nBinRapEff < thisEff->GetYaxis()->FindBin( y0 ) ) return 0.;
  if( thisEff->GetYaxis()->FindBin( y0 ) < 1 ) return 0.;

  int nBinPhiEff = thisEff->GetNbinsZ();
  double binWidthPhiEff = thisEff->GetZaxis()->GetBinWidth( 1 );

  double effi = 0.;

  double ypnr = (y0 + std::abs(rapMinEff))/binWidthRapEff + 1.;                                                                                                                                                                                                                 
  ypnr = ypnr-.5 + 1e-20;
  int iyb = TMath::Floor(ypnr);
  double fyf = ypnr-iyb;
  double fyb = 1.-fyf;
  int iyf = iyb+1;
  iyb = std::max(0,iyb);
  iyf = std::min(nBinRapEff,iyf);

  double ptr  = pt/binWidthPtEff + 1.;
  ptr  = ptr-.5 + 1e-20;
  int iptb = TMath::Floor(ptr);
  double fptf = ptr-iptb;
  double fptb = 1.-fptf;
  int iptf = iptb+1;
  iptb = std::max(0,iptb);
  iptf = std::min(nBinPtEff,iptf);

  double phim  = phiLab + 180.;
  double phir  = phim/binWidthPhiEff + 1.;
  phir  = phir-.5 + 1e-20;
  int iphib = TMath::Floor(phir);
  double fphif = phir-iphib;
  double fphib = 1.-fphif;
  int iphif = iphib+1;
  if( iphib < 1 ) iphib = nBinPhiEff;
  if( iphif > nBinPhiEff ) iphif = 1;


  double efiBBB = thisEff->GetBinContent(iptb,iyb,iphib);
  double efiBFB = thisEff->GetBinContent(iptb,iyf,iphib);
  double efiFBB = thisEff->GetBinContent(iptf,iyb,iphib);
  double efiFFB = thisEff->GetBinContent(iptf,iyf,iphib);
  double efiBBF = thisEff->GetBinContent(iptb,iyb,iphif);
  double efiBFF = thisEff->GetBinContent(iptb,iyf,iphif);
  double efiFBF = thisEff->GetBinContent(iptf,iyb,iphif);
  double efiFFF = thisEff->GetBinContent(iptf,iyf,iphif);

  effi = fphib*(fptb*(fyb*efiBBB + fyf*efiBFB)
               +fptf*(fyb*efiFBB + fyf*efiFFB))
        +fphif*(fptb*(fyb*efiBBF + fyf*efiBFF)
               +fptf*(fyb*efiFBF + fyf*efiFFF));

  return effi;
}

void SimulationTools::ConstructTMfromQ()
{
  auto hProj = hist_totalq->ProjectionX("", 1, nBinQphi-1);
  for( int ireal=0; ireal<nBinQphi-1; ireal++ )
  {
    double content = hProj->GetBinContent( ireal+1 ); // Had been integrated through projection

    if( content!=0. )
      for( int imeas=0; imeas<nBinQphi-1; imeas++ )
        hist_totalq->SetBinContent( ireal+1, imeas+1, hist_totalq->GetBinContent(ireal+1, imeas+1)/content );
  }


  for( int ipid=0; ipid<nParticleType; ipid++ )
  {
    auto hproj = fHistogramRestored[ipid]->Project3D("XY");
    //double mass = particle_mass[ipid];
    //double AA   = particle_a[ipid];
    for( int ipt=0; ipt<nBinPt; ipt++ )
    {
      double pt_center = binWidthPt[ipid]*(ipt+0.5);
      for( int irap=0; irap<nBinRap; irap++ )
      {
        // ------------------------------------------------------------
        // If there is no data in 2D distribution (pt vs y)
        double val_2d = hproj->GetBinContent( irap+1, ipt+1 );
        if( !val_2d ) continue; // Then skip to save the time
        // ------------------------------------------------------------

        double rap_center = hist_pt_rap_phi_rest[ipid]->GetYaxis()->GetBinCenter( irap+1 );

        int nsamp = 40;
        for( int isamp=0; isamp<nsamp; isamp++ )
        {
          // Real phi angle
          for( int imeas=0; imeas<nBinPhi-1; imeas++ )
          {
            // Particle angle rel. to rp.
            //double phi_meas = hist_phi_diff[ipid][ipt][irap]->GetXaxis()->GetBinCenter( imeas+1 ) + gRandom->Uniform(-binWidthPhy/2., binWidthPhy/2.);
            double phi_meas = ( imeas-nBinPhy )*binWidthPhy + gRandom->Uniform(-binWidthPhy/2., binWidthPhy/2.);
            if( phi_meas < -TMath::Pi() ) phi_meas = phi_meas + TMath::TwoPi();
            if( phi_meas >  TMath::Pi() ) phi_meas = phi_meas - TMath::TwoPi();

            for( int imeas_rp=0; imeas_rp<nBinQphi-1; imeas_rp++ )
            {
              //double phi_meas_rp = hist_central->GetBinCenter( imeas_rp+1 ) + gRandom->Uniform(-binWidthQphy/2., binWidthQphy/2.);
              double phi_meas_rp = ( imeas_rp-nBinQphy )*binWidthQphy + gRandom->Uniform(-binWidthQphy/2., binWidthQphy/2.);
              if( phi_meas_rp < -TMath::Pi() ) phi_meas_rp = phi_meas_rp + TMath::TwoPi();
              if( phi_meas_rp >  TMath::Pi() ) phi_meas_rp = phi_meas_rp - TMath::TwoPi();

              double phi_lab = phi_meas + phi_meas_rp;
              if( phi_lab < -TMath::Pi() ) phi_lab = phi_lab + TMath::TwoPi();
              if( phi_lab >  TMath::Pi() ) phi_lab = phi_lab - TMath::TwoPi();

              // Sampling random pt and rap
              //double pt  = pt_center  + gRandom->Uniform( -binWidthPt[ipid]/2.,  binWidthPt[ipid]/2. );
              double pt  = 0.;
              double rap = rap_center + gRandom->Uniform( -binWidthRap/2., binWidthRap/2. );


              int IPE1 = ipt+1-1;
              int IPE21 = ipt+1+IPE1;
              int IPE1K = IPE1*IPE1;
              double PTK = binWidthPt[ipid]*binWidthPt[ipid] * ( IPE1K + IPE21 * gRandom->Uniform() );
              pt = std::sqrt( PTK );

              double effi = 0.;
              //int thisBin = 0;
              //thisBin = fEfficiency_vector[ipid]->FindFixBin( pt*1e+3, rap, phi_lab*TMath::RadToDeg() );
              //effi = fEfficiency_vector[ipid]->GetEfficiency( thisBin );
              effi = GetEfficiency( ipid, pt*1e+3, rap, phi_lab*TMath::RadToDeg() );

              if( effi!=0. )
              {
                for( int ireal_rp=0; ireal_rp<nBinQphi-1; ireal_rp++ )
                {
                  //double phi_real_rp = hist_central->GetBinCenter( ireal_rp+1 ) + gRandom->Uniform(-binWidthQphy/2., binWidthQphy/2.);
                  double phi_real_rp = ( ireal_rp-nBinQphy )*binWidthQphy + gRandom->Uniform(-binWidthQphy/2., binWidthQphy/2.);
                  if( phi_real_rp < -TMath::Pi() ) phi_real_rp = phi_real_rp + TMath::TwoPi();
                  if( phi_real_rp >  TMath::Pi() ) phi_real_rp = phi_real_rp - TMath::TwoPi();

                  double phi_real = phi_lab - phi_real_rp;
                  if( phi_real < -TMath::Pi() ) phi_real = phi_real + TMath::TwoPi();
                  if( phi_real >  TMath::Pi() ) phi_real = phi_real - TMath::TwoPi();

                  int phiBin = hist_phi_diff[ipid][ipt][irap]->GetXaxis()->FindBin( phi_real );
                  if( phiBin==nBinPhi ) phiBin=1;
                  double prev_val = hist_phi_diff[ipid][ipt][irap]->GetBinContent( phiBin, imeas+1 );

                  double total_val = hist_totalq->GetBinContent( ireal_rp+1, imeas_rp+1 );
                  double val = prev_val + total_val*effi;
                  hist_phi_diff[ipid][ipt][irap]->SetBinContent( phiBin, imeas+1, val );
                }
              }

            }
          }
        }

        hist_phi_diff[ipid][ipt][irap]->Scale( 1./(nBinQphy*2 * nsamp) );

        for( int iphi=0; iphi<nBinPhi; iphi++ )
        {
          double val = hist_phi_diff[ipid][ipt][irap]->GetBinContent( iphi+1, 1 );
          hist_phi_diff[ipid][ipt][irap]->SetBinContent( iphi+1, nBinPhi, val );

          val = hist_phi_diff[ipid][ipt][irap]->GetBinContent( 1, iphi+1 );
          hist_phi_diff[ipid][ipt][irap]->SetBinContent( nBinPhi, iphi+1, val );
        }

      }
    }
  }

}

void SimulationTools::NormalizeEstimatedDistribution()
{
  for( int ipid=0; ipid<nParticleType; ipid++ )
    for( int ipt=0; ipt<nBinPt; ipt++ )
      for( int iy=0; iy<nBinRap; iy++ )
        for( int iphi=0; iphi<nBinPhi; iphi++ )
        {
          double err = fHistogramMeasured[ipid]->GetBinError( ipt+1, iy+1, iphi+1 );
          if( err )
          {
            double val = fHistogramMeasured[ipid]->GetBinContent( ipt+1, iy+1, iphi+1 );
            double ANS    = val/err;
            double AN     = ANS*ANS;
            double factor = val/AN;
            int multi     = gRandom->Poisson( AN );
            hist_pt_rap_phi_meas[ipid]->SetBinContent( ipt+1, iy+1, iphi+1, factor*multi );
          }
        }
}

void SimulationTools::EstimaedError( int nSample )
{
  for( int ipid=0; ipid<nParticleType; ipid++ )
    for( int ipt=0; ipt<nBinPt; ipt++ )
      for( int iy=0; iy<nBinRap; iy++ )
        for( int iphi=0; iphi<nBinPhi; iphi++ )
        {
          double err = hist_pt_rap_phi_rest[ipid]->GetBinError( ipt+1, iy+1, iphi+1 );
          double val = hist_pt_rap_phi_rest[ipid]->GetBinContent( ipt+1, iy+1, iphi+1 );
          double val_org = fHistogramRestored[ipid]->GetBinContent( ipt+1, iy+1, iphi+1 );
          err = err*err * nSample + (val-val_org)*(val-val_org);
          err = TMath::Sqrt( err / (nSample+1) );

          hist_pt_rap_phi_rest[ipid]->SetBinError( ipt+1, iy+1, iphi+1, err );
          if( nSample%4==0 ) hist_pt_rap_phi_rest[ipid]->SetBinContent( ipt+1, iy+1, iphi+1, val_org ); // To stablize?
        }
}

void SimulationTools::RLdeblurring()
{
  bool reg_on = 1;
  double gamma_factor = 1.00;

  for( int ipid=0; ipid<nParticleType; ipid++ )
    for( int ipt=0; ipt<nBinPt; ipt++ )
      for( int iy=0; iy<nBinRap; iy++ )
        for( int iphi=0; iphi<nBinPhi; iphi++ )
        {
          double val = hist_pt_rap_phi_rest[ipid]->GetBinContent( ipt+1, iy+1, iphi+1 );
          val = TMath::Max( val, 1e-10 );

          hist_pt_rap_phi_rest[ipid]->SetBinContent( ipt+1, iy+1, iphi+1, val );
        }


  int nIteration = 100;

  for( int IT=0; IT<nIteration; IT++ )
  {
    // ----------------------- START OF nIteration -----------------------------
    for( int ipid=0; ipid<nParticleType; ipid++ )
    {
      for( int iBinPt=0; iBinPt<nBinPt; iBinPt++ )
      {
        for( int iBinRap=0; iBinRap<nBinRap; iBinRap++ )
        {

          for( int imeas=0; imeas<nBinPhi; imeas++ )
          {
            double nrj = 0;
            for( int ireal=0; ireal<nBinPhi-1; ireal++ )
            {
              double lhs = hist_phi_diff[ipid][iBinPt][iBinRap]->GetBinContent( ireal+1, imeas+1 );
              double rhs = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1, iBinRap+1, ireal+1 );

              nrj += lhs*rhs;
            }

            // n^r_j, acting P_ji on N^r
            hist_pt_rap_phi_blur[ipid]->SetBinContent( iBinPt+1, iBinRap+1, imeas+1, nrj );
          }
        }
      }
    }

    for( int ipid=0; ipid<nParticleType; ipid++ )
    {
      for( int iBinPt=0; iBinPt<nBinPt; iBinPt++ )
      {
        for( int iBinRap=0; iBinRap<nBinRap; iBinRap++ )
        {
          double vm1, v_0, vp1;
          for( int ireal=0; ireal<nBinPhi; ireal++ )
          {
            double ari = 0;
            double weighted_prob = 0;
            for( int imeas=0; imeas<nBinPhi-1; imeas++ )
            {
              double lhs1 = 0;
              lhs1 = hist_pt_rap_phi_meas[ipid]->GetBinContent( iBinPt+1, iBinRap+1, imeas+1 ); // Estimated
              double lhs2 = hist_pt_rap_phi_blur[ipid]->GetBinContent( iBinPt+1, iBinRap+1, imeas+1 ); // One more blurred
              double lhs = 0;
              if( lhs2 != 0 ) lhs = lhs1 / lhs2;
              double rhs = hist_phi_diff[ipid][iBinPt][iBinRap]->GetBinContent( ireal+1, imeas+1 );

              ari += lhs*rhs;
              if( lhs2 != 0 ) weighted_prob += rhs;
            }

            if( ari!=0 && weighted_prob!=0 )
            {
              double regul = 1.;

              int binPhiM1 = ireal-1 + 1;
              int binPhi0  =   ireal + 1;
              int binPhiP1 = ireal+1 + 1;

              if( binPhiM1 < 1 )         binPhiM1 = nBinPhi - 1;
              if( nBinPhi  <= binPhiP1 ) binPhiP1 = 1;

              // Switch to Pawels modulo
              binPhiM1 = modulo( ireal-nBinPhy -1+nBinPhy, 2*nBinPhy ) - nBinPhy;
              binPhiP1 = modulo( ireal-nBinPhy +1+nBinPhy, 2*nBinPhy ) - nBinPhy;
              binPhiM1 = binPhiM1 + nBinPhy + 1;
              binPhiP1 = binPhiP1 + nBinPhy + 1;

              // Consider 3D regularization
                double scaleLambda = 0.;

                // Phi regularization
                vm1 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1, iBinRap+1, binPhiM1 );
                v_0 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1, iBinRap+1, binPhi0 );
                vp1 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1, iBinRap+1, binPhiP1 );
                scaleLambda += (v_0 < vm1) ? -1 : 1;
                scaleLambda += (v_0 < vp1) ? -1 : 1;

                // rapidity regularization
                if( 1 < iBinRap && iBinRap < nBinRap-1 )
                {
                  vm1 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1, iBinRap+1 -1, binPhi0 );
                  v_0 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1, iBinRap+1   , binPhi0 );
                  vp1 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1, iBinRap+1 +1, binPhi0 );
                  scaleLambda += (v_0 < vm1) ? -0.5 : 0.5;
                  scaleLambda += (v_0 < vp1) ? -0.5 : 0.5;
                }

                // momentum regularization
                if( iBinPt < nBinPt )
                {
                  v_0 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1   , iBinRap+1, binPhi0 );
                  vp1 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1 +1, iBinRap+1, binPhi0 );
                  scaleLambda += (v_0 < vp1) ? -0.5 : 0.5;

                  if( 0 < iBinPt )
                  {
                    vm1 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1 -1, iBinRap+1, binPhi0 );
                    scaleLambda += (v_0 < vm1) ? -0.5 : 0.5;
                  }
                  else
                  {
                    // Switch to Pawels modulo
                    int binComp = (ireal - nBinPhy);
                    int thisPhiBin = modulo( binComp, nBinPhy*2 ) - nBinPhy;
                    thisPhiBin = thisPhiBin + nBinPhy + 1;
                    vm1 = hist_pt_rap_phi_rest[ipid]->GetBinContent( iBinPt+1, iBinRap+1, thisPhiBin );
                    scaleLambda += (v_0 < vm1) ? -0.5 : 0.5;
                  }
                }

                if( reg_on ) 
                {
                  regul = 1. / (1. + 0.5*scaleLambda*fLambda);
                  if( ireal==nBinPhy || ireal==0 || ireal==(nBinPhi-1) ) regul = 1. / (1. + 0.25*scaleLambda*fLambda);
                }


              double final_val = 0.;
              if( weighted_prob>0. ) final_val = TMath::Power( ari/(weighted_prob+1e-12), gamma_factor ) * regul * v_0;

              hist_pt_rap_phi_temp[ipid]->SetBinContent( iBinPt+1, iBinRap+1, ireal+1, final_val );
            }
          }

        } // Rapidity Loop Over
      } // Pt loop
    } // pID loop

    for( int ipid=0; ipid<nParticleType; ipid++ )
    {
      for( int iBinPt=0; iBinPt<nBinPt; iBinPt++ )
      {
        for( int iBinRap=0; iBinRap<nBinRap; iBinRap++ )
        {
          hist_pt_rap_phi_temp[ipid]->SetBinContent( iBinPt+1, iBinRap+1, nBinPhi, hist_pt_rap_phi_temp[ipid]->GetBinContent(iBinPt+1, iBinRap+1, 1) );
          for( int i=0; i<nBinPhi; i++ ) hist_pt_rap_phi_rest[ipid]->SetBinContent( iBinPt+1, iBinRap+1, i+1, hist_pt_rap_phi_temp[ipid]->GetBinContent(iBinPt+1, iBinRap+1, i+1) );
          for( int i=0; i<nBinPhi; i++ ) hist_pt_rap_phi_rest[ipid]->SetBinError( iBinPt+1, iBinRap+1, i+1, 0 );
        }
      }
    }

  } // Total iteration loop

}

void SimulationTools::SetupOptimal()
{
  // If you are trying to initially calculate the optimal weight,
  if( fDataTree || fModel || fUpdateInternal )
  {
    fHistogramMeasured = hist_pt_rap_phi_meas;
    fHistogramRestored = hist_pt_rap_phi_rest;
  }

  // Optimal weight calculations
  optWeight.resize( nParticleType );
  for( int ipid=0; ipid<nParticleType; ipid++ ) optWeight.at(ipid).resize( nBinRap ); // Slice pT bins

  for( int ipid=0; ipid<nParticleType; ipid++ )
  {
    for( int iy=0; iy<nBinRap; iy++ )
      for( int iphi=0; iphi<nBinPhi-1; iphi++ ) // The both end bins are the same bin
      {
        int ipt_wipe = 0;
        double val_bk = fHistogramRestored[ipid]->GetBinContent( 1, iy+1, iphi+1 );
        for( int ipt=1; ipt<nBinPt; ipt++ )
        {
          double val_inc = fHistogramRestored[ipid]->GetBinContent( ipt+1, iy+1, iphi+1 );
          if( 3*val_bk < val_inc ) ipt_wipe=ipt;
          else val_bk = val_inc;
        }

        for( int ipt=ipt_wipe; 0&& ipt_wipe && ipt<nBinPt; ipt++ )
          fHistogramRestored[ipid]->SetBinContent( ipt+1, iy+1, iphi+1, 0. );
      }
  }

  for( int ipid=0; ipid<nParticleType; ipid++ ) 
  {
    double totalSum = 0;
    std::vector< double > dnSumMeas;
    dnSumMeas.resize(nBinRap);
    for( int iy=0; iy<nBinRap; iy++ )
    {
      double dnSum   = 0;
      double dpxSum  = 0;
      double dpt2Sum = 0;

      double dnSum_meas   = 0;
      for( int ipt=0; ipt<nBinPt; ipt++ )
      {
        double pt = fHistogramRestored[ipid]->GetXaxis()->GetBinCenter( ipt+1 );
        double ptf = pt + 0.5*binWidthPt[ipid];
        double ptb = pt - 0.5*binWidthPt[ipid];

        double fapd   =  ( TMath::Power(ptf,2.) - TMath::Power(ptb,2.) )/2.;
        double fapdpt =  ( TMath::Power(ptf,3.) - TMath::Power(ptb,3.) )/3.;
        double fapdptk = ( TMath::Power(ptf,4.) - TMath::Power(ptb,4.) )/4.;

        //int nPhiCut = 0;
        double dnSumPhi = 0;
        double dphiSum   = 0;

        double dnSumPhi_meas = 0;
        for( int iphi=0; iphi<nBinPhi-1; iphi++ ) // The both end bins are the same bin
        {
          double phi = fHistogramRestored[ipid]->GetZaxis()->GetBinCenter( iphi+1 );
          double val = fHistogramRestored[ipid]->GetBinContent( ipt+1, iy+1, iphi+1 );
          if( dnCut < val )
          {
            //nPhiCut++;

            double phif = double( iphi-nBinPhy+0.5 ) * binWidthPhy;
            double phib = double( iphi-nBinPhy-0.5 ) * binWidthPhy;
            double cosint = TMath::Sin( phif ) - TMath::Sin( phib );

            dnSumPhi = dnSumPhi + val;
            dphiSum  = dphiSum + val*cosint;

            double measVal = fHistogramMeasured[ipid]->GetBinContent( ipt+1, iy+1, iphi+1 );
            dnSumPhi_meas = dnSumPhi_meas + measVal;
          }
        }

        dnSum   = dnSum  + dnSumPhi*binWidthPhy*fapd;;
        dpxSum  = dpxSum + dphiSum*fapdpt;
        dpt2Sum = dpt2Sum + dnSumPhi*binWidthPhy*fapdptk;

        // For cut (avg)
        dnSum_meas = dnSum_meas + dnSumPhi_meas*binWidthPhy*fapd;
      }

      if( dnSum )
      {
        double pxr = dpxSum / dnSum;
        double pt2r = dpt2Sum / dnSum;
        double den = pt2r - pxr*pxr; // <pt>^2 - <px>^2

        if( den )
        {
          double wgr = pxr / den;
          optWeight[ipid][iy] = wgr;
        }
      }

      if( (nBinRap-1)/2-1 <= iy ) totalSum = totalSum + dnSum;
      dnSumMeas[iy] = dnSum_meas;
    }

    double avgSum = totalSum / ( (nBinRap-1)/2 );
    double avgCut = avgSum / 40.;

    for( int iy=0; iy<nBinRap; iy++ )
      if( dnSumMeas[iy] < avgCut ) optWeight[ipid][iy] = 0;

    for( int iy=0; iy<nBinRap; iy++ ) cout << ipid << " " << iy << " " << optWeight[ipid][iy] << endl;

  }
}

void SimulationTools::SaveOutputFiles()
{
  auto file = new TFile( fOutputName, "RECREATE" );
  file->cd();
  hist_central->Write();
  hist_totalq->Write();

  for( int ipid=0; ipid<nParticleType; ipid++ ) for( int irap=0; irap<nBinRap; irap++ ) for( int ipt=0; ipt<nBinPt; ipt++ ) hist_phi_diff[ipid][ipt][irap]->Write();
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_real[ipid]->Write();
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_meas[ipid]->Write();
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_rest[ipid]->Write();
  hist_rp_phi_diff->Write();

  file->Close();
}


