double particle_mass[5] = {0};

void PrepareMass()
{
  double be = .008;
  double mass_proton     = .9383; 
  double mass_neutron    = .9396; 
  double mass_avg        = .5*(mass_proton+mass_neutron);
  double mass_minus_be   = mass_avg - be;
  //mass_avg_square = mass_avg*mass_avg;

  //  DEUTERON
  double mass_avg2       = mass_avg + mass_avg;
  double be_d            = .002225; 
  double mass_a2 = mass_avg + mass_avg - be_d;

  //  A=3
  double mass_avg3         = mass_avg2 + mass_avg;
  double be_t              = .0086; 
  double be_3he            = .0080; 
  double mass_minus_be_t   = mass_proton + mass_neutron + mass_neutron - be_t;
  double mass_minus_be_3he = mass_proton + mass_proton + mass_neutron - be_3he;
  double mass_a3           = .5*(mass_minus_be_t+mass_minus_be_3he);

  //  ALPHA
  double mass_avg4 = mass_avg3 + mass_avg;
  double be_4he    = .0286; 
  double mass_a4   = mass_avg4 - be_4he;

  particle_mass[0] = mass_avg;
  particle_mass[1] = mass_a2;
  particle_mass[2] = mass_a3;
  particle_mass[3] = mass_a3;
  particle_mass[4] = mass_a4;
}

void draw()
{
  PrepareMass();
  double y_cm = 0.372156;

  auto file1 = new TFile( "eff_Sn108_midcentral.root", "READ" );
  auto file2 = new TFile( "eff_Sn108_midcentral.root_backup", "READ" );

  auto e1 = (TEfficiency*) file1->Get( "e3_alpha_HighRes" );
  auto e2 = (TEfficiency*) file2->Get( "e3_alpha_HighRes" );

  auto h1 = (TH3D*) e1->GetPassedHistogram();
  auto h2 = (TH3D*) e2->GetPassedHistogram();

  auto hist_angle = new TH2D( "hist_angle", "", 400, -200, 200, 200, 0, 200 );
  auto hist_phase = new TH2D( "hist_phase", "", 400, -2*y_cm  , 2*y_cm  , 200, 0, 2   );
  auto hist_angle_det = new TH2D( "hist_angle_det", "", 400, -200, 200, 200, 0, 200 );
  auto hist_phase_det = new TH2D( "hist_phase_det", "", 400, -2*y_cm  , 2*y_cm  , 200, 0, 2   );

  int nEvent = 1e+6;
  for( int ieve=0; ieve<nEvent; ieve++ )
  {
    double mass = particle_mass[4];
    double pt  =  gRandom->Uniform( 0., 2. ); // GeV/c
    double rap = gRandom->Uniform( -2., 2. ); // scaled
    double phi = gRandom->Uniform( -180., 180. ); // deg

    double rap_lab = (rap+1.) * y_cm;
    double mt = TMath::Sqrt( pt*pt + mass*mass );
    double pz_lab = mt * TMath::SinH( rap_lab );
    double theta = TMath::ATan2( pt, pz_lab );

    hist_phase->Fill( rap*y_cm, pt );
    hist_angle->Fill( phi, theta*TMath::RadToDeg() );

    int thisBin = e2->FindFixBin( pt*1e+3, rap, phi );
    double eff = e2->GetEfficiency( thisBin );

    if( gRandom->Uniform() < eff )
    {
      hist_phase_det->Fill( rap*y_cm, pt );
      hist_angle_det->Fill( phi, theta*TMath::RadToDeg() );
    }
  }

  hist_phase_det->Divide( hist_phase );
  hist_angle_det->Divide( hist_angle );

  /*
  auto c = new TCanvas( "c", "", 1400, 600 );
  c->Divide( 2, 1 );
  c->cd( 1 );
  hist_phase_det->Smooth();
  hist_phase_det->Smooth();
  hist_phase_det->Draw( "colz" );
  c->cd( 2 );
  hist_angle_det->Smooth();
  hist_angle_det->Smooth();
  hist_angle_det->Draw( "colz" );
  */

  auto outputFile = new TFile( "eff.root", "RECREATE" );
  hist_phase_det->Write();
  hist_angle_det->Write();


  /*
     auto c = new TCanvas( "c", "", 1400, 600 );
     c->Divide( 2, 1 );
     c->cd( 1 );
     h1->Project3D("XY")->Draw("colz");
     c->cd( 2 );
     h2->Project3D("XY")->Draw("colz");

     return;
     gStyle->SetOptStat( 0 );

     auto hist1 = (TH2D*) h1->Project3D("XY");
     auto hist2 = (TH2D*) h2->Project3D("YZ");

     auto c1 = new TCanvas( "c1", "", 700, 600 );
     c1->SetLeftMargin( 0.16 );
     c1->SetRightMargin( 0.04 );
     c1->SetBottomMargin( 0.16 );
     c1->SetTopMargin( 0.04 );
     hist1->SetTitle( "" );
     hist1->GetXaxis()->SetTitle( "y_{0}" );
     hist1->GetYaxis()->SetTitle( "p_{#perp}" );
     hist1->Draw( "col" );

     auto c2 = new TCanvas( "c2", "", 700, 600 );
     c2->SetLeftMargin( 0.16 );
     c2->SetRightMargin( 0.04 );
     c2->SetBottomMargin( 0.16 );
     c2->SetTopMargin( 0.04 );
     hist2->SetTitle( "" );
     hist2->GetXaxis()->SetTitle( "#phi_{lab}" );
     hist2->GetYaxis()->SetTitle( "y_{0}" );
     hist2->Draw( "col" );
   */
}
