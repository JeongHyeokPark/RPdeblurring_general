TString particleName[6] = { "proton", "deuteron", "triton", "3He", "alpha", "neutron" };

bool saveFigure = 1;

bool drawFlow = 1;
bool v_one = 0;

bool ystar = 1;

bool expdata = 1;

//double xpa[21] = { -0.600 , -0.500 , -0.400 , -0.300 , -0.200 , -0.100 , 0.000 , 0.100 , 0.200 , 0.300 , 0.400 , 0.500 , 0.600 , 0.700 , 0.800 , 0.900 , 1.000 , 1.100 , 1.200 , 1.300 , 1.400 };// Even fine binning
//double ypa[21] = { -0.2826, 0.0110, -0.0290, -0.0515, -0.0763, -0.0832, -0.0843, -0.0972, -0.0878, -0.0682, -0.0467, -0.0074, 0.0382, 0.0732, 0.1154, 0.1318, 0.1190, 0.1394, 0.1305, 0.2014, 0.323};// Even fine binning
//double xpa[21] = { -2.000 , -1.800 , -1.600 , -1.400 , -1.200 , -1.000 , -0.800 , -0.600 , -0.400 , -0.200 , 0.000 , 0.200 , 0.400 , 0.600 , 0.800 , 1.000 , 1.200 , 1.400 , 1.600 , 1.800 , 2.000 }; // Uneven wider binning
//double ypa[21] = { 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, -0.0081, -0.0560, -0.0923, -0.1099, -0.0851, -0.0468, 0.0355, 0.0886, 0.1109, 0.1180, 0.9836, 0.0000, 0.0000, 0.0000 }; // Uneven wider binning
double xpa[21] = { -2.000, -1.800, -1.600, -1.400, -1.200, -1.000, -0.800, -0.600, -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000, 1.200, 1.400, 1.600, 1.800, 2.000 }; // Even wider binning
double ypa[21] = { 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0073, -0.0414, -0.0792, -0.0927, -0.0861, -0.0418, 0.0324, 0.0909, 0.1048, 0.1271, 0.5826, 0.0000, 0.0000, 0.0000 }; // Even wider binning

auto gpa = new TGraph( 21, xpa, ypa );


// numUp: The updated number of (r) for N^(r)
// nPidTest: pid number             for preview of the restored distribution
// nPtTest: number of pt bin        for preview of the restored distribution
// nRapTest: number of rapidity bin for preview of the restored distribution
void draw( int numUp=0, int nPidTest=4, int nPtTest=2, int nRapTest=18 )
{
  static const int nParticleType = 5;

  auto file1 = new TFile( "./input.root", "READ" );
  if( !saveFigure ) file1 = new TFile( "./output_err.root", "READ" );
  auto file_exp = new TFile( "/home/Symmetry/jhpark/spirit/ExpFlowData_Sn108_midcentral_opt.root", "READ" );

  TH3D* hist_pt_rap_phi_real[nParticleType];
  TH3D* hist_pt_rap_phi_meas_thefirst[nParticleType];
  TH3D* hist_pt_rap_phi_rest[nParticleType];

  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_real[ipid] = (TH3D*) file1->Get( Form("hist_pt_rap_phi_real_%s", particleName[ipid].Data()) );
  if( expdata ) 
    for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_real[ipid] = (TH3D*) file_exp->Get( Form("hist_pt_rap_phi_corr_%s", particleName[ipid].Data()) );
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_meas_thefirst[ipid] = (TH3D*) file1->Get( Form("hist_pt_rap_phi_meas_%s", particleName[ipid].Data()) );
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_rest[ipid] = (TH3D*) file1->Get( Form("hist_pt_rap_phi_rest_%s", particleName[ipid].Data()) );

  int nBinPt = hist_pt_rap_phi_real[0]->GetNbinsX();
  double binWidthPt = hist_pt_rap_phi_real[0]->GetXaxis()->GetBinWidth(1);
  double ptMin = hist_pt_rap_phi_real[0]->GetXaxis()->GetXmin();
  double ptMax = hist_pt_rap_phi_real[0]->GetXaxis()->GetXmax();

  int nBinRap = hist_pt_rap_phi_real[0]->GetNbinsY();
  double binWidthY = hist_pt_rap_phi_real[0]->GetYaxis()->GetBinWidth(1);
  double yyMin = hist_pt_rap_phi_real[0]->GetYaxis()->GetXmin();
  double yyMax = hist_pt_rap_phi_real[0]->GetYaxis()->GetXmax();

  int nBinPhi = hist_pt_rap_phi_real[0]->GetNbinsZ();
  double binWidthPhi = hist_pt_rap_phi_real[0]->GetZaxis()->GetBinWidth(1);
  double phiMin = hist_pt_rap_phi_real[0]->GetZaxis()->GetXmin();
  double phiMax = hist_pt_rap_phi_real[0]->GetZaxis()->GetXmax();
  int nBinPhi_half = ( nBinPhi - 1 )/2;


  double rapTest = -(0.2 + 0.001);
  nRapTest = hist_pt_rap_phi_rest[0]->GetYaxis()->FindBin( rapTest ) - 1;



  auto c = new TCanvas("c", "", 700, 600);


  auto g1 = new TGraph( hist_pt_rap_phi_real[nPidTest]->ProjectionZ("real_ang", nPtTest+1, nPtTest+1, nRapTest+1, nRapTest+1) );
  if( expdata ) 
  {
    auto hist_corr = (TH3D*) file_exp->Get( Form("hist_pt_rap_phi_corr_%s", particleName[nPidTest].Data()) );
    g1 = new TGraph( hist_corr->ProjectionZ("corr_ang", nPtTest+1, nPtTest+1, nRapTest+1, nRapTest+1) );
  }
  auto g2 = new TGraphErrors( hist_pt_rap_phi_meas_thefirst[nPidTest]->ProjectionZ("esti_ang", nPtTest+1, nPtTest+1, nRapTest+1, nRapTest+1) );
  auto g3 = new TGraphErrors( hist_pt_rap_phi_rest[nPidTest]->ProjectionZ("rest_ang", nPtTest+1, nPtTest+1, nRapTest+1, nRapTest+1) );
  cout << g1->GetMean(2) << " " << g2->GetMean(2) << " " << g3->GetMean(2) << endl;

  double g1Mean = 0;
  for( int i=0; i<g1->GetN()-1; i++ ) g1Mean += g1->GetY()[i]/(g1->GetN()-1);
  double g2Mean = 0;
  for( int i=0; i<g2->GetN()-1; i++ ) g2Mean += g2->GetY()[i]/(g2->GetN()-1);
  double g3Mean = 0;
  for( int i=0; i<g3->GetN()-1; i++ ) g3Mean += g3->GetY()[i]/(g3->GetN()-1);

  for( int i=0; i<g2->GetN(); i++ ) g2->SetPointError(i, 0, g2->GetEY()[i]);
  for( int i=0; i<g3->GetN(); i++ ) g3->SetPointError(i, 0, g3->GetEY()[i]);

  double val1 = g1Mean;
  double val3 = g3Mean;
  double val = val3/val1;
  for( int i=0; i<g1->GetN(); i++ ) g1->SetPoint( i, g1->GetX()[i], g1->GetY()[i]*val );

  double val2 = g2Mean;
  val = val3/val2;
  for( int i=0; i<g2->GetN(); i++ ) g2->SetPoint( i, g2->GetX()[i], g2->GetY()[i]*val );

  g1->SetMarkerColor( kRed );
  g2->SetMarkerColor( kGreen );
  g3->SetMarkerColor( kBlue );
  g1->SetMarkerStyle( 24 );
  g2->SetMarkerStyle( 26 );
  g3->SetMarkerStyle( 28 );

  double yMin = g1Mean * (1-0.40);
  double yMax = g1Mean * (1+0.40);
  g2->GetYaxis()->SetRangeUser( 0, g2->GetMean(2)*2. );
  //g1->GetYaxis()->SetRangeUser( yMin, yMax );
  g2->SetTitle( Form("ptcle = %s, Update = %d", particleName[nPidTest].Data(), numUp) );
  g2->GetXaxis()->SetTitle( "#phi (#phi_{ptcle} - #Psi_{RP})" );
  g2->GetYaxis()->SetTitle( "d^{3}N/p_{T}dp_{T}dyd#phi" );
  g2->Draw("AP");
  g1->Draw("Psame");
  g3->Draw("Psame");

  auto ll = new TLegend( 0.10, 0.15, 0.50, 0.3 );
  ll->SetBorderSize( 0 );
  ll->SetFillStyle( 0 );
  if( expdata ) ll->AddEntry( g1, "Corrected", "P" );
  else ll->AddEntry( g1, "Real", "P" );
  ll->AddEntry( g2, "Measured", "P" );
  ll->AddEntry( g3, "Restored", "P" );
  ll->Draw();
  double thisPt = hist_pt_rap_phi_rest[nPidTest]->GetXaxis()->GetBinCenter( nPtTest+1 );
  double thisRap = hist_pt_rap_phi_rest[nPidTest]->GetYaxis()->GetBinCenter( nRapTest+1 );
  ll->SetHeader( Form("p_{T} = %.2lf (GeV/c) /  y_{0} = %.2lf", thisPt, thisRap), "C" );
  ll->Draw();

  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_real[ipid]->SetName( Form("hist_pt_rap_phi_real_%s", particleName[ipid].Data()) );
  if( saveFigure ) c->SaveAs( Form("~/public_html/files/deblurring_spirit_eff/new_eff/update_%d.png", numUp) );




  if( drawFlow )
  {
    double y_cm = 0.372156;

    int nBinPt = hist_pt_rap_phi_real[nPidTest]->GetNbinsX();
    double binWidthPt = hist_pt_rap_phi_real[nPidTest]->GetXaxis()->GetBinWidth(1);
    double ptMin = hist_pt_rap_phi_real[nPidTest]->GetXaxis()->GetXmin();
    double ptMax = hist_pt_rap_phi_real[nPidTest]->GetXaxis()->GetXmax();

    int nBinRap = hist_pt_rap_phi_real[nPidTest]->GetNbinsY();
    double binWidthRap = hist_pt_rap_phi_real[nPidTest]->GetYaxis()->GetBinWidth(1);
    double yyMin = hist_pt_rap_phi_real[nPidTest]->GetYaxis()->GetXmin();
    double yyMax = hist_pt_rap_phi_real[nPidTest]->GetYaxis()->GetXmax();

    int nBinPhi = hist_pt_rap_phi_real[nPidTest]->GetNbinsZ();
    double binWidthPhy = hist_pt_rap_phi_real[nPidTest]->GetZaxis()->GetBinWidth(1);
    double phiMin = hist_pt_rap_phi_real[nPidTest]->GetZaxis()->GetXmin();
    double phiMax = hist_pt_rap_phi_real[nPidTest]->GetZaxis()->GetXmax();

    int momMin = 0;
    int momMax = nBinPt;

    auto cc = new TCanvas();

    auto gReal1 = new TGraphErrors();                                                                                                                                                                                                                 
    auto gReal2 = new TGraphErrors();
    auto gReal3 = new TGraphErrors();
    for( int irap=0; irap<nBinRap; irap++ )
    {
      double weight = 0;
      double thisV1 = 0;
      double thisV2 = 0;
      double thisV3 = 0;
      for( int imom=momMin; imom<momMax; imom++ )
      {
        for( int iphi=1; iphi<nBinPhi; iphi++ )
        {
          //if( hist_pt_rap_phi_real[nPidTest]->GetBinContent( imom+1, irap+1, iphi+1 ) < 0.15 ) continue;
          double val = hist_pt_rap_phi_real[nPidTest]->GetBinContent( imom+1, irap+1, iphi+1 ) * hist_pt_rap_phi_real[nPidTest]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
          double phi = hist_pt_rap_phi_real[nPidTest]->GetZaxis()->GetBinCenter( iphi+1 );
          weight += val;

          thisV1 += TMath::Cos(    phi )*val;
          thisV2 += TMath::Cos( 2.*phi )*val;
          thisV3 += TMath::Cos( 3.*phi )*val;
        }
      }
      thisV1 /= weight;
      thisV2 /= weight;
      thisV3 /= weight;

      double rap = hist_pt_rap_phi_real[nPidTest]->GetYaxis()->GetBinCenter( irap+1 );
      if( ystar ) rap *= y_cm;

      if( weight != 0 ) gReal1->SetPoint( gReal1->GetN(), rap, thisV1 );
      if( expdata ) gReal1->SetPoint( gReal1->GetN()-1, rap, thisV1/0.748 );
      if( weight != 0 ) gReal2->SetPoint( gReal2->GetN(), rap, thisV2 );
      if( expdata ) gReal2->SetPoint( gReal2->GetN()-1, rap, thisV2/0.416 );
      if( weight != 0 ) gReal3->SetPoint( gReal3->GetN(), rap, thisV3 );
    }



    nBinPt = hist_pt_rap_phi_rest[nPidTest]->GetNbinsX();
    binWidthPt = hist_pt_rap_phi_rest[nPidTest]->GetXaxis()->GetBinWidth(1);
    ptMin = hist_pt_rap_phi_rest[nPidTest]->GetXaxis()->GetXmin();
    ptMax = hist_pt_rap_phi_rest[nPidTest]->GetXaxis()->GetXmax();

    nBinRap = hist_pt_rap_phi_rest[nPidTest]->GetNbinsY();
    binWidthRap = hist_pt_rap_phi_rest[nPidTest]->GetYaxis()->GetBinWidth(1);
    yyMin = hist_pt_rap_phi_rest[nPidTest]->GetYaxis()->GetXmin();
    yyMax = hist_pt_rap_phi_rest[nPidTest]->GetYaxis()->GetXmax();

    nBinPhi = hist_pt_rap_phi_rest[nPidTest]->GetNbinsZ();
    binWidthPhy = hist_pt_rap_phi_rest[nPidTest]->GetZaxis()->GetBinWidth(1);
    phiMin = hist_pt_rap_phi_rest[nPidTest]->GetZaxis()->GetXmin();
    phiMax = hist_pt_rap_phi_rest[nPidTest]->GetZaxis()->GetXmax();

    momMin = 0;
    momMax = nBinPt;


    auto gMeas1 = new TGraphErrors();                                                                                                                                                                                                                 
    auto gMeas2 = new TGraphErrors();
    auto gMeas3 = new TGraphErrors();
    for( int irap=0; irap<nBinRap; irap++ )
    {
      double weight = 0;
      double thisV1 = 0;
      double thisV2 = 0;
      double thisV3 = 0;
      for( int imom=momMin; imom<momMax; imom++ )
      {
        for( int iphi=1; iphi<nBinPhi; iphi++ )
        {
          if( hist_pt_rap_phi_meas_thefirst[nPidTest]->GetBinContent( imom+1, irap+1, iphi+1 ) < 0.15 ) continue;
          double val = hist_pt_rap_phi_meas_thefirst[nPidTest]->GetBinContent( imom+1, irap+1, iphi+1 ) * hist_pt_rap_phi_meas_thefirst[nPidTest]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
          double phi = hist_pt_rap_phi_meas_thefirst[nPidTest]->GetZaxis()->GetBinCenter( iphi+1 );
          weight += val;

          thisV1 += TMath::Cos(    phi )*val;
          thisV2 += TMath::Cos( 2.*phi )*val;
          thisV3 += TMath::Cos( 3.*phi )*val;
        }
      }
      thisV1 /= weight;
      thisV2 /= weight;
      thisV3 /= weight;

      double rap = hist_pt_rap_phi_meas_thefirst[nPidTest]->GetYaxis()->GetBinCenter( irap+1 );
      if( ystar ) rap *= y_cm;

      if( weight != 0 ) gMeas1->SetPoint( gMeas1->GetN(), rap, thisV1 );
      if( weight != 0 ) gMeas2->SetPoint( gMeas2->GetN(), rap, thisV2 );
      if( weight != 0 ) gMeas3->SetPoint( gMeas3->GetN(), rap, thisV3 );
    }



    auto gRest1 = new TGraphErrors();                                                                                                                                                                                                                 
    auto gRest2 = new TGraphErrors();
    auto gRest3 = new TGraphErrors();
    for( int irap=0; irap<nBinRap; irap++ )
    {
      double weight = 0;
      double thisV1 = 0;
      double thisV2 = 0;
      double thisV3 = 0;
      for( int imom=momMin; imom<momMax; imom++ )
      {
        for( int iphi=1; iphi<nBinPhi; iphi++ )
        {
          if( hist_pt_rap_phi_rest[nPidTest]->GetBinContent( imom+1, irap+1, iphi+1 ) < 0.15 ) continue;
          double val = hist_pt_rap_phi_rest[nPidTest]->GetBinContent( imom+1, irap+1, iphi+1 ) * hist_pt_rap_phi_rest[nPidTest]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
          double phi = hist_pt_rap_phi_rest[nPidTest]->GetZaxis()->GetBinCenter( iphi+1 );
          weight += val;

          thisV1 += TMath::Cos(    phi )*val;
          thisV2 += TMath::Cos( 2.*phi )*val;
          thisV3 += TMath::Cos( 3.*phi )*val;
        }
      }
      thisV1 /= weight;
      thisV2 /= weight;
      thisV3 /= weight;

      double rap = hist_pt_rap_phi_rest[nPidTest]->GetYaxis()->GetBinCenter( irap+1 );
      if( ystar ) rap *= y_cm;

      if( weight != 0 ) gRest1->SetPoint( gRest1->GetN(), rap, thisV1 );
      if( weight != 0 ) gRest2->SetPoint( gRest2->GetN(), rap, thisV2 );
      if( weight != 0 ) gRest3->SetPoint( gRest3->GetN(), rap, thisV3 );
    }

    gReal1->GetXaxis()->SetTitle( "y_{0} = y_{lab}/y^{NN}_{cm} - 1" );

    if( ystar ) gReal1->GetXaxis()->SetTitle( "y^{*}" );
    if( ystar ) gMeas1->GetXaxis()->SetTitle( "y^{*}" );
    if( ystar ) gReal2->GetXaxis()->SetTitle( "y^{*}" );
    if( ystar ) gMeas2->GetXaxis()->SetTitle( "y^{*}" );

    gReal1->GetYaxis()->SetTitle( "v1" );
    gReal1->SetMarkerStyle( 25 );
    gReal2->SetMarkerStyle( 27 );
    gReal3->SetMarkerStyle( 29 );
    gReal1->SetMarkerSize( 2 );
    gReal2->SetMarkerSize( 2 );
    gReal3->SetMarkerSize( 2 );
    gReal1->SetMarkerColor( kRed  );
    gReal2->SetMarkerColor( kRed );
    gReal3->SetMarkerColor( kRed );
    gReal1->GetYaxis()->SetRangeUser( -0.7, 0.7 );

    gMeas1->SetMarkerStyle( 25 );
    gMeas2->SetMarkerStyle( 27 );
    gMeas3->SetMarkerStyle( 29 );
    gMeas1->SetMarkerSize( 2 );
    gMeas2->SetMarkerSize( 2 );
    gMeas3->SetMarkerSize( 2 );
    gMeas1->SetMarkerColor( kBlue  );
    gMeas2->SetMarkerColor( kBlue );
    gMeas3->SetMarkerColor( kBlue );

    gRest1->SetMarkerStyle( 25 );
    gRest2->SetMarkerStyle( 27 );
    gRest3->SetMarkerStyle( 29 );
    gRest1->SetMarkerSize( 2 );
    gRest2->SetMarkerSize( 2 );
    gRest3->SetMarkerSize( 2 );
    gRest1->SetMarkerColor( kGreen  );
    gRest2->SetMarkerColor( kGreen );
    gRest3->SetMarkerColor( kGreen );

    gReal1->SetLineColor( kRed );
    gMeas1->SetLineColor( kBlue );
    gRest1->SetLineColor( kGreen );
    gReal2->SetLineColor( kRed );
    gMeas2->SetLineColor( kBlue );
    gRest2->SetLineColor( kGreen );

    gReal1->SetLineStyle( 2 );
    gMeas1->SetLineStyle( 2 );
    gRest1->SetLineStyle( 2 );
    gReal2->SetLineStyle( 2 );
    gMeas2->SetLineStyle( 2 );
    gRest2->SetLineStyle( 2 );


    cc->SetGridx();
    cc->SetGridy();
    cc->cd();
    if( v_one )
    {
      gReal1->GetXaxis()->SetLimits( -1, 1 );
      gReal1->GetXaxis()->SetRangeUser( -1, 1 );
      gReal1->GetYaxis()->SetRangeUser( -0.60, 0.60 );
      gReal1->SetTitle( Form( "ptcle = %s / Update = %d", particleName[nPidTest].Data(), numUp ) );
      if( ystar ) gReal1->GetXaxis()->SetRangeUser( -0.2 , 0.4  );
      if( ystar ) gReal1->GetYaxis()->SetRangeUser( -0.25, 0.45 );
      gReal1->Draw("APL");
      gMeas1->Draw("PLsame");
      gRest1->Draw("PLsame");
    }
    else
    {
      //gReal2->GetXaxis()->SetLimits( -1, 1 );
      gReal2->GetXaxis()->SetRangeUser( -0.6, 1.2 );
      gReal2->GetYaxis()->SetRangeUser( -0.14, 0.14 );
      gReal2->SetTitle( Form( "ptcle = %s / Update = %d", particleName[nPidTest].Data(), numUp ) );
      if( ystar ) gReal2->GetXaxis()->SetRangeUser( -0.2 , 0.4  );
      if( ystar ) gReal2->GetYaxis()->SetRangeUser( -0.18, 0.28 );
      gReal2->Draw("APL");
      gMeas2->Draw("PLsame");
      gRest2->Draw("PLsame");

      if( ystar ) for( int i=0; i<gpa->GetN(); i++ ) gpa->SetPoint( i, gpa->GetX()[i]*y_cm, gpa->GetY()[i] );
      gpa->SetLineColor( kMagenta );
      gpa->SetLineStyle( 1 );
      gpa->SetLineWidth( 3 );
      gpa->Draw("Lsame");
    }

    auto l = new TLegend( 0.15, 0.7, 0.4, 0.9);                                                                                                                                                                                                 
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    if( expdata ) l->AddEntry( gReal1, "v1 corr", "P" );
    else l->AddEntry( gReal1, "v1 real", "P" );
    l->AddEntry( gMeas1, "v1 meas", "P" );
    l->AddEntry( gRest1, "v1 rest", "P" );
    //l->AddEntry( gReal2, "v2 real", "P" );
    l->Draw();
    if( saveFigure ) cc->SaveAs( Form("~/public_html/files/deblurring_spirit_eff/new_eff/flow.png") );
  }

}
