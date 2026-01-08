void test()
{
  auto file1 = new TFile( "binning_even/input.root", "READ" );
  auto file2 = new TFile( "binning_uneven/input.root", "READ" );

  auto hist1 = (TH3D*) file1->Get( "hist_pt_rap_phi_rest_alpha" );
  auto hist2 = (TH3D*) file2->Get( "hist_pt_rap_phi_rest_alpha" );

  int nBinPt = hist1->GetNbinsX();
  double binWidthPt = hist1->GetXaxis()->GetBinWidth(1);
  double ptMin = hist1->GetXaxis()->GetXmin();
  double ptMax = hist1->GetXaxis()->GetXmax();

  int nBinRap = hist1->GetNbinsY();
  double binWidthRap = hist1->GetYaxis()->GetBinWidth(1);
  double yyMin = hist1->GetYaxis()->GetXmin();
  double yyMax = hist1->GetYaxis()->GetXmax();

  int nBinPhi = hist1->GetNbinsZ();
  double binWidthPhy = hist1->GetZaxis()->GetBinWidth(1);
  double phiMin = hist1->GetZaxis()->GetXmin();
  double phiMax = hist1->GetZaxis()->GetXmax();

    int momMin = 0;
    int momMax = nBinPt;

  auto gRest1 = new TGraphErrors();                                                                                                                                                                                                                 
  for( int irap=0; irap<nBinRap; irap++ )
  {
    double weight = 0;
    double thisV2 = 0;
    for( int imom=momMin; imom<momMax; imom++ )
    {
      for( int iphi=1; iphi<nBinPhi; iphi++ )
      {
        if( hist1->GetBinContent( imom+1, irap+1, iphi+1 ) < 0.15 ) continue;
        double val = hist1->GetBinContent( imom+1, irap+1, iphi+1 ) * hist1->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
        double phi = hist1->GetZaxis()->GetBinCenter( iphi+1 );
        weight += val;

        thisV2 += TMath::Cos(  2*phi )*val;
      }
    }
    thisV2 /= weight;

    double rap = hist1->GetYaxis()->GetBinCenter( irap+1 );
    //if( ystar ) rap *= y_cm;

    if( weight != 0 ) gRest1->SetPoint( gRest1->GetN(), rap, thisV2 );
  }


  nBinPt = hist2->GetNbinsX();
  binWidthPt = hist2->GetXaxis()->GetBinWidth(1);
  ptMin = hist2->GetXaxis()->GetXmin();
  ptMax = hist2->GetXaxis()->GetXmax();

  nBinRap = hist2->GetNbinsY();
  binWidthRap = hist2->GetYaxis()->GetBinWidth(1);
  yyMin = hist2->GetYaxis()->GetXmin();
  yyMax = hist2->GetYaxis()->GetXmax();

  nBinPhi = hist2->GetNbinsZ();
  binWidthPhy = hist2->GetZaxis()->GetBinWidth(1);
  phiMin = hist2->GetZaxis()->GetXmin();
  phiMax = hist2->GetZaxis()->GetXmax();

  auto gRest2 = new TGraphErrors();                                                                                                                                                                                                                 
  for( int irap=0; irap<nBinRap; irap++ )
  {
    double weight = 0;
    double thisV2 = 0;
    for( int imom=momMin; imom<momMax; imom++ )
    {
      for( int iphi=1; iphi<nBinPhi; iphi++ )
      {
        if( hist2->GetBinContent( imom+1, irap+1, iphi+1 ) < 0.15 ) continue;
        double val = hist2->GetBinContent( imom+1, irap+1, iphi+1 ) * hist2->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
        double phi = hist2->GetZaxis()->GetBinCenter( iphi+1 );
        weight += val;

        thisV2 += TMath::Cos(  2*phi )*val;
      }
    }
    thisV2 /= weight;

    double rap = hist2->GetYaxis()->GetBinCenter( irap+1 );
    //if( ystar ) rap *= y_cm;

    if( weight != 0 ) gRest2->SetPoint( gRest2->GetN(), rap, thisV2 );
  }


  gRest1->SetMarkerStyle( 25 );
  gRest2->SetMarkerStyle( 27 );
  gRest1->SetMarkerSize( 2 );
  gRest2->SetMarkerSize( 2 );
  gRest1->SetMarkerColor( kGreen  );
  gRest2->SetMarkerColor( kGreen );

  gRest1->SetLineColor( kGreen );
  gRest2->SetLineColor( kGreen );

  gRest1->SetLineStyle( 2 );
  gRest2->SetLineStyle( 2 );


  gRest1->GetXaxis()->SetRangeUser( -0.6, 1.2 );
  gRest1->GetYaxis()->SetRangeUser( -0.14, 0.14 );
  //gRest1->SetTitle( Form( "ptcle = %s / Update = %d", particleName[nPidTest].Data(), numUp ) );
  gRest1->Draw("APL");
  gRest2->Draw("PLsame");
}
