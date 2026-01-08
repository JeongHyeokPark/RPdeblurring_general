int isys = 1;

bool ystar = 1;
bool v_one = 0;
bool expdata = 1;

void drawOnePanel()
{
  int fileN = 50;

  TString sysName = isys ? "108" : "132";
  //auto file1 = new TFile( "./binning_even/output_err.root", "READ" );
  auto file1 = new TFile( "./input.root", "READ" );
  // ---------------------- Prepare 3D histograms ---------------------- 

  double y_cm = 0.372156;

  int nParticleType = 5;

  // No Real Distribution
  TString particleName[6] = { "proton", "deuteron", "triton", "3He", "alpha", "neutron" };
  std::vector <TH3D* > hist_pt_rap_phi_real; hist_pt_rap_phi_real.resize( 6 );
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_real[ipid] = (TH3D*) file1->Get( Form("hist_pt_rap_phi_real_%s", particleName[ipid].Data()) );
  if( expdata )
  {
    auto file_exp = new TFile( "/home/Symmetry/jhpark/spirit/ExpFlowData_Sn108_midcentral_opt.root", "READ" );
    for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_real[ipid] = (TH3D*) file_exp->Get( Form("hist_pt_rap_phi_corr_%s", particleName[ipid].Data()) );
  }
  std::vector <TH3D* > hist_pt_rap_phi_meas; hist_pt_rap_phi_meas.resize( 6 );
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_meas[ipid] = (TH3D*) file1->Get( Form("hist_pt_rap_phi_meas_%s", particleName[ipid].Data()) );
  std::vector <TH3D* > hist_pt_rap_phi_rest; hist_pt_rap_phi_rest.resize( 6 );
  for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_rest[ipid] = (TH3D*) file1->Get( Form("hist_pt_rap_phi_rest_%s", particleName[ipid].Data()) );
  //for( int ipid=0; ipid<nParticleType; ipid++ ) hist_pt_rap_phi_rest[ipid]->RebinY();

  int nBinPt = hist_pt_rap_phi_rest[0]->GetNbinsX();
  double binWidthPt = hist_pt_rap_phi_rest[0]->GetXaxis()->GetBinWidth(1);
  double ptMin = hist_pt_rap_phi_rest[0]->GetXaxis()->GetXmin();
  double ptMax = hist_pt_rap_phi_rest[0]->GetXaxis()->GetXmax();

  double momMin = 0;
  double momMax = nBinPt;

  int nBinRap = hist_pt_rap_phi_rest[0]->GetNbinsY();
  double binWidthRap = hist_pt_rap_phi_rest[0]->GetYaxis()->GetBinWidth(1);
  double yyMin = hist_pt_rap_phi_rest[0]->GetYaxis()->GetXmin();
  double yyMax = hist_pt_rap_phi_rest[0]->GetYaxis()->GetXmax();

  int nBinPhi = hist_pt_rap_phi_rest[0]->GetNbinsZ();
  double binWidthPhy = hist_pt_rap_phi_rest[0]->GetZaxis()->GetBinWidth(1);
  double phiMin = hist_pt_rap_phi_rest[0]->GetZaxis()->GetXmin();
  double phiMax = hist_pt_rap_phi_rest[0]->GetZaxis()->GetXmax();

  // ---------------------- Prepare 3D histograms ---------------------- 

  for( int nTestPid=0; nTestPid<5; nTestPid++ )
  {

    double cutRapidity = -99.;
    if( nTestPid==3 ) cutRapidity = -0.3;
    else cutRapidity = -0.5;

    auto c = new TCanvas( Form("c_%s", particleName[nTestPid].Data()), "", 700, 1000 );
    c->SetLeftMargin( 0.16 );
    c->SetRightMargin( 0.04 );
    c->SetTopMargin( 0.00 );
    c->SetBottomMargin( 0.16 );

    auto p1 = new TPad("p1", "", 0.00, 0.54, 1.0, 1.0);
    p1->SetBottomMargin(0);
    p1->SetTopMargin(0.04);
    p1->SetLeftMargin(0.16);
    p1->SetRightMargin(0.04);
    p1->Draw();
    auto p2 = new TPad("p2", "", 0.00, 0.00, 1.0, 0.54);
    p2->SetTopMargin(0);
    p2->SetLeftMargin(0.16);
    p2->SetBottomMargin(0.20);
    p2->SetRightMargin(0.04);
    p2->Draw();


    nBinPt = hist_pt_rap_phi_real[nTestPid]->GetNbinsX();
    binWidthPt = hist_pt_rap_phi_real[nTestPid]->GetXaxis()->GetBinWidth(1);
    ptMin = hist_pt_rap_phi_real[nTestPid]->GetXaxis()->GetXmin();
    ptMax = hist_pt_rap_phi_real[nTestPid]->GetXaxis()->GetXmax();
    momMin = 0;
    momMax = nBinPt;

    nBinRap = hist_pt_rap_phi_real[nTestPid]->GetNbinsY();
    binWidthRap = hist_pt_rap_phi_real[nTestPid]->GetYaxis()->GetBinWidth(1);
    yyMin = hist_pt_rap_phi_real[nTestPid]->GetYaxis()->GetXmin();
    yyMax = hist_pt_rap_phi_real[nTestPid]->GetYaxis()->GetXmax();

    nBinPhi = hist_pt_rap_phi_real[nTestPid]->GetNbinsZ();
    binWidthPhy = hist_pt_rap_phi_real[nTestPid]->GetZaxis()->GetBinWidth(1);
    phiMin = hist_pt_rap_phi_real[nTestPid]->GetZaxis()->GetXmin();
    phiMax = hist_pt_rap_phi_real[nTestPid]->GetZaxis()->GetXmax();

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
          //if( hist_pt_rap_phi_real[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 ) < 0.15 ) continue;
          double val = hist_pt_rap_phi_real[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 ) * hist_pt_rap_phi_real[nTestPid]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
          double phi = hist_pt_rap_phi_real[nTestPid]->GetZaxis()->GetBinCenter( iphi+1 );
          weight += val;

          thisV1 += TMath::Cos(    phi )*val;
          thisV2 += TMath::Cos( 2.*phi )*val;
          thisV3 += TMath::Cos( 3.*phi )*val;
        }
      }

      thisV1 /= weight;
      thisV2 /= weight;
      thisV3 /= weight;

      double rap = hist_pt_rap_phi_real[nTestPid]->GetYaxis()->GetBinCenter( irap+1 );
      if( rap < cutRapidity ) continue;
      if( ystar )rap *= y_cm;

      if( weight != 0 ) gReal1->SetPoint( gReal1->GetN(), rap, thisV1 );
      if( !isys&& expdata ) gReal1->SetPoint( gReal1->GetN()-1, rap, thisV1/0.747553 );
      if(  isys&& expdata ) gReal1->SetPoint( gReal1->GetN()-1, rap, thisV1/0.742754 );
      if( weight != 0 ) gReal2->SetPoint( gReal2->GetN(), rap, thisV2 );
      if( !isys&& expdata ) gReal2->SetPoint( gReal2->GetN()-1, rap, thisV2/0.416267 );
      if(  isys&& expdata ) gReal2->SetPoint( gReal2->GetN()-1, rap, thisV2/0.41033  );
      if( weight != 0 ) gReal3->SetPoint( gReal3->GetN(), rap, thisV3 );
    }







    nBinPt = hist_pt_rap_phi_rest[nTestPid]->GetNbinsX();
    binWidthPt = hist_pt_rap_phi_rest[nTestPid]->GetXaxis()->GetBinWidth(1);
    ptMin = hist_pt_rap_phi_rest[nTestPid]->GetXaxis()->GetXmin();
    ptMax = hist_pt_rap_phi_rest[nTestPid]->GetXaxis()->GetXmax();

  nBinRap = hist_pt_rap_phi_rest[nTestPid]->GetNbinsY();
  binWidthRap = hist_pt_rap_phi_rest[nTestPid]->GetYaxis()->GetBinWidth(1);
  yyMin = hist_pt_rap_phi_rest[nTestPid]->GetYaxis()->GetXmin();
  yyMax = hist_pt_rap_phi_rest[nTestPid]->GetYaxis()->GetXmax();

  nBinPhi = hist_pt_rap_phi_rest[nTestPid]->GetNbinsZ();
  binWidthPhy = hist_pt_rap_phi_rest[nTestPid]->GetZaxis()->GetBinWidth(1);
  phiMin = hist_pt_rap_phi_rest[nTestPid]->GetZaxis()->GetXmin();
  phiMax = hist_pt_rap_phi_rest[nTestPid]->GetZaxis()->GetXmax();

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
      double weightErr = 0;
      double thisErrV1 = 0;
      double thisErrV2 = 0;
      double thisErrV3 = 0;
      for( int imom=momMin; imom<momMax; imom++ )
      {
        for( int iphi=1; iphi<nBinPhi; iphi++ )
        {
          double phi = hist_pt_rap_phi_meas[nTestPid]->GetZaxis()->GetBinCenter( iphi+1 );

          if( hist_pt_rap_phi_meas[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 )>0.15 ) 
          {
            double val = hist_pt_rap_phi_meas[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 ) * hist_pt_rap_phi_meas[nTestPid]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
            weight += val;

            thisV1 += TMath::Cos(    phi )*val;
            thisV2 += TMath::Cos( 2.*phi )*val;
            thisV3 += TMath::Cos( 3.*phi )*val;
          }

          if( hist_pt_rap_phi_meas[nTestPid]->GetBinError( imom+1, irap+1, iphi+1 )>0.15 ) 
          {
            double err = hist_pt_rap_phi_meas[nTestPid]->GetBinError( imom+1, irap+1, iphi+1 ) * hist_pt_rap_phi_meas[nTestPid]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
            weightErr += err;

            thisErrV1 += TMath::Cos(    phi )*err;
            thisErrV2 += TMath::Cos( 2.*phi )*err;
            thisErrV3 += TMath::Cos( 3.*phi )*err;
          }
        }
      }
      thisV1 /= weight;
      thisV2 /= weight;
      thisV3 /= weight;

      thisErrV1 /= weightErr;
      thisErrV2 /= weightErr;
      thisErrV3 /= weightErr;

      double rap = hist_pt_rap_phi_meas[nTestPid]->GetYaxis()->GetBinCenter( irap+1 );
      if( rap < cutRapidity ) continue;
      if( ystar ) rap *= y_cm;

      if( weight != 0 ) gMeas1->SetPoint( gMeas1->GetN(), rap, thisV1 );
      if( weight != 0 ) gMeas2->SetPoint( gMeas2->GetN(), rap, thisV2 );
      if( weight != 0 ) gMeas3->SetPoint( gMeas3->GetN(), rap, thisV3 );
      //if( weightErr != 0 ) gMeas1->SetPointError( gMeas1->GetN()-1, 0., thisErrV1 );
      //if( weightErr != 0 ) gMeas2->SetPointError( gMeas2->GetN()-1, 0., thisErrV2 );
      //if( weightErr != 0 ) gMeas3->SetPointError( gMeas3->GetN()-1, 0., thisErrV3 );
    }



    bool sampleError = 1;
    auto gRest1 = new TGraphErrors();                                                                                                                                                                                                                 
    auto gRest2 = new TGraphErrors();
    auto gRest3 = new TGraphErrors();
    for( int irap=0; irap<nBinRap; irap++ )
    {
      double weight = 0;
      double thisV1 = 0;
      double thisV2 = 0;
      double thisV3 = 0;
      double weightErr = 0;
      double thisErrV1 = 0;
      double thisErrV2 = 0;
      double thisErrV3 = 0;
      for( int imom=momMin; imom<momMax; imom++ )
      {
        for( int iphi=1; iphi<nBinPhi; iphi++ )
        {
          double phi = hist_pt_rap_phi_rest[nTestPid]->GetZaxis()->GetBinCenter( iphi+1 );

          if( hist_pt_rap_phi_rest[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 )>0.15 ) 
          {
            double val = hist_pt_rap_phi_rest[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 ) * hist_pt_rap_phi_rest[nTestPid]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
            weight += val;

            thisV1 += TMath::Cos(    phi )*val;
            thisV2 += TMath::Cos( 2.*phi )*val;
            thisV3 += TMath::Cos( 3.*phi )*val;

            double err = hist_pt_rap_phi_rest[nTestPid]->GetBinError( imom+1, irap+1, iphi+1 ) * hist_pt_rap_phi_rest[nTestPid]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
            weightErr += err;

            thisErrV1 += err*TMath::Cos(    phi );
            thisErrV2 += err*TMath::Cos( 2.*phi );
            thisErrV3 += err*TMath::Cos( 3.*phi );
          }
        }
      }
      /*
      for( int imom=momMin; imom<momMax; imom++ )
      {
        for( int iphi=1; iphi<nBinPhi; iphi++ )
        {
          double phi = hist_pt_rap_phi_rest[nTestPid]->GetZaxis()->GetBinCenter( iphi+1 );

          if( hist_pt_rap_phi_rest[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 )>0.15 ) 
          {
            double err = hist_pt_rap_phi_rest[nTestPid]->GetBinError( imom+1, irap+1, iphi+1 ) * hist_pt_rap_phi_rest[nTestPid]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;

            thisErrV1 += err*err * TMath::Power((TMath::Cos(1.*phi)*weight-thisV1)/weight/weight, 2.);
            thisErrV2 += err*err * TMath::Power((TMath::Cos(2.*phi)*weight-thisV2)/weight/weight, 2.);
            thisErrV3 += err*err * TMath::Power((TMath::Cos(3.*phi)*weight-thisV3)/weight/weight, 2.);
          }
        }
      }
      thisErrV1 = TMath::Sqrt( thisErrV1 );
      thisErrV2 = TMath::Sqrt( thisErrV2 );
      thisErrV3 = TMath::Sqrt( thisErrV3 );
      */
      thisV1 /= weight;
      thisV2 /= weight;
      thisV3 /= weight;

      thisErrV1 /= weightErr; thisErrV2 /= weightErr; thisErrV3 /= weightErr;

      double rap = hist_pt_rap_phi_rest[nTestPid]->GetYaxis()->GetBinCenter( irap+1 );
      if( rap < cutRapidity ) continue;
      if( ystar ) rap *= y_cm;

      if( !sampleError&& weight != 0 ) gRest1->SetPoint( gRest1->GetN(), rap, thisV1 );
      if( !sampleError&& weight != 0 ) gRest2->SetPoint( gRest2->GetN(), rap, thisV2 );
      if( !sampleError&& weight != 0 ) gRest3->SetPoint( gRest3->GetN(), rap, thisV3 );
      if( !sampleError&& weight != 0 ) gRest1->SetPointError( gRest1->GetN()-1, 0., thisErrV1 );
      if( !sampleError&& weight != 0 ) gRest2->SetPointError( gRest2->GetN()-1, 0., thisErrV2 );
      if( !sampleError&& weight != 0 ) gRest3->SetPointError( gRest3->GetN()-1, 0., thisErrV3 );
    }



    for( int irap=0; sampleError&& irap<nBinRap; irap++ )
    {
      double rap = hist_pt_rap_phi_rest[nTestPid]->GetYaxis()->GetBinCenter( irap+1 );
      if( rap < cutRapidity ) continue;
      if( ystar ) rap *= y_cm;

      double weight = 0;
      double thisV1 = 0;
      double thisV2 = 0;
      double thisV3 = 0;
      double weightErr = 0;
      double thisErrV1 = 0;
      double thisErrV2 = 0;
      double thisErrV3 = 0;

      thisV1 = 0;
      thisV2 = 0;
      thisV3 = 0;
      for( int imom=momMin; imom<momMax; imom++ )
      {
        for( int iphi=1; iphi<nBinPhi; iphi++ )
        {
          double phi = hist_pt_rap_phi_rest[nTestPid]->GetZaxis()->GetBinCenter( iphi+1 );

          if( hist_pt_rap_phi_rest[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 )>0.15 ) 
          {
            double val = (hist_pt_rap_phi_rest[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 ))
            * hist_pt_rap_phi_rest[nTestPid]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
            weight += val;

            thisV1 += TMath::Cos(    phi )*val;
            thisV2 += TMath::Cos( 2.*phi )*val;
            thisV3 += TMath::Cos( 3.*phi )*val;
          }

        }
      }
      thisV1 /= weight;
      thisV2 /= weight;
      thisV3 /= weight;

      for( int isamp=0; isamp<50; isamp++ )
      {
        double fweight = 0.;
        double fthisV1 = 0;
        double fthisV2 = 0;
        double fthisV3 = 0;
        for( int imom=momMin; imom<momMax; imom++ )
        {
          for( int iphi=1; iphi<nBinPhi; iphi++ )
          {
            double phi = hist_pt_rap_phi_rest[nTestPid]->GetZaxis()->GetBinCenter( iphi+1 );

            if( hist_pt_rap_phi_rest[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 )>0.15 ) 
            {
              double val = (hist_pt_rap_phi_rest[nTestPid]->GetBinContent( imom+1, irap+1, iphi+1 ) + hist_pt_rap_phi_rest[nTestPid]->GetBinError( imom+1, irap+1, iphi+1 )*gRandom->Gaus(0,1))
                * hist_pt_rap_phi_rest[nTestPid]->GetXaxis()->GetBinCenter( imom+1 ) * binWidthPt * binWidthRap * binWidthPhy;
              if( hist_pt_rap_phi_rest[nTestPid]->GetBinError( imom+1, irap+1, iphi+1 )==0 ) val = 0;
              fweight += val;

              fthisV1 += TMath::Cos(    phi )*val;
              fthisV2 += TMath::Cos( 2.*phi )*val;
              fthisV3 += TMath::Cos( 3.*phi )*val;
            }

          }
        }
        fthisV1 /= fweight;
        fthisV2 /= fweight;
        fthisV3 /= fweight;

        thisErrV1 += TMath::Power(thisV1-fthisV1, 2.);
        thisErrV2 += TMath::Power(thisV2-fthisV2, 2.);
        thisErrV3 += TMath::Power(thisV3-fthisV3, 2.);
      }
      thisErrV1 = TMath::Sqrt( thisErrV1/50. );
      thisErrV2 = TMath::Sqrt( thisErrV2/50. );
      thisErrV3 = TMath::Sqrt( thisErrV3/50. );

      if( weight != 0 ) gRest1->SetPoint( gRest1->GetN(), rap, thisV1 );
      if( weight != 0 ) gRest2->SetPoint( gRest2->GetN(), rap, thisV2 );
      if( weight != 0 ) gRest3->SetPoint( gRest3->GetN(), rap, thisV3 );
      if( weight != 0 ) gRest1->SetPointError( gRest1->GetN()-1, 0., thisErrV1 );
      if( weight != 0 ) gRest2->SetPointError( gRest2->GetN()-1, 0., thisErrV2 );
      if( weight != 0 ) gRest3->SetPointError( gRest3->GetN()-1, 0., thisErrV3 );
    }




    if( ystar ) gReal1->GetXaxis()->SetTitle( "y^{*}" );
    if( ystar ) gMeas1->GetXaxis()->SetTitle( "y^{*}" );
    if( ystar ) gReal2->GetXaxis()->SetTitle( "y^{*}" );
    if( ystar ) gMeas2->GetXaxis()->SetTitle( "y^{*}" );







    gRest1->SetMarkerStyle( 34 );
    gRest2->SetMarkerStyle( 34 );
    gRest3->SetMarkerStyle( 34 );
    gRest1->SetMarkerSize( 3 );
    gRest2->SetMarkerSize( 3 );
    gRest3->SetMarkerSize( 3 );
    gRest1->SetMarkerColor( kRed  );
    gRest2->SetMarkerColor( kRed );
    gRest3->SetMarkerColor( kRed );
    gRest1->GetYaxis()->SetRangeUser( -0.7, 0.7 );

    gMeas1->SetMarkerStyle( 56 );
    gMeas2->SetMarkerStyle( 56 );
    gMeas3->SetMarkerStyle( 56 );
    gMeas1->SetMarkerSize( 3 );
    gMeas2->SetMarkerSize( 3 );
    gMeas3->SetMarkerSize( 3 );
    gMeas1->SetMarkerColor( kAzure+1 );
    gMeas2->SetMarkerColor( kAzure+1 );
    gMeas3->SetMarkerColor( kAzure+1 );

    gReal1->SetMarkerStyle( 74 );
    gReal2->SetMarkerStyle( 74 );
    gReal3->SetMarkerStyle( 74 );
    //gReal1->SetMarkerStyle( 33 );
    //gReal2->SetMarkerStyle( 33 );
    //gReal3->SetMarkerStyle( 33 );
    gReal1->SetMarkerSize( 4 );
    gReal2->SetMarkerSize( 4 );
    gReal3->SetMarkerSize( 3 );
    gReal1->SetMarkerColor( kGreen+3  );
    gReal2->SetMarkerColor( kGreen+3 );
    gReal3->SetMarkerColor( kGreen+3 );
    gReal1->SetMarkerColorAlpha( kGreen+3, 0.60 );
    gReal2->SetMarkerColorAlpha( kGreen+3, 0.60 );
    gReal3->SetMarkerColorAlpha( kGreen+3, 0.60 );

    gReal1->SetLineColor( kGreen+3 );
    gMeas1->SetLineColor( kBlue );
    gRest1->SetLineColor( kRed );
    gReal2->SetLineColor( kGreen+3 );
    gMeas2->SetLineColor( kBlue );
    gRest2->SetLineColor( kRed );

    gReal1->SetLineStyle( 2 );
    gMeas1->SetLineStyle( 2 );
    gRest1->SetLineStyle( 2 );
    gReal2->SetLineStyle( 2 );
    gMeas2->SetLineStyle( 2 );
    gRest2->SetLineStyle( 2 );


    TString collSysName = "";
    if( isys ) collSysName = "^{108}Sn+^{112}Sn @ 270 AMeV";
    else collSysName = "^{132}Sn+^{124}Sn @ 270 AMeV";

    TString impSys = "1.6 fm < b < 3.9 fm";



    if( ystar ) gMeas1->GetXaxis()->SetRangeUser( -0.25, 0.45 );
    p1->cd();

    auto g = new TMultiGraph();


    g->GetXaxis()->SetTitle( "y_{0} = y_{lab}/y^{NN}_{cm} - 1" );
    g->GetYaxis()->SetTitle( "v_{1}" );

    if( ystar ) g->GetXaxis()->SetTitle( "y^{*}" );


    g->GetYaxis()->CenterTitle();
    g->GetYaxis()->SetNdivisions(505);
    g->GetYaxis()->SetTitleSize(38);
    g->GetYaxis()->SetTitleFont(45);
    g->GetYaxis()->SetTitleOffset(1.2);
    g->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g->GetYaxis()->SetLabelSize(30);

    // X axis ratio plot settings
    g->GetXaxis()->SetTitleSize(37);
    g->GetXaxis()->SetTitleFont(43);
    g->GetXaxis()->SetTitleOffset(4);
    g->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g->GetXaxis()->SetLabelSize(30);


    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();

    //g->Add( gMeas1 );
    g->Add( gRest1 );
    g->Add( gReal1 );
    g->GetXaxis()->SetLimits( -0.25/y_cm, 0.45/y_cm );
    g->GetXaxis()->SetRangeUser( -0.25/y_cm, 0.45/y_cm );
    g->GetYaxis()->SetRangeUser( -0.65, 0.65 );
    if( ystar ) g->GetXaxis()->SetLimits( -0.25, 0.45 );
    if( ystar ) g->GetXaxis()->SetRangeUser( -0.25, 0.45 );
    g->Draw("APLE");



    auto gl = new TGraph();
    gl->SetPoint( gl->GetN(), -100, 0 );
    gl->SetPoint( gl->GetN(),  100, 0 );
    gl->SetLineStyle( 2 );
    gl->Draw( "Lsame" );


    p2->cd();
    gReal2->GetXaxis()->SetLimits( -1, 2 );
    gReal2->GetXaxis()->SetRangeUser( -1, 1.2 );
    gReal2->GetYaxis()->SetRangeUser( -0.15, 0.35 );

    if( ystar ) gReal2->GetXaxis()->SetRangeUser( -0.25, 0.45 );

    auto g2 = new TMultiGraph();
    g2->GetXaxis()->SetTitle( "y_{0} = y_{lab}/y^{NN}_{cm} - 1" );
    g2->GetYaxis()->SetTitle( "v_{2}" );
    if( ystar ) g2->GetXaxis()->SetTitle( "y^{*} = y_{lab} - y^{NN}_{cm}" );
    if( ystar ) g2->GetXaxis()->SetTitle( "y^{*}" );


    g2->GetXaxis()->CenterTitle();
    g2->GetYaxis()->CenterTitle();

    //g2->Add( gMeas2 );
    g2->Add( gRest2 );
    g2->Add( gReal2 );

    g2->GetXaxis()->SetLimits( -1, 1.2 );
    g2->GetXaxis()->SetRangeUser( -1, 1.2 );
    g2->GetYaxis()->SetRangeUser( -0.16, 0.16 );

    g2->GetXaxis()->SetLimits( -0.25/y_cm, 0.45/y_cm );
    g2->GetXaxis()->SetRangeUser( -0.25/y_cm, 0.45/y_cm );
    if( ystar ) g2->GetXaxis()->SetLimits( -0.25, 0.45 );
    if( ystar ) g2->GetXaxis()->SetRangeUser( -0.25, 0.45 );


    g2->GetYaxis()->CenterTitle();
    //g2->GetYaxis()->SetNdivisions(505);
    g2->GetYaxis()->SetTitleSize(38);
    g2->GetYaxis()->SetTitleFont(45);
    g2->GetYaxis()->SetTitleOffset(1.2);
    g2->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g2->GetYaxis()->SetLabelSize(30);

    // X axis ratio plot settings
    g2->GetXaxis()->SetTitleSize(38);
    g2->GetXaxis()->SetTitleFont(43);
    g2->GetXaxis()->SetTitleOffset(1.2);
    g2->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    g2->GetXaxis()->SetLabelSize(30);

    g2->Draw("APLE");




    auto l = new TLegend( 0.17, 0.65, 0.5, 0.95);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextSize( 0.065 );

    //l->AddEntry( gMeas3, "Estimated", "P" );
    //if( expdata ) l->AddEntry( gReal1, "Corrected", "P" );
    //l->AddEntry( gMeas1, "w/o corr", "P" );
    //if( expdata ) l->AddEntry( gReal1, "w/   corr", "P" );
       //l->AddEntry( gRest1, "<cos n#phi>", "P" );
    //l->AddEntry( gMeas1, "#LT cos #it{n#phi'} #GT ", "P" );
    //l->AddEntry( gReal1, "#LT cos #it{n#phi}  #GT ", "P" );

    l->AddEntry( gReal3, "Mizuki", "P" );
    l->AddEntry( gRest1, "Deblurred", "P" );

    p1->cd();

    auto ll = new TLegend( 0.50, 0.25, 1.00, 0.35 );
    //ll->SetHeader( Form("%s", particleName[nTestPid].Data()) );
    ll->SetHeader( collSysName );
    ll->SetFillStyle(0);
    ll->SetBorderSize(0);
    ll->SetTextSize( 0.06 );
    ll->Draw("same");
    auto lll = new TLegend( 0.50, 0.20, 1.00, 0.25 );
    lll->SetHeader( impSys );
    lll->SetFillStyle(0);
    lll->SetBorderSize(0);
    lll->SetTextSize( 0.06 );
    lll->Draw("same");
    auto llll = new TLegend( 0.50, 0.13, 1.00, 0.20 );
    llll->SetHeader( particleName[nTestPid] );
    llll->SetFillStyle(0);
    llll->SetBorderSize(0);
    llll->SetTextSize( 0.06 );
    llll->Draw("same");


    p1->cd();
    l->Draw();

    TString figName = "plots/";
    figName += particleName[nTestPid] + "_";
    figName += "all.pdf";

    c->SaveAs( figName );
  }

}
