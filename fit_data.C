// ---------------------------------------------------------------------------
// General parameters, we only consider alphas for the test

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
double mass = mass_a4;

double dn_cut = 0.15;

// General parameters, we only consider alphas for the test
// ---------------------------------------------------------------------------

//TEfficiency* teff = nullptr;
double rap = -99;

TH2D* teff = nullptr;
double ThermalDist( double pt, double phi, double d2Nmax, double vx, double temp, double t2 );
double modelFunction( double *x, double *par );


void fit_data()
{
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Simplex");

  //auto file_eff = new TFile( "eff_Sn108_midcentral.root", "READ" );
  //teff = (TEfficiency*) file_eff->Get( "e3_alpha_HighRes" );
  auto file_eff = new TFile( "data_midcentral.root", "READ" );
  teff = (TH2D*) file_eff->Get( "h2PtYEff_108Sn_4He_mbin0_iter0" );

  auto file = new TFile( "output_err.root", "READ" );
  //auto file = new TFile( "input.root", "READ" );
  auto hist = (TH3D*) file->Get( "hist_pt_rap_phi_rest_alpha" );

  int nBinPt = hist->GetNbinsX();
  double binWidthPt = hist->GetXaxis()->GetBinWidth(1);
  double ptMin = hist->GetXaxis()->GetXmin();
  double ptMax = hist->GetXaxis()->GetXmax();

  int nBinRap = hist->GetNbinsY();
  double binWidthRap = hist->GetYaxis()->GetBinWidth(1);
  double yyMin = hist->GetYaxis()->GetXmin();
  double yyMax = hist->GetYaxis()->GetXmax();

  int nBinPhi = hist->GetNbinsZ();
  double binWidthPhy = hist->GetZaxis()->GetBinWidth(1);
  double phiMin = hist->GetZaxis()->GetXmin();
  double phiMax = hist->GetZaxis()->GetXmax();

  auto fModel = new TF2( "fModel", modelFunction, ptMin, ptMax, phiMin, phiMax, 4 );
  auto hFit = new TH2D( "hFit", "", nBinPt, ptMin, ptMax, nBinPhi, phiMin, phiMax );

  auto hh = new TH2D( "hh", "", nBinRap, yyMin, yyMax, 100, -2, 2 );

  auto gtest = new TGraph();

  // Prepare parameters, (ANORM-Max. d3N{ipid,y}, vx-px/m, 
  // Temp-{TempIn-px2/m, TempOut-py2/m}, t2)

  //For every rapidity bins, (considering only one particle)
  for( int iy=0; iy<nBinRap; iy++ )
  {
    hFit->Reset();

    rap = hist->GetYaxis()->GetBinCenter( iy+1 );

    double weight  = 0.;
    double avg_px  = 0.;
    double avg_px2 = 0.;
    double avg_py  = 0.;
    double avg_py2 = 0.;
    double avg_v1  = 0.;
    double avg_v2  = 0.;

    double d2Nmax = 0.; // This is the first parameter, ANORM
    for( int ip=0; ip<nBinPt; ip++ )
    {
      double pt = hist->GetXaxis()->GetBinCenter( ip+1 );
      for( int iphi=0; iphi<nBinPhi-1; iphi++ ) // Note the both end bins are compensated each other
      {
        double phi = hist->GetZaxis()->GetBinCenter( iphi+1 );

        double content = hist->GetBinContent( ip+1, iy+1, iphi+1 );
        if( content < dn_cut ) 
        {
          hist->SetBinContent( ip+1, iy+1, iphi+1, 0. );
          hist->SetBinError( ip+1, iy+1, iphi+1, 0. );
          continue;
        }
        if( d2Nmax < content ) d2Nmax = content;

        // Normalize back to count
        content = hist->GetBinContent( ip+1, iy+1, iphi+1 )*pt*binWidthPt*binWidthRap*binWidthPhy ;

        double px  = pt*TMath::Cos( phi );
        double py  = pt*TMath::Sin( phi );
        double px2 = px*px;
        double py2 = py*py;

        weight  += content;
        avg_px  += content * px; avg_px2 += content * px2;
        avg_py  += content * py; avg_py2 += content * py2;
        avg_v1  += content * TMath::Cos( 1.*phi );
        avg_v2  += content * TMath::Cos( 2.*phi );
      }
    }

    if( weight ) // To prevent segmentation due to dividing by zero
    {
      avg_px  /= weight; avg_px2 /= weight;
      avg_py  /= weight; avg_py2 /= weight;
      avg_v1  /= weight; avg_v2  /= weight;
    }

    double vx_src  = avg_px / mass; // <px>/m
    double tempIn  = ( avg_px2 - avg_px*avg_px ) / mass; // ( <px^2> - <px>^2 )/m
    double tempOut = ( avg_py2                 ) / mass; // ( <py^2> )/m
    double temp    = 2. / ( 1./tempIn + 1./tempOut );
    //double t2      = ( tempOut - tempIn )/( tempOut + tempIn ); // Alternative: temp * ( 1./tempIn - 1./tempOut )/2.;
    double t2      = temp * ( 1./tempIn - 1./tempOut )/2.;




    // Now the parameters are: d2Nmax, vx_src, temp, t2 for fitting d2N(pt, phi)

    for( int ip=0; ip<nBinPt; ip++ )
    {
      for( int iphi=0; iphi<nBinPhi; iphi++ )
      {
        double val = hist->GetBinContent( ip+1, iy+1, iphi+1 );
        double err = hist->GetBinError( ip+1, iy+1, iphi+1 );
        if( dn_cut < val ) 
        {
          hFit->SetBinContent( ip+1, iphi+1, val );
          hFit->SetBinError( ip+1, iphi+1, err );
        }
        else 
        {
          hFit->SetBinContent( ip+1, iphi+1, 0 );
          hFit->SetBinError( ip+1, iphi+1, 0 );
        }

      }
    }

    if( hFit->GetEntries()==0 ) continue;



    cout << "RAP: " << rap << endl;
    cout << "PARS: " << d2Nmax << " " << vx_src << " " << temp << " " << t2 << endl;
    fModel->SetParameters( d2Nmax, vx_src, temp, t2 ); // Rapidity par for efficiency calc
    //hFit->Fit( fModel, "QN", "" );
    hFit->Fit( fModel, "N", "" );
    cout << endl;

    //if( iy==12 )
    {
      auto ge = new TGraphErrors();
      auto gf = new TGraph();
      auto gf1 = new TGraph();

      for( int ip=0; ip<nBinPt; ip++ )
      {
        double pt = hist->GetXaxis()->GetBinCenter( ip+1 );
        double p1 = hist->GetBinContent( ip+1, iy+1, hist->GetZaxis()->FindBin(0.) );
        double p2 = hist->GetBinContent( ip+1, iy+1, hist->GetZaxis()->FindBin(-TMath::Pi()) );
        double p1Err = hist->GetBinError( ip+1, iy+1, hist->GetZaxis()->FindBin(0.) );
        double p2Err = hist->GetBinError( ip+1, iy+1, hist->GetZaxis()->FindBin(-TMath::Pi()) );

        ge->SetPoint( ge->GetN(),  pt/mass, p1 );
        ge->SetPointError( ge->GetN()-1, 0, p1Err );
        ge->SetPoint( ge->GetN(), -pt/mass, p2 );
        ge->SetPointError( ge->GetN()-1, 0, p2Err );
      }

        for( int i=0; i<40; i++ )
        {
          double px  = -2.0 + 0.1*i;
          double pt  = std::abs(px);
          double phi = px>0. ? 0 : TMath::Pi();
          double p1 = fModel->Eval( pt, phi );
          //p1 = ThermalDist( pt, phi, d2Nmax, vx_src, temp, t2 );

          gf->SetPoint( gf->GetN(), px/mass, p1 );
          //gf1->SetPoint( gf->GetN(), px/mass, 
          //ThermalDist(pt, phi, fModel->GetParameter(1), fModel->GetParameter(2), fModel->GetParameter(3), fModel->GetParameter(4)) );
        }

      ge->SetMarkerStyle( 24 );
      ge->SetMarkerColor( kRed );
      gf->SetLineColor( kRed );
      //gf1->SetLineColor( kBlue );

      auto c = new TCanvas();
      c->SetLogy();
      ge->SetTitle( Form("rap = %.4lf", rap) );
      ge->Draw( "APE" );
      gf->Draw( "Lsame" );
      //gf1->Draw( "Lsame" );
    }

  }
}


double modelFunction( double *x, double *par )
{
  double pt  = x[0];
  double phi = x[1];

  //double rap    = par[0];
  double d2Nmax = par[0];
  double vx     = par[1];
  double temp   = std::abs(par[2]);
  double t2     = par[3];

  double dN = ThermalDist( pt, phi, d2Nmax, vx, temp, t2 );

  //int thisBin = teff->FindFixBin( pt*1e+3, rap, phi*TMath::RadToDeg() );
  //double eff = teff->GetEfficiency( thisBin ); // Need to be implemented
  //int thisBin = teff->FindFixBin( rap, pt*1e+3 );
  //double eff = teff->GetBinContent( thisBin ); // Need to be implemented

  double eff = 1.;

  double result = 0.;
  if( 0. < dN ) result = dN*eff;

  return result;
}


double ThermalDist( double pt, double phi, double d2Nmax, double vx, double temp, double t2 )
{
  double gam  = 1./( 1. + vx*vx );
  double gamb = vx*gam;

  double pt2 = pt*pt;
  double et = TMath::Sqrt( mass*mass + pt2 );

  double px = pt * TMath::Cos( phi );
  double py = pt * TMath::Sin( phi );
  double py2 = py*py;

  double px_loc = gam*px - gamb*et;
  double px2_loc = px_loc*px_loc;
  double pt2_loc = px2_loc + py2;
  double et_loc = TMath::Sqrt( mass*mass + pt2_loc );

  double phi_loc = 0;
  if( py!= 0. || px_loc!= 0. )
    phi_loc = TMath::ATan2( py, px_loc );

  double den = TMath::Max( std::abs(temp), 1e-6 );
  double pden = (1.+t2*TMath::Cos(2*phi_loc)) / den;
  double ag = pt2_loc*pden / (et_loc + mass);

  double res = 0;
  if( std::abs(ag) < 60. )
    res = d2Nmax*TMath::Exp(-ag) * et_loc / et;

  return res;
}
