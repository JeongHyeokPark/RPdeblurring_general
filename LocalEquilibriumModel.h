#include <TROOT.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <stdio.h>

class LocalEquilibriumModel {

  public: 
    void PrepareParameters( double ke, double avg_charge );
    void BESKN( double xx, int NU, int KODE, int N, double &yy, int &NZ );
    double BESI01( double xx, int NU, int KODE );
    double BESK01( double xx, int NU, int KODE, int &MZ ); // No MZ USED
    void goto96( int &NZ, int N, int ND, std::vector< double > &Y );
    void POISSON( double AV, int &M );
    void GetMomentum( double &PX, double &PY, double &PZ, double &P, double AM, double T );
    void LOREN( double BGX, double BGY, double BGZ, double G, double PCX, double PCY, double PCZ, double EC, double &PX, double &PY, double &PZ, double &E );
    void EmitParticlesInCM( std::vector<int> *fZ, std::vector<int> *fA, std::vector<double> *fMass, std::vector<TLorentzVector> *fVector );
    void CalcAvgMulti();
    double GetNuclearMass( int Z, int A );

  private:
    double ke_beam = -1;
    double avg_multi = -1;

    double avg_multi_p; double avg_multi_d; double avg_multi_t;
    double avg_multi_h3; double avg_multi_h4;

    int multi_max; int multi_max2;

    const int nParticleType = 5;
    double particle_mass[6];
    int particle_a[6]       = { 1, 2, 3, 3, 4, 1 };
    int particle_z[6]       = { 1, 1, 1, 2, 2, 0 };
    double mass_proton; double mass_neutron; double mass_avg;
    double mass_minus_be; double mass_a2;
    double be;

    double mass_avg2, be_d;
    double mass_avg3, be_t, mass_minus_be_t, be_3he, mass_minus_be_3he, mass_a3;
    double mass_avg4, be_4he, mass_a4;
};

void LocalEquilibriumModel::PrepareParameters( double ke, double avg_charge )
{
  ke_beam = ke;
  avg_multi = avg_charge;
  multi_max = 200;
  multi_max2 = 2*multi_max;

  // ---------------- contents in flowsym_mod.for -------------------------
  //  PARTICLE DATA
  //  NUCLEONS
  be = .008;
  mass_proton     = .9383; 
  mass_neutron    = .9396; 
  mass_avg        = .5*(mass_proton+mass_neutron);
  mass_minus_be   = mass_avg - be;
  //mass_avg_square = mass_avg*mass_avg;

  //  DEUTERON
  mass_avg2       = mass_avg + mass_avg;
  be_d            = .002225; 
  mass_a2 = mass_avg + mass_avg - be_d;

  //  A=3
  mass_avg3         = mass_avg2 + mass_avg;
  be_t              = .0086; 
  be_3he            = .0080; 
  mass_minus_be_t   = mass_proton + mass_neutron + mass_neutron - be_t;
  mass_minus_be_3he = mass_proton + mass_proton + mass_neutron - be_3he;
  mass_a3           = .5*(mass_minus_be_t+mass_minus_be_3he);

  //  ALPHA
  mass_avg4 = mass_avg3 + mass_avg;
  be_4he    = .0286; 
  mass_a4   = mass_avg4 - be_4he;

  particle_mass[0] = mass_avg;
  particle_mass[1] = mass_a2;
  particle_mass[2] = mass_a3;
  particle_mass[3] = mass_a3;
  particle_mass[4] = mass_a4;
  particle_mass[5] = mass_avg;
  // ---------------- contents in flowsym_mod.for -------------------------
}

void LocalEquilibriumModel::BESKN( double xx, int NU, int KODE, int N, double &yy, int &NZ )
{
  int NULIM[2] = { 35 , 70 };
  double ELIM = 667.0;
  NZ++;

  std::vector< double > Y;
  Y.resize( NZ );
  for( int i=0; i<NZ; i++ ) Y[i] = 0;
  Y[0]=yy;

  // ---------- ERRCHK --------------
  // TEST INPUT ARGUMENTS
  if(KODE<1 || KODE>2) return; // goto 90
  if(NU<0) return; // goto 91
  if(xx<=0.) return; // goto 92
  if(N<1) return; // goto 93
  // ---------- ERRCHK --------------

  int ETX = KODE-1;
  // NUD AND ND ARE DUMMY VALUES FOR NU AND N
  // NZ = NUMBER OF UNDERFLOWS ON KODE=1
  int ND=N;
  int NUD=NU;

  double S1=0;
  double S2=0;
  NZ=0;
  int MZ=0;

  int NN=std::min(2,ND); 
  double FNU=NUD;
  double FN=NUD+ND-1;
  double FNN=FN;
  bool goto300 = false;
  if(FN<2.0) goto300 = true;

  // OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
  // FOR THE LAST ORDER, NU+N-1.GE.NULIM
  double ZN=xx/FN;
  double RTZ=std::sqrt( 1.+ZN*ZN );
  double GLN=std::log( (1.+RTZ)/ZN );
  double T=RTZ*( 1.-ETX )+ETX/( ZN+RTZ );
  double CN=-FN*( T-GLN );
  if( CN>ELIM ) return;
  bool goto200 = false;
  if( NUD<NULIM[NN-1] ) goto200 = true;
  if( NN!=1 && !goto200 ) // if nn==1 goto 10 / No calculation, go while loop
  {
    FN=FNU;
    ZN=xx/FN;
    RTZ=std::sqrt( 1.+ZN*ZN );
    GLN=std::log( (1.+RTZ)/ZN );
    T=RTZ*(1.-ETX)+ETX/(ZN+RTZ);
    CN=-FN*(T-GLN);
  }

  while( 5==5 ) // 5 CONTINUE
  {
    // UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
    // FOR THE FIRST ORDER, NU.GE.NULIM

    if( !goto300 )
    {
      double TRX = 0.;
      double TM = 0.;
      // label 10
      bool goto250 = false;
      if( !goto200 )
      {
        if( CN<-ELIM ) // goto 95
        {
          while( 95==95 )
          {
            NUD=NUD+1;
            ND=ND-1;
            if( ND==0 ) return; // goto 96 shortcut
            NN=std::min(2,ND);
            FNU=NUD;
            if( FNN<2. ) continue; // goto 95
            if( NUD<NULIM[NN-1] ) break;// goto 200
          }
        }
        else
        {
          // label 20 NOT USED
          // ASYMPTOTIC EXPANSION FOR ORDERS NU AND NU+1.GE.NULIM
          //                            ASKBES( KODE,FNU,NN,xx,RTZ,CN,Y ); //
          //                            No need to implement
          // goto (96,27), NN
          if( NN==1 )
          {
            goto96( NZ, N, ND, Y );
            yy = Y[0];
            return;
          }
          if( NN==2 )
          {
            TRX=2./xx;
            TM=TRX*(FNU+1.);
            goto250 = true;
          }
        }
      }

      // label 200
      while( 200==200 )
      {
        bool goto251 = false;
        if( !goto250 )
        {
          if( KODE!=2 ) // if KODE==2 goto 201, No calculation below
          {
            if( xx>ELIM ) // label 204
            {
              while( 95==95 )
              {
                NUD=NUD+1;
                ND=ND-1;
                if( ND==0 ) return; // goto 96
                NN=std::min(2,ND);
                FNU=NUD;
                if( FNN<2. ) continue;// goto 95
                if( NUD<NULIM[NN-1] ) break;// goto 200
              }
            }
          }
          // label 201
          S1=BESK01( xx,0,KODE,MZ ); // No MZ USED
          S2=BESK01( xx,1,KODE,MZ ); // No MZ USED
          TRX=2./xx;
          TM=TRX;
          // HERE NU=0 OR 1 IMPLIES N.GE.2. OR NU.GE.2 IMPLIES N.GE.1
          int IN=NUD-1;
          // FORWARD RECUR FROM 0 TO NU TO GET Y(1)
          double S = 0.;
          if( IN>0 )
          {
            // label 203
            for( int i=0; i<IN; i++ )
            {
              S=S2;
              S2=TM*S2+S1;
              S1=S;
              TM=TM+TRX;
            }
          }
          // FORWARD RECUR FROM NU TO NU+1 TO GET Y(2)
          if( IN>=0 )
          {
            // label 211
            if( ND==1 ) goto251 = true;
            if( ND!=1 ) // if ND==1 goto 251
            {
              S=S2;
              S2=TM*S2+S1;
              S1=S;
              TM=TM+TRX;
              Y[0]=S1;
              Y[1]=S2;
            }
          }
          else
          {
            // label 220
            // FORWARD RECUR FROM NU+2 TO NU+N-1
            Y[0]=S1;
            Y[1]=S2;
          }
        }

        if( ND==2 ) 
        {
          goto96( NZ, N, ND, Y );
          yy = Y[0];
          return;
        }
        int NM2=ND-2;
        for( int i=0; i<NM2; i++ )
        {
          Y[i+2]=TM*Y[i+1]+Y[i];
          TM=TM+TRX;
        }
        if( goto251 ) Y[0]=S2;
        goto96( NZ, N, ND, Y );
        yy = Y[0];
        return;

      }

    }

    // FN.LT.2 IMPLIES THE CASES NU=0 , N=1 OR 2
    // NU=1 , N=1
    // label 300
    if( KODE==1 )
    {
      // label 309
      if( xx>ELIM )
      {
        while( 95==95 )
        {
          NUD=NUD+1;
          ND=ND-1;
          if( ND==0 ) 
          {
            goto96( NZ, N, ND, Y );
            yy = Y[0];
            return; // goto 96
          }
          NN=std::min(2,ND);
          FNU=NUD;
          if( FNN<2. ) continue; // goto 95
          if( NUD<NULIM[NN-1] ) break;// goto 200
        }
      }
    }

    // label 310
    int NUM1=NUD-1;
    for( int i=0; i<ND; i++ )
    {
      int J=NUM1+i;
      double ANS=BESK01( xx,J,KODE,MZ ); // No MZ USED
      Y[i]=ANS;
    }
    goto96( NZ, N, ND, Y );
    yy = Y[0];
    return; // goto 96

  }

}

void LocalEquilibriumModel::goto96( int &NZ, int N, int ND, std::vector< double > &Y )
{
  NZ=N-ND;
  if(NZ==0) return;
  if(ND!=0) 
  {
    for( int i=0; i<ND; i++ )
    {
      int J=N-i+1;
      int K=ND-i+1;
      Y[J]=Y[K];
    }
  }
  for( int i=0; i<NZ; i++ ) Y[i]=0;
}


double LocalEquilibriumModel::BESI01( double xx, int NU, int KODE )
{
  double ELIM = 667e0;
  int N1 = 18; int N2 = 19;
  int M1 = 16; int M2 = 17;
  double AI0[18] = { 4.09597733866566e+00, 4.61798815304736e+00,
    1.93239151962331e+00, 5.12072928951984e-01, 1.19244431087252e-01,
    2.04767667749564e-02, 3.29580433784138e-03, 4.17924020345661e-04,
    5.13658909026708e-05, 5.16036527683877e-06, 5.12886756184602e-07,
    4.26546791785515e-08, 3.55823254642485e-09, 2.52433036418873e-10,
    1.81419096716532e-11, 1.12207674472741e-12, 7.08316386761487e-14,
    3.88308114508344e-15 };
  double BI0[18] = { 1.70828618312063e-01,-3.11054153985637e-02,
    4.27652251957381e-03,-6.60652186435765e-04, 1.08431306444222e-04,
    -1.84385649051770e-05, 3.18457963076919e-06,-5.48115448669346e-07,
    9.24224168162523e-08,-1.50625048871349e-08, 2.35026532228336e-09,
    -3.49041942588856e-10, 4.91824029507119e-11,-6.56708591744622e-12,
    8.30939305382694e-13,-9.97106489731714e-14, 1.13614195790750e-14,
    -1.23110947393343e-15 };
  double CI0[19] = { 4.02245205507054e-01, 3.36911647825569e-03,
    6.88975834691682e-05, 2.89137052083476e-06, 2.04891858946905e-07,
    2.26666899049816e-08, 3.39623202570920e-09, 4.94060238821134e-10,
    1.18891471071382e-11,-3.14991652807488e-11,-1.32158118398615e-11,
    -1.79417853277908e-12, 7.18012446526379e-13, 3.85277837074653e-13,
    1.54008635387642e-14,-4.15056915510600e-14,-9.55484713964049e-15,
    3.81168204982708e-15, 1.77256039316563e-15 };
  double AI1[18] = { 1.11729135052299e+00, 8.95198877808018e-01,
    3.36971231909350e-01, 7.24336630339257e-02, 1.53300913475520e-02,
    2.26237207566113e-03, 3.36089476871761e-04, 3.78329005043256e-05,
    4.34367081734334e-06, 3.95545350342659e-07, 3.70599583606154e-08,
    2.83473074207374e-09, 2.24494698516683e-10, 1.48067027975197e-11,
    1.01588206538587e-12, 5.89015292859435e-14, 3.56565058329046e-15,
    1.84456468936977e-16 };
  double BI1[18] = { 1.54236847772347e-01,-2.20987741718608e-02,
    2.19839683870564e-03,-2.07055998880429e-04, 1.22621848833534e-05,
    1.43239687327040e-06,-7.91717404214453e-07, 2.15368545162913e-07,
    -4.70664774819712e-08, 9.04215441510569e-09,-1.57608499172581e-09,
    2.52964136818147e-10,-3.77065425810734e-11, 5.25015808110040e-12,
    -6.85889299579377e-13, 8.43830741608493e-14,-9.80745080362827e-15,
    1.07990041545005e-15 };
  double CI1[19] = { 3.89288117509140e-01,-9.76109749136147e-03,
    -1.10588938762624e-04,-3.88256480887769e-06,-2.51223623787022e-07,
    -2.63146884688956e-08,-3.83538038596299e-09,-5.58974346220626e-10,
    -1.89749581225807e-11, 3.25260358287477e-11, 1.41258074355730e-11,
    2.03562854598022e-12,-7.19855176877244e-13,-4.08355112068096e-13,
    -2.10154176162218e-14, 4.27244003763443e-14, 1.04202766179381e-14,
    -3.81440283617941e-15,-1.88035459733009e-15 };

  double result = -999.;
  if( KODE<1 || KODE>2 ) return result;
  if( NU<0 || NU>1 ) return result;
  double AX = std::abs( xx );
  int IK=NU+1;

  if( IK==1 )
  {
    //     I/SUB(0)/(X) BESSEL FUNCTION
    //IF(AX.GT.4D0) GO TO 110
    if( AX<=4e0 )
    {
      double TT=AX-2e0;
      double T=TT*.5e0;
      //int J=N1; // FORTRAN C++
      int J=N1-1; // For C++
      double F1=AI0[J];
      double F2=0e0;
      for( int i=0; i<M1; i++ )
      {
        J=J-1;
        double TEMP1=F1;
        F1=TT*F1-F2+AI0[J];
        F2=TEMP1;
      }
      double ANS=T*F1-F2+AI0[0];
      //GO TO (106,107), KODE
      if( KODE==2 )
      {
        // label 107 
        result=ANS*std::exp( -AX );
        return result;
      }
      if( KODE==1 )
      {
        // label 106 
        result=ANS;
        return result;
      }
    }
    else
    {
      // labe 110
      // if(AX.GT.8.) GO TO 120
      if( AX<=8. )
      {
        double TT=AX-6e0;
        double T=TT*.5e0;
        //int J=N1; // FORTARN
        int J=N1-1; // For C++
        double F1=BI0[J];
        double F2=0e0;
        for( int i=0; i<M1; i++ )
        {
          J=J-1;
          double TEMP1=F1;
          F1=TT*F1-F2+BI0[J];
          F2=TEMP1;
        }
        double ANS=T*F1-F2+BI0[0];
        //GO TO (116,106),KODE
        if( KODE==1 )
        {
          result=ANS*std::exp( AX );
          return result;
        }
        if( KODE==2 )
        {
          result=ANS;
          return result;
        }
      }
      else
      {
        //labe 120
        double T=16e0/AX-1e0;
        double TT=T+T;
        //int J=N2; // FORTRAN
        int J=N2-1; // For C++
        double F1=CI0[J];
        double F2=0e0;
        for( int i=0; i<M2; i++ )
        {
          J=J-1;
          double TEMP1=F1;
          F1=TT*F1-F2+CI0[J];
          F2=TEMP1;
        }
        double ANS=(T*F1-F2+CI0[1])/std::sqrt( AX );
        //GO TO (126,106),KODE
        if( KODE==1 )
        {
          if( AX>ELIM ) return -999;
          result=ANS*std::exp( AX );
          return result;
        }
        if( KODE==2 )
        {
          result=ANS;
          return result;
        }
      }

    }

  }
  if( IK==2 )
  {
    //     I/SUB(1)/(X) BESSEL FUNCTION
    //IF(AX.GT.4D0) GO TO 210
    if( AX<=4e0 )
    {
      double TT=AX-2e0;
      double T=TT*.5e0;
      //int J=N1; // FORTRAN C++
      int J=N1-1; // For C++
      double F1=AI1[J];
      double F2=0e0;
      for( int i=0; i<M1; i++ )
      {
        J=J-1;
        double TEMP1=F1;
        F1=TT*F1-F2+AI1[J];
        F2=TEMP1;
      }
      double ANS=T*F1-F2+AI1[0];
      //GO TO (106,107), KODE
      if( KODE==2 )
      {
        // label 107 
        result=ANS*std::exp( -AX );
        return result;
      }
      if( KODE==1 )
      {
        // label 106 
        result=ANS;
        return result;
      }
    }
    else
    {
      // labe 210
      // if(AX.GT.8.) GO TO 120
      if( AX<=8. )
      {
        double TT=AX-6e0;
        double T=TT*.5e0;
        //int J=N1; // FORTARN
        int J=N1-1; // For C++
        double F1=BI1[J];
        double F2=0e0;
        for( int i=0; i<M1; i++ )
        {
          J=J-1;
          double TEMP1=F1;
          F1=TT*F1-F2+BI1[J];
          F2=TEMP1;
        }
        double ANS=T*F1-F2+BI1[0];
        if( KODE==1 )
        {
          // label 216
          result=ANS*std::exp( AX );
          return result;
        }
        if( KODE==2 )
        {
          // label 217
          result=ANS;
          return result;
        }
      }
      else
      {
        //labe 220
        double T=16e0/AX-1e0;
        double TT=T+T;
        //int J=N2; // FORTRAN
        int J=N2-1; // For C++
        double F1=CI1[J];
        double F2=0e0;
        for( int i=0; i<M2; i++ )
        {
          J=J-1;
          double TEMP1=F1;
          F1=TT*F1-F2+CI1[J];
          F2=TEMP1;
        }
        double ANS=(T*F1-F2+CI1[0])/std::sqrt( AX );
        if( KODE==1 )
        {
          // label 226
          if( AX>ELIM ) return -999;
          result=ANS*std::exp( AX );
          return result;
        }
        if( KODE==2 )
        {
          // label 217
          result=ANS;
          return result;
        }
      }

    }

  }
  return result;
}

double LocalEquilibriumModel::BESK01( double xx, int NU, int KODE, int &MZ )
{
  double ELIM = 667e0;
  int N1=15; int N2=21; int N3=14; int N4=22;
  int M1=13; int M2=19; int M3=12; int M4=20;

  double EMLT = -1.15931515658412e-1;
  double  AK0[15] = { 4.89766239821850e-01, 6.85071478635282e-01,
    2.20263147944787e-01, 2.94143920648988e-02, 4.78389179259074e-03,
    3.65095673281410e-04, 3.96552105234993e-05, 2.13571244350405e-06,
    1.74123741003729e-07, 7.26649791536202e-09, 4.74231032016171e-10,
    1.61749721604899e-11, 8.80111420783973e-13, 2.53991582114493e-14,
    1.18507018244790e-15 };
  double  BK0[21] = { 1.21101810253138e+00, 1.66889519207676e-02,
    -3.32975981412374e-03, 6.70802900717429e-04,-1.36272865384902e-04,
    2.78860869644180e-05,-5.74300240638451e-06, 1.18942294775776e-06,
    -2.47573075834394e-07, 5.17613906627857e-08,-1.08653503852187e-08,
    2.28899509247299e-09,-4.83794289746110e-10, 1.02556286108303e-10,
    -2.17990101893024e-11, 4.64499519650372e-12,-9.92020238192219e-13,
    2.12307482474184e-13,-4.55250947249326e-14, 9.77948196627543e-15,
    -2.10428853630363e-15 };
  double  CK0[14] = { 1.23878593820170e+00,-1.41756153741679e-02,
    3.37809034043106e-04,-1.39201655756935e-05, 7.92815326397551e-07,
    -5.64584651140153e-08, 4.75381608069133e-09,-4.56892133422380e-10,
    4.89340910679015e-11,-5.74009881252931e-12, 7.27938173682122e-13,
    -9.88017112659449e-14, 1.42381749918418e-14,-2.16443591571975e-15 };
  double  AK1[15] = { 4.82486240842478e-01, 6.05181915156673e-01,
    1.53253828510709e-01, 3.35117014379959e-02, 3.36879962666501e-03,
    4.41042696291495e-04, 2.83189625853098e-05, 2.64548793222978e-06,
    1.25637058646602e-07, 9.12550891350504e-09, 3.44783302431681e-10,
    2.04890221867815e-11, 6.43563932230815e-13, 3.23636474094122e-14,
    8.70458673285249e-16 };
  double  BK1[22] = { 1.38930840300441e+00,-5.69753920597669e-02,
    1.20249909184072e-02,-2.55387533195591e-03, 5.45258387639594e-04,
    -1.16936329705121e-04, 2.51746904352026e-05,-5.43777151950107e-06,
    1.17797109704500e-06,-2.55830051061306e-07, 5.56856483833337e-08,
    -1.21451129231835e-08, 2.65359402057287e-09,-5.80715571337876e-10,
    1.27268495790715e-10,-2.79285412641216e-11, 6.13612949711377e-12,
    -1.34963109251970e-12, 2.97146532041432e-13,-6.54827012807628e-14,
    1.44428209792469e-14,-3.18800891899413e-15 };
  double  CK1[14] = { 1.29837185686353e+00, 4.44491508992785e-02,
    -5.87146375780839e-04, 2.02849692864934e-05,-1.05926607637730e-06,
    7.16046166444819e-08,-5.82233926413546e-09, 5.45714706062071e-10,
    -5.73482672655767e-11, 6.62769161077494e-12,-8.30462152737046e-13,
    1.11605017034540e-13,-1.59497634434667e-14, 2.40744099380385e-15 };

  double result = -999.;
  if( KODE<1 || KODE>2 ) return result;
  if( NU<0 || NU>1 ) return result;
  if( xx<=0. ) return result;

  int IK=NU+1;
  MZ=0;

  // goto(100,200) IK
  if( IK==1 )
  {
    //     K/SUB(0)/(X)  BESSEL FUNCTION
    //if(X.GT.2D0) GO TO 110
    if( xx<=2e0 ) 
    {
      double T=xx-1.;
      double TT=T+T;
      //int J=N1; // FORTRAN
      int J=N1-1; // For C++
      double F1=AK0[J];
      double F2=0.;
      for( int i=0; i<M1; i++ )
      {
        J=J-1;
        double TEMP1=F1;
        F1=TT*F1-F2+AK0[J];
        F2=TEMP1;
      }
      double BNS=BESI01( xx,0,1 );
      double ANS=-(EMLT+std::log( xx ))*BNS+(T*F1-F2+AK0[0]);
      if( KODE==1 )
      {
        // label 106
        result=ANS;
        return result;
      }
      if( KODE==2 )
      {
        // label 107
        result=ANS*std::exp( xx );
        return result;
      }
    }
    else
    {
      // label 110
      //if( X.GT.5e0 ) GO TO 120
      if( xx<=5e0 ) 
      {
        double T=(xx+xx-7e0)/3e0;
        double TT=T+T;
        //int J=N2; // FORTRAN
        int J=N2-1; // For C++
        double F1=BK0[J];
        double F2=0.;
        for( int i=0; i<M2; i++ )
        {
          J=J-1;
          double TEMP1=F1;
          F1=TT*F1-F2+BK0[J];
          F2=TEMP1;
        }
        double ANS=(T*F1-F2+BK0[0])/std::sqrt( xx );
        //GO TO (116,106),KODE
        if( KODE==1 )
        {
          // label 116
          result=ANS*std::exp( -xx );
          return result;
        }
        if( KODE==2 )
        {
          // label 106
          result=ANS;
          return result;
        }
      }
      else
      {
        // label 120
        double T=10e0/xx-1e0;
        double TT=T+T;
        //int J=N3; // FORTRAN
        int J=N3-1; // For C++ 
        double F1=CK0[J];
        double F2=0e0;
        for( int i=0; i<M3; i++ )
        {
          J=J-1;
          double TEMP1=F1;
          F1=TT*F1-F2+CK0[J];
          F2=TEMP1;
        }
        double ANS=(T*F1-F2+CK0[0])/std::sqrt( xx );
        //GO TO (126,106),KODE
        if( KODE==1 )
        {
          // label 126
          if( xx>ELIM ) 
          {
            MZ=0;
            return 0.;
          }
          result=ANS*std::exp( -xx );
          return result;
        }
        if( KODE==2 )
        {
          // label 106
          result=ANS;
          return result;
        }
      }

    }
  }
  if( IK==2 )
  {
    //     K/SUB(1)/(X)  BESSEL FUNCTION

    //if(X.GT.2D0) GO TO 210
    if( xx<=2e0 ) 
    {
      double T=xx-1.;
      double TT=T+T;
      //int J=N1; // FORTRAN
      int J=N1-1; // For C++
      double F1=AK1[J];
      double F2=0.;
      for( int i=0; i<M1; i++ )
      {
        J=J-1;
        double TEMP1=F1;
        F1=TT*F1-F2+AK0[J];
        F2=TEMP1;
      }
      double BNS=BESI01( xx,1,1 );
      double ANS=-(EMLT+std::log( xx ))*BNS+(T*F1-F2+AK1[0]);
      if( KODE==1 )
      {
        // label 106
        result=ANS;
        return result;
      }
      if( KODE==2 )
      {
        // label 107
        result=ANS*std::exp( xx );
        return result;
      }
    }
    else
    {
      // label 210
      //if( X.GT.5e0 ) GO TO 120
      if( xx<=5e0 ) 
      {
        double T=(xx+xx-7e0)/3e0;
        double TT=T+T;
        //int J=N4; // FORTRAN
        int J=N4-1; // For C++
        double F1=BK1[J];
        double F2=0.;
        for( int i=0; i<M4; i++ )
        {
          J=J-1;
          double TEMP1=F1;
          F1=TT*F1-F2+BK1[J];
          F2=TEMP1;
        }
        double ANS=(T*F1-F2+BK1[0])/std::sqrt( xx );
        //GO TO (116,106),KODE
        if( KODE==1 )
        {
          // label 116
          result=ANS*std::exp( -xx );
          return result;
        }
        if( KODE==2 )
        {
          // label 106
          result=ANS;
          return result;
        }
      }
      else
      {
        // label 120
        double T=10e0/xx-1e0;
        double TT=T+T;
        //int J=N3; // FORTRAN
        int J=N3-1; // For C++ 
        double F1=CK1[J];
        double F2=0e0;
        for( int i=0; i<M3; i++ )
        {
          J=J-1;
          double TEMP1=F1;
          F1=TT*F1-F2+CK1[J];
          F2=TEMP1;
        }
        double ANS=(T*F1-F2+CK1[0])/std::sqrt( xx );
        //GO TO (226,106),KODE
        if( KODE==1 )
        {
          // label 226
          if( xx>ELIM ) 
          {
            MZ=0;
            return 0.;
          }
          result=ANS*std::exp( -xx );
          return result;
        }
        if( KODE==2 )
        {
          // label 106
          result=ANS;
          return result;
        }
      }

    }
  }

  return result;
}

void LocalEquilibriumModel::POISSON( double AV, int &M )
{
  // --------------- FORTRAN -----------------
  M = 0;
  if( AV<=0. ) return;

  double ran = double(gRandom->Uniform());
  double COMP=ran*std::exp( AV );
  double SUM=1e0;
  double DSUM=1e0;

  while( 10==10 )
  {
    if( COMP<=SUM ) return;
    M=M+1;
    DSUM=DSUM*AV/M;
    SUM=SUM+DSUM;
  }

  // --------------- FORTRAN -----------------

  // --------------- ROOT -----------------
  // --------------- ROOT -----------------
}


void LocalEquilibriumModel::LOREN( double BGX, double BGY, double BGZ, double G, double PCX, double PCY, double PCZ, double EC, double &PX, double &PY, double &PZ, double &E )
{
  double BGPC=BGX*PCX+BGY*PCY+BGZ*PCZ;
  double EG=EC+BGPC/(G+1.);
  PX=PCX+EG*BGX;
  PY=PCY+EG*BGY;
  PZ=PCZ+EG*BGZ;

  E=G*EC+BGPC;
}

double LocalEquilibriumModel::GetNuclearMass( int Z, int A )
{
  double mass = 0.;
  double mass_p = 938.272;
  double mass_n = 939.565;

  int N = A - Z;

  // Weizsacker-Bethe
  double Av = 16;
  double As = 17;
  double Ac = 0.7;
  double Asym = 23;

  double BE = Av * A
    - As * TMath::Power( double ( A ), 2./3. )
    - Ac * Z*Z/TMath::Power( double ( A ), 1./3. )
    - Asym * ( N - Z )* ( N - Z ) / A;

  mass = Z * mass_p + N * mass_n - BE;

  return mass;
}

void LocalEquilibriumModel::EmitParticlesInCM( std::vector<int> *fZ, std::vector<int> *fA, std::vector<double> *fMass, std::vector<TLorentzVector> *fVector )
{
  fVector->resize( 0 );
  fZ->resize( 0 );
  fA->resize( 0 );
  fMass->resize( 0 );


  // Usual case
  double e_beam = mass_minus_be + ke_beam;
  double p_beam = std::sqrt( (e_beam+mass_minus_be) * (e_beam-mass_minus_be) );
  double y_beam = std::log( (e_beam+p_beam) / mass_minus_be );
  double y_cm   = .5*y_beam; // beam in cm


  //  BEAM-AXIS FLOW
  double yLongCut = 0.81*y_cm; // cutoff longitudinal flow Rectangle Convoluted
  double rapidity_long_width  = .10*y_cm; // gaussian overlay f/longitudinal w/Gaussian
  double del = .75;
  double ys = yLongCut / TMath::ACosH(1.+del);
  double facyr = TMath::SinH( yLongCut/ys );
  //yLongCut = .4*y_cm; // Original Flow
  //rapidity_long_width  = .22*y_cm; // Original Flow

  //  TRANSVERSE FLOW
  double yTransverse = .22*y_cm;  // amplitude
  //double yTransverse = .0*y_cm;  // NO FLOW
  double yTanSaturation = .40*y_cm; // saturation scale Hyperbolic Tangent?

  //  RADIAL VELOCITY
  double rapidity_radial = .60*y_cm;
  double elliptic_asymmetry = 0.16; // Elliptic asymmetry
  //double elliptic_asymmetry = 0.; // NO FLOW asymmetry
  double rapidity_x = (1.-.5*elliptic_asymmetry)*rapidity_radial;
  double rapidity_y = (1.+.5*elliptic_asymmetry)*rapidity_radial; // Transverse Radial Velocity Convolution

  double ke_beam_cm = mass_minus_be*std::pow( std::sinh(y_cm), 2 )/( 1.+std::cosh(y_cm) ); //      Beam Energy in CM
  double temp = (2./3.)*ke_beam_cm*.45 + .007; // estimate Coarse

  double temp_max = ke_beam_cm;
  double temp_min = 0.008;


  int multi_tot = 0;
  int multi_p, multi_d, multi_t, multi_h3, multi_h4;

  POISSON( avg_multi_p, multi_p );
  multi_p = std::min( multi_p, multi_max );      // capacity f/one species

  POISSON( avg_multi_d, multi_d );
  multi_d = std::min( multi_d, multi_max );      // capacity f/one species

  POISSON( avg_multi_t, multi_t );
  multi_t = std::min( multi_t, multi_max );      // capacity f/one species

  POISSON( avg_multi_h3, multi_h3 );
  multi_h3 = std::min( multi_h3, multi_max );      // capacity f/one species

  POISSON( avg_multi_h4, multi_h4 );    // selection of real multiplicities for one event
  multi_h4 = std::min( multi_h4, multi_max );      // capacity f/one species

  if( multi_tot>multi_max2 ) // trimming multiplicity if no space
  {
    multi_p  = multi_p  * multi_max2/multi_tot;
    multi_d  = multi_d  * multi_max2/multi_tot;
    multi_t  = multi_t  * multi_max2/multi_tot;
    multi_h3 = multi_h3 * multi_max2/multi_tot;
    multi_h4 = multi_h4 * multi_max2/multi_tot;
  }
  multi_tot = multi_p + multi_d + multi_t + multi_h3 + multi_h4;

  std::vector< int > particle_multi;
  particle_multi.resize(nParticleType);
  particle_multi[0] = multi_p;
  particle_multi[1] = multi_d;
  particle_multi[2] = multi_t;
  particle_multi[3] = multi_h3;
  particle_multi[4] = multi_h4;

  for( int particle_type=0; particle_type<nParticleType; particle_type++ )
  {
    double mass  = particle_mass[particle_type];
    int multi = particle_multi[particle_type];
    int IAI   = particle_a[particle_type];
    int IZI   = particle_z[particle_type];

    int IP = -1;
    double mass_square = mass*mass;
    for( int I=0; I<multi; I++ ) // line 235
    {
      //  BOOKKEEPING
      IP=IP+1;


      //  radial + elliptic on top of local
      //  RAPIDITY
      double yLongitude = ys * TMath::ASinH( gRandom->Uniform(-1., 1.) * facyr ); // Modified flow
      yLongitude = yLongitude + rapidity_long_width*gRandom->Gaus(0., 1.); // Modified flow


      double y_ratio = std::abs( yLongitude )/y_cm;
      double fag = (1.-1./3.)*y_ratio/(1.-1./3.*y_ratio);
      double FACT = 0;
      if( y_ratio < 1. )
        FACT = .5*( std::cos(fag*TMath::Pi()) + 1. );
      else
        FACT = 0.;

      temp = temp_min + ( temp_max - temp_min )*FACT;
      temp = std::max( temp_min, temp );

      double p_tot  = 0;
      double p_x = 0;
      double p_y = 0;
      double p_z = 0;
      GetMomentum( p_x, p_y, p_z, p_tot, mass, temp );


      double e_cm = std::sqrt( p_tot*p_tot+mass_square ); // local momentum
      //  RADIAL/ELLIPTIC
      /*
      double phir = gRandom->Uniform( 0, TMath::TwoPi() );
      double cos_phi = std::cos( phir );
      double sin_phi = std::sin( phir );

      double yRadial = 1./std::sqrt( (cos_phi/rapidity_x)*(cos_phi/rapidity_x) + (sin_phi/rapidity_y)*(sin_phi/rapidity_y) );
      //double betaRadial = std::tanh( yRadial ); // No use
      double gammaRadial = std::cosh( yRadial );
      double betaGammaRadial = std::sinh( yRadial );
      double betaGammaX = betaGammaRadial*cos_phi;
      double betaGammaY = betaGammaRadial*sin_phi;

      LOREN( betaGammaX, betaGammaY, 0., gammaRadial, p_x, p_y, p_z, e_cm, pr_x, pr_y, pr_z, er );
      */

      double YCOLW_MAX = y_cm*0.33;
      double YCOLW = YCOLW_MAX * FACT;
      double YCOL = YCOLW * gRandom->Gaus( 0., 1. );
      double GCOL = std::cosh(YCOL);
      double BGCOL = std::sinh(YCOL);

      double pr_x, pr_y, pr_z, er;
      //LOREN( betaGammaX, betaGammaY, 0., gammaRadial, p_x, p_y, p_z, e_cm, pr_x, pr_y, pr_z, er );
      LOREN( 0., BGCOL, 0., GCOL, p_x, p_y, p_z, e_cm, pr_x, pr_y, pr_z, er );





      double gammaLong  = std::cosh( yLongitude );
      double betaGammaLong = std::sinh( yLongitude );

      double pl_x, pl_y, pl_z, el;
      LOREN( 0., 0., betaGammaLong, gammaLong, pr_x, pr_y, pr_z, er, pl_x, pl_y, pl_z, el );

      //  longitudinal added
      //  SIDEWARD
      double ySide = yTransverse * std::tanh( yLongitude/yTanSaturation );
      double gammaSide = std::cosh( ySide );
      double betaGammaSide = std::sinh( ySide );

      double final_px = 0;
      double final_py = 0;
      double final_pz = 0;
      double final_e = 0;
      LOREN( betaGammaSide, 0., 0., gammaSide, pl_x, pl_y, pl_z, el, final_px, final_py, final_pz, final_e );

      double ycm = TMath::Log( (final_e+final_pz) / TMath::Sqrt( mass*mass + final_px*final_px + final_py*final_py ) );
      double y0 = ycm/y_cm; // Let it be y0, yCM/yBeam
      double pt = TMath::Sqrt( final_px*final_px + final_py*final_py );
      double mt = TMath::Sqrt( pt*pt + mass*mass );
      final_e  = mt * TMath::CosH( y0 );
      final_pz = mt * TMath::SinH( y0 );

      TLorentzVector thisVector( final_px, final_py, final_pz, final_e );
      fVector->push_back( thisVector );

      fZ->push_back( IZI );
      fA->push_back( IAI );
      fMass->push_back( mass );

    } // end of loop over specific species
  }

}

void LocalEquilibriumModel::CalcAvgMulti()
{
  double ratio_to_breakup   = 1./6.;  // break-up density in units of normal
  double saturation_density = .16;
  double hbc                = .19733; // H -> hbc / h bar c

  double AG = 0, BESK2N = 0, BESK2D = 0, BESK2TE = 0, BESK24E = 0;

  double e_beam = mass_minus_be + ke_beam;
  double p_beam = std::sqrt( (e_beam+mass_minus_be) * (e_beam-mass_minus_be) );
  double y_beam = std::log( (e_beam+p_beam) / mass_minus_be );
  double y_cm   = .5*y_beam; // beam in cm

  //  BEAM-AXIS FLOW
  double yLongCut = .4*y_cm; // cutoff longitudinal flow Rectangle Convoluted
  double rapidity_long_width  = .22*y_cm; // gaussian overlay f/longitudinal w/Gaussian

  //  TRANSVERSE FLOW
  double yTransverse = .2*y_cm;  // amplitude
  double yTanSaturation = .45*y_cm; // saturation scale Hyperbolic Tangent?

  //  RADIAL VELOCITY
  double rapidity_radial = .38*y_cm;
  double elliptic_asymmetry = 0.16; // Elliptic asymmetry
  double rapidity_x = (1.-.5*elliptic_asymmetry)*rapidity_radial;
  double rapidity_y = (1.+.5*elliptic_asymmetry)*rapidity_radial; // Transverse Radial Velocity Convolution

  //  temp ESTIMATE
  double ke_beam_cm = mass_minus_be*std::pow( std::sinh(y_cm), 2 )/( 1.+std::cosh(y_cm) ); //      Beam Energy in CM
  double temp = (2./3.)*ke_beam_cm*.45 + .007; // estimate Coarse

  //  ESTIMATE OF CHEMICAL POT
  double AN = saturation_density*ratio_to_breakup;   // breakup density in 1/fm^3

  int NZ = 0;
  AG = mass_avg/temp;             // going f/relativistic
  BESKN( AG, 2, 2, 1, BESK2N, NZ );
  double EMIT = AN/4.*2. * TMath::Pi()*TMath::Pi() * hbc*hbc*hbc / (mass_avg*mass_avg * temp*BESK2N);    // see notes: exp((mu-m)/T) relat

  AG = mass_a2/temp;             // deuteron now
  BESKN( AG, 2, 2, 1, BESK2D, NZ );
  double ratio_d_p = 3./2.*EMIT*std::exp(be/temp) * (mass_a2/mass_avg)*(mass_a2/mass_avg) * BESK2D/BESK2N;

  //  DEUTERON/PROTON RATIO
  AG = mass_a3/temp;             // triton/helion
  BESKN( AG, 2, 2, 1, BESK2TE, NZ );
  double ratio_t_p = 2./2. * EMIT*EMIT * std::exp(be_3he/temp) * (mass_a3/mass_avg)*(mass_a3/mass_avg) * BESK2TE/BESK2N;

  //  TRITON/PROTON AND HELION/PROTON RATIO
  AG = mass_a4/temp;             // alpha
  BESKN( AG, 2, 2, 1, BESK24E, NZ );
  double ratio_a_p = 1./2. * EMIT*EMIT*EMIT * std::exp(be_4he/temp) * (mass_a4/mass_avg)*(mass_a4/mass_avg) * BESK24E/BESK2N;

  double charges = 1. + ratio_d_p + 2.*ratio_t_p + ratio_t_p + 2.*ratio_a_p; //net average charge assuming p is 1
  double charge_density_p  = 1./charges;
  double charge_density_d  = ratio_d_p  * charge_density_p;
  double charge_density_t  = ratio_t_p  * charge_density_p;
  double charge_density_h = 2.*charge_density_t;
  double charge_density_a = 2.*ratio_a_p * charge_density_p;

  avg_multi_p  = avg_multi * charge_density_p;
  avg_multi_d  = avg_multi * charge_density_d;
  avg_multi_t  = avg_multi * charge_density_t;
  avg_multi_h3 = avg_multi * charge_density_h/2.;
  avg_multi_h4 = avg_multi * charge_density_a/2.; //individual multiplicities


  // Rough multiplicity correction -- J Park
  avg_multi_p  = avg_multi * charge_density_p * 1.2;
  avg_multi_d  = avg_multi * charge_density_d;
  avg_multi_t  = avg_multi * charge_density_t * 3.11;
  avg_multi_h3 = avg_multi * charge_density_h/2. * 0.79;
  avg_multi_h4 = avg_multi * charge_density_a/2. * 6.8; //individual multiplicities
}

void LocalEquilibriumModel::GetMomentum( double &PX, double &PY, double &PZ, double &P, double AM, double T )
{
  double PI=4.*std::atan(1.);

  double TM=T/AM;
  double TMK=TM*TM;

  double PMXK=2.*TM*(TM+std::sqrt( TMK+1. ));
  double PMX=std::sqrt( PMXK );

  double PMK=PMX*PMX;
  double FX=PMK*std::exp( -PMK/(TM*(1.+std::sqrt( PMK+1. ))) );
  double PML=std::sqrt( FX );
  double PMH=PMX*1.8;

  double AIL=FX*PML/3.;
  double AIX=FX*(PMH-PML);

  double FGH=0;
  double AFGPH=0;
  //CALL LFGENP(PMH,FGH,AFGPH);
  {
    double SPMK=std::sqrt( PMH*PMH+1. );
    AFGPH=2./PMH+PMH/(TM*(SPMK+1.))*(-2.+(PMH*PMH)/(SPMK*(SPMK+1.)));
    FGH= PMH*PMH*std::exp( -PMH*PMH/(TM*(1.+std::sqrt( PMH*PMH+1. ))) );
  }
  double AIH=-FGH/AFGPH;

  double AI=AIL+AIX+AIH;
  double RL=AIL/AI;
  double RX=AIX/AI;
  double RH=AIH/AI;

  double FM=0;
  double PM=0;
  while( 100==100 )
  {
    double RINS=double(gRandom->Uniform());
    if( RINS<=RL )
    {
      if(RINS==0.) PM=0.;
      else PM=PML* std::pow( (RINS/RL), (1./3.) );
      FM=PM*PM;
    }
    else if( RINS<=RL+RX )
    {
      PM=PML+(PMH-PML)*(RINS-RL)/RX;
      FM=FX;
    }
    else
    {
      if(RINS>=1.) continue;
      PM=PMH+std::log( (1.-RINS)/RH )/AFGPH;
      FM=FGH*std::exp( AFGPH*(PM-PMH) );
    }
    //if( FM*double(gRandom->Uniform())>FGEN(PM) ) continue; // goto 100
    double FGEN = PM*PM*std::exp( -PM*PM/(TM*(1.+std::sqrt( PM*PM+1. ))) );
    if( !(FM*double(gRandom->Uniform())>FGEN) ) break;
  }

  P=PM*AM;

  double ST=0;
  double CT=2.*(double(gRandom->Uniform())-.5);
  if( std::abs(CT)<1. ) ST=std::sqrt(1.-CT*CT);
  else ST=0.;
  double PHI=gRandom->Uniform( 0, TMath::TwoPi() );

  PZ=P*CT;
  double PT=P*ST;
  PX=PT*std::cos(PHI);
  PY=PT*std::sin(PHI);

}
