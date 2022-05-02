#include "VertexFitterSimple.h"
#include <iostream>

#include "TFile.h"
#include "TString.h"

#include  "MCParticle.h"

using namespace VertexFitterSimple;

bool debug = false;


TVector3 VertexFitterSimple::ParToP(TVectorD Par){
  double fB = 2;  // 2 Tesla
  
  Double_t C    = Par(2);
  Double_t phi0 = Par(1);
  Double_t ct   = Par(4);
  //
  TVector3 Pval;
  Double_t pt = fB*0.2998 / TMath::Abs(2 * C);
  Pval(0) = pt*TMath::Cos(phi0);
  Pval(1) = pt*TMath::Sin(phi0);
  Pval(2) = pt*ct;
  //
  return Pval;
}


TVectorD VertexFitterSimple::XPtoPar(TVector3 x, TVector3 p, Double_t Q){

  double fB = 2;  // 2 Tesla
  
  //
  TVectorD Par(5);
  // Transverse parameters
  Double_t a = -Q*fB*0.2998;                      // Units are Tesla, GeV and meters
  Double_t pt = p.Pt();
  Double_t C = a / (2 * pt);                      // Half curvature
  //cout << "ObsTrk::XPtoPar: fB = " << fB << ", a = " << a << ", pt = " << pt << ", C = " << C << endl;
  Double_t r2 = x.Perp2();
  Double_t cross = x(0)*p(1) - x(1)*p(0);
  Double_t T = TMath::Sqrt(pt*pt - 2 * a*cross + a*a*r2);
  Double_t phi0 = TMath::ATan2((p(1) - a*x(0)) / T, (p(0) + a*x(1)) / T); // Phi0
  Double_t D;                                                     // Impact parameter D
  if (pt < 10.0) D = (T - pt) / a;
  else D = (-2 * cross + a*r2) / (T + pt);
  //
  Par(0) = D;             // Store D
  Par(1) = phi0;  // Store phi0
  Par(2) = C;             // Store C
  //Longitudinal parameters
  Double_t B = C*TMath::Sqrt(TMath::Max(r2 - D*D,0.0) / (1 + 2 * C*D));
  Double_t st = TMath::ASin(B) / C;
  Double_t ct = p(2) / pt;
  Double_t z0 = x(2) - ct*st;
  //
  Par(3) = z0;            // Store z0
  Par(4) = ct;            // Store cot(theta)
  //
  return Par;
}


//
//TH1F* hTry;
//
Double_t VertexFitterSimple::FastRv(TVectorD p1, TVectorD p2){
  //
  // Find radius of intersection between two tracks in the transverse plane
  //
  // p = (D,phi, C)
  //
  // Solving matrix
  TMatrixDSym H(2);
  H(0, 0) = -TMath::Cos(p2(1));
  H(0, 1) =  TMath::Cos(p1(1));
  H(1, 0) = -TMath::Sin(p2(1));
  H(1, 1) =  TMath::Sin(p1(1));
  Double_t Det = TMath::Sin(p2(1) - p1(1));
  H *= 1.0 / Det;
  //
  // Convergence parameters
  Int_t Ntry = 0;
  Int_t NtryMax = 100;
  Double_t eps = 1000.;
  Double_t epsMin = 1.0e-6;
  //
  // Vertex finding loop
  //
  TVectorD cterm(2);
  cterm(0) = p1(0);
  cterm(1) = p2(0);
  TVectorD xv(2);
  Double_t R = 1000.;
  while (eps > epsMin)
    {
      xv = H * cterm;
      Ntry++;
      if (Ntry > NtryMax)
	{
	  std::cout << "FastRv: maximum number of iteration reached" << std::endl;
	  break;
	}
      Double_t Rnew = TMath::Sqrt(xv(0) * xv(0) + xv(1) * xv(1));
      eps = Rnew - R;
      R = Rnew;
      cterm(0) = p1(2) * R * R;
      cterm(1) = p2(2) * R * R;
    }
  //
  return R;
}
TMatrixDSym VertexFitterSimple::RegInv3(TMatrixDSym &Smat0){
  //
  // Regularized inversion of symmetric 3x3 matrix with positive diagonal elements
  //
  TMatrixDSym Smat = Smat0;
  Int_t N = Smat.GetNrows();
  if (N != 3){
    std::cout << "RegInv3 called with  matrix size != 3. Abort & return standard inversion." << std::endl;
    return Smat.Invert();
  }
  TMatrixDSym D(N); D.Zero();
  Bool_t dZero = kTRUE;	// No elements less or equal 0 on the diagonal
  for (Int_t i = 0; i < N; i++) if (Smat(i, i) <= 0.0)dZero = kFALSE;
  if (dZero){
    for (Int_t i = 0; i < N; i++) D(i, i) = 1.0 / TMath::Sqrt(Smat(i, i));
    TMatrixDSym RegMat = Smat.Similarity(D);
    TMatrixDSym Q(2);
    for (Int_t i = 0; i < 2; i++){
      for (Int_t j = 0; j < 2; j++)Q(i, j) = RegMat(i, j);
    }
    Double_t Det = 1 - Q(0, 1)*Q(1, 0);
    TMatrixDSym H(2);
    H = Q;
    H(0, 1) = -Q(0, 1);
    H(1, 0) = -Q(1, 0);
    TVectorD p(2);
    p(0) = RegMat(0, 2);
    p(1) = RegMat(1, 2);
    Double_t pHp = H.Similarity(p);
    Double_t h = pHp-Det;
    //
    TMatrixDSym pp(2); pp.Rank1Update(p);
    TMatrixDSym F = (h*H) - pp.Similarity(H);
    F *= 1.0 / Det;
    TVectorD b = H*p;
    TMatrixDSym InvReg(3);
    for (Int_t i = 0; i < 2; i++)
      {
	InvReg(i, 2) = b(i);
	InvReg(2, i) = b(i);
	for (Int_t j = 0; j < 2; j++) InvReg(i, j) = F(i, j);
      }
    InvReg(2, 2) = -Det;
    //
    InvReg *= 1.0 / h;
    //
    //
    return InvReg.Similarity(D);
  }
  else
    {
      //std::cout << "RegInv3: found negative elements in diagonal. Return standard inversion." << std::endl;
      return Smat.Invert();
    }
}
//
//
//
TMatrixD VertexFitterSimple::Fill_A(TVectorD par, Double_t phi){
  //
  // Derivative of track 3D position vector with respect to track parameters at constant phase 
  //
  // par = vector of track parameters
  // phi = phase
  //
  TMatrixD A(3, 5);
  //
  // Decode input arrays
  //
  Double_t D = par(0);
  Double_t p0 = par(1);
  Double_t C = par(2);
  Double_t z0 = par(3);
  Double_t ct = par(4);
  //
  // Fill derivative matrix dx/d alpha
  // D
  A(0, 0) = -TMath::Sin(p0);
  A(1, 0) = TMath::Cos(p0);
  A(2, 0) = 0.0;
  // phi0
  A(0, 1) = -D*TMath::Cos(p0) + (TMath::Cos(phi + p0) - TMath::Cos(p0)) / (2 * C);
  A(1, 1) = -D*TMath::Sin(p0) + (TMath::Sin(phi + p0) - TMath::Sin(p0)) / (2 * C);
  A(2, 1) = 0.0;
  // C
  A(0, 2) = -(TMath::Sin(phi + p0) - TMath::Sin(p0)) / (2 * C*C);
  A(1, 2) = (TMath::Cos(phi + p0) - TMath::Cos(p0)) / (2 * C*C);
  A(2, 2) = -ct*phi / (2 * C*C);
  // z0
  A(0, 3) = 0.0;
  A(1, 3) = 0.0;
  A(2, 3) = 1.0;
  // ct = lambda
  A(0, 4) = 0.0;
  A(1, 4) = 0.0;
  A(2, 4) = phi / (2 * C);
  //
  return A;
}

//
TVectorD VertexFitterSimple::Fill_a(TVectorD par, Double_t phi){
  //
  // Derivative of track 3D position vector with respect to phase at constant track parameters
  //
  // par = vector of track parameters
  // phi = phase
  //
  TVectorD a(3);
  //
  // Decode input arrays
  //
  Double_t D = par(0);
  Double_t p0 = par(1);
  Double_t C = par(2);
  Double_t z0 = par(3);
  Double_t ct = par(4);
  //
  a(0) = TMath::Cos(phi + p0) / (2 * C);
  a(1) = TMath::Sin(phi + p0) / (2 * C);
  a(2) = ct / (2 * C);
  //
  return a;
}
//

TVectorD VertexFitterSimple::Fill_x0(TVectorD par){
  //
  // Calculate track 3D position at R = |D| (minimum approach to z-axis)
  //
  TVectorD x0(3);
  //
  // Decode input arrays
  //
  Double_t D = par(0);
  Double_t p0 = par(1);
  Double_t C = par(2);
  Double_t z0 = par(3);
  Double_t ct = par(4);
  //
  x0(0) = -D *TMath::Sin(p0);
  x0(1) = D*TMath::Cos(p0);
  x0(2) = z0;
  //
  return x0;
}

//
TVectorD VertexFitterSimple::Fill_x(TVectorD par, Double_t phi){
  //
  // Calculate track 3D position for a given phase, phi
  //
  TVectorD x(3);
  //
  // Decode input arrays
  //
  Double_t D = par(0);
  Double_t p0 = par(1);
  Double_t C = par(2);
  Double_t z0 = par(3);
  Double_t ct = par(4);
  //
  TVectorD x0 = Fill_x0(par);
  x(0) = x0(0) + (TMath::Sin(phi + p0) - TMath::Sin(p0)) / (2 * C);
  x(1) = x0(1) - (TMath::Cos(phi + p0) - TMath::Cos(p0)) / (2 * C);
  x(2) = x0(2) + ct*phi / (2 * C);
  //
  return x;
}



VertexingUtils::FCCAnalysesVertex  VertexFitterSimple::VertexFitter( int Primary, 
								     ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
								     ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
								     bool BeamSpotConstraint,
								     double bsc_sigmax, double bsc_sigmay, double bsc_sigmaz, 
                                                                     double bsc_x, double bsc_y, double bsc_z )  {



  // input = a collection of recoparticles (in case one want to make associations to RecoParticles ?)
  // and thetracks = the collection of all TrackState in the event
  
  VertexingUtils::FCCAnalysesVertex thevertex;
  
  // retrieve the tracks associated to the recoparticles
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks = ReconstructedParticle2Track::getRP2TRK( recoparticles, thetracks );
  
  // and run the vertex fitter
  
  //FCCAnalysesVertex thevertex = VertexFitter_Tk( Primary, tracks, thetracks) ;
  thevertex = VertexFitter_Tk( Primary, tracks,
			       BeamSpotConstraint, bsc_sigmax, bsc_sigmay, bsc_sigmaz, bsc_x, bsc_y, bsc_z );

  //fill the indices of the tracks
  thevertex.reco_ind = get_reco_ind(recoparticles,thetracks);
  
  return thevertex;
}



VertexingUtils::FCCAnalysesVertex  VertexFitterSimple::VertexFitter_Tk( int Primary, 
									ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
                                                                        bool BeamSpotConstraint, 
                                                                        double bsc_sigmax, double bsc_sigmay, double bsc_sigmaz,
									double bsc_x, double bsc_y, double bsc_z )  {
  
  // Units for the beam-spot : mum
  // See https://github.com/HEP-FCC/FCCeePhysicsPerformance/tree/master/General#generating-events-under-realistic-fcc-ee-environment-conditions

  if(debug) std::cout << "Starting VertexFitter_Tk!" << std::endl;


  // final results :
  VertexingUtils::FCCAnalysesVertex TheVertex;
  
  edm4hep::VertexData result;
  ROOT::VecOps::RVec<float> reco_chi2;
  ROOT::VecOps::RVec< TVectorD >  updated_track_parameters;
  ROOT::VecOps::RVec<int> reco_ind;
  ROOT::VecOps::RVec<float> final_track_phases;
  ROOT::VecOps::RVec< TVector3 >  updated_track_momentum_at_vertex;
  
  TheVertex.vertex = result;
  TheVertex.reco_chi2 = reco_chi2;
  TheVertex.updated_track_parameters = updated_track_parameters;
  TheVertex.reco_ind = reco_ind;
  TheVertex.final_track_phases = final_track_phases;
  TheVertex.updated_track_momentum_at_vertex = updated_track_momentum_at_vertex;
  
  
  int Ntr = tracks.size();
  TheVertex.ntracks = Ntr; 
  if ( Ntr <= 1) return TheVertex;   // can not reconstruct a vertex with only one track...
  

  // if a beam-spot constraint is required :
  TMatrixDSym BeamSpotCovI(3);
  TVectorD BeamSpotPos(3);
  if (BeamSpotConstraint) {   // fill in the inverse of the covariance matrix. Convert the units into meters
     BeamSpotCovI(0,0) = 1./pow( bsc_sigmax * 1e-6, 2) ;   // mum to m
     BeamSpotCovI(1,1) = 1./pow( bsc_sigmay * 1e-6, 2) ;   
     BeamSpotCovI(2,2) = 1./pow( bsc_sigmaz * 1e-6, 2) ; 
     BeamSpotPos(0) = bsc_x * 1e-6;
     BeamSpotPos(1) = bsc_y * 1e-6 ;
     BeamSpotPos(2) = bsc_z * 1e-6 ;
  }

  Double_t *final_chi2 = new Double_t[Ntr];
  Double_t *final_phases = new Double_t[Ntr];
  std::vector< TVectorD > final_delta_alpha ;
  TVectorD dummy(5);
  for (int i=0; i < Ntr; i++) {
    final_delta_alpha.push_back( dummy );
  }
  
  
  //
  // Vertex fit (units are meters)
  //
  // Initial variable definitions
  TVectorD x0(3); for (Int_t v = 0; v < 3; v++)x0(v) = 100.; // set to large value
  Double_t Chi2 = 0;
  //
  
  TVectorD x(3);
  TMatrixDSym covX(3);
  
  
  // Stored quantities
  Double_t *fi = new Double_t[Ntr];				// Phases 
  TVectorD **x0i = new TVectorD*[Ntr];			// Track expansion point
  TVectorD **ai = new TVectorD*[Ntr];				// dx/dphi
  Double_t *a2i = new Double_t[Ntr];				// a'Wa
  TMatrixDSym **Di = new TMatrixDSym*[Ntr];		// W-WBW
  TMatrixDSym **Wi = new TMatrixDSym*[Ntr];		// (ACA')^-1
  TMatrixDSym **Winvi = new TMatrixDSym*[Ntr];	// ACA'
  TMatrixD  **Ai = new TMatrixD*[Ntr];            // A
  TMatrixDSym **Covi = new TMatrixDSym*[Ntr];     // Cov matrix of the track parameters
  
  //
  // vertex radius approximation
  // Maximum impact parameter
  Double_t Rd = 0;
  for (Int_t i = 0; i < Ntr; i++)
    {
      //ObsTrk* t = tracks[i];
      //TVectorD par = t->GetObsPar();
      edm4hep::TrackState t = tracks[i] ;
      TVectorD par = VertexingUtils::get_trackParam( t ) ;
      Double_t Dabs = TMath::Abs(par(0));
      if (Dabs > Rd)Rd = Dabs;
    }
  //
  // Find track pair with largest phi difference
  Int_t isel; Int_t jsel; // selected track indices
  Double_t dphiMax = -9999.;	// Max phi difference 
  for (Int_t i = 0; i < Ntr-1; i++)
    {
      //ObsTrk* ti = tracks[i];
      //TVectorD pari = ti->GetObsPar();
      edm4hep::TrackState ti = tracks[i] ;
      TVectorD pari = VertexingUtils::get_trackParam( ti );
      Double_t phi1 = pari(1);
      
      for (Int_t j = i+1; j < Ntr; j++)
	{
	  //ObsTrk* tj = tracks[j];
	  //TVectorD parj = tj->GetObsPar();
	  edm4hep::TrackState tj = tracks[j];
	  TVectorD parj = VertexingUtils::get_trackParam( tj );
	  Double_t phi2 = parj(1);
	  Double_t dphi = TMath::Abs(phi2 - phi1);
	  if (dphi > TMath::Pi())dphi = TMath::TwoPi() - dphi;
	  if (dphi > dphiMax)
	    {
	      isel = i; jsel = j;
	      dphiMax = dphi;
	    }
	}
    }
  //
  // 
  //ObsTrk* t1 = tracks[isel];
  //TVectorD p1 = t1->GetObsPar();
  edm4hep::TrackState t1 = tracks[isel];
  TVectorD p1 = VertexingUtils::get_trackParam( t1 );
  //ObsTrk* t2 = tracks[jsel];
  //TVectorD p2 = t2->GetObsPar();
  edm4hep::TrackState t2 = tracks[jsel];
  TVectorD p2 = VertexingUtils::get_trackParam( t2 );
  Double_t R = FastRv(p1, p2);
  if (R > 1.0) R = Rd;
  R = 0.5 * (R + Rd);
  //
  // Iteration properties
  //
  Int_t Ntry = 0;
  Int_t TryMax = 100;
  if (BeamSpotConstraint) TryMax = TryMax * 5;
  Double_t eps = 1.0e-9; // vertex stability
  Double_t epsi = 1000.;
  //
  while (epsi > eps && Ntry < TryMax)		// Iterate until found vertex is stable
    {
      x.Zero();
      TVectorD cterm(3); TMatrixDSym H(3); TMatrixDSym DW1D(3);
      covX.Zero();	// Reset vertex covariance
      cterm.Zero();	// Reset constant term
      H.Zero();		// Reset H matrix
      DW1D.Zero();
      // 
      for (Int_t i = 0; i < Ntr; i++)
	{
	  // Get track helix parameters and their covariance matrix 
	  //ObsTrk *t = tracks[i];
	  //TVectorD par = t->GetObsPar();
	  //TMatrixDSym Cov = t->GetCov(); 
	  edm4hep::TrackState t = tracks[i] ;
	  TVectorD par = VertexingUtils::get_trackParam( t ) ;
	  TMatrixDSym Cov = VertexingUtils::get_trackCov( t );
	  Covi[i] = new TMatrixDSym(Cov);         // Store matrix
	  Double_t fs;
	  if (Ntry <= 0)	// Initialize all phases on first pass
	    {
	      Double_t D = par(0);
	      Double_t C = par(2);
	      Double_t arg = TMath::Max(1.0e-6, (R*R - D*D) / (1 + 2 * C*D));
	      fs = 2 * TMath::ASin(C*TMath::Sqrt(arg));
	      fi[i] = fs;
	    }
	  //
	  // Starting values
	  //
	  fs = fi[i];								// Get phase
	  TVectorD xs = Fill_x(par, fs);
	  x0i[i] = new TVectorD(xs);				// Start helix position
	  // W matrix = (A*C*A')^-1; W^-1 = A*C*A'
	  TMatrixD A = Fill_A(par, fs);			// A = dx/da = derivatives wrt track parameters
	  Ai[i]  = new TMatrixD(A);       // Store matrix
	  TMatrixDSym Winv = Cov.Similarity(A);	// W^-1 = A*C*A'
	  Winvi[i] = new TMatrixDSym(Winv);		// Store W^-1 matrix
	  TMatrixDSym W = RegInv3(Winv);			// W = (A*C*A')^-1
	  Wi[i] = new TMatrixDSym(W);				// Store W matrix
	  TVectorD a = Fill_a(par, fs);			// a = dx/ds = derivatives wrt phase
	  ai[i] = new TVectorD(a);				// Store a
	  Double_t a2 = W.Similarity(a);
	  a2i[i] = a2;							// Store a2
	  // Build D matrix
	  TMatrixDSym B(3); 
	  B.Rank1Update(a, 1.0);
	  B *= -1. / a2;
	  B.Similarity(W);
	  TMatrixDSym Ds = W+B;					// D matrix
	  Di[i] = new TMatrixDSym(Ds);			// Store D matrix
	  TMatrixDSym DsW1Ds = Winv.Similarity(Ds);	// Service matrix to calculate covX
	  DW1D += DsW1Ds;								
	  // Update hessian
	  H += Ds;
	  // update constant term
	  cterm += Ds * xs;
	}				// End loop on tracks
      //

      TMatrixDSym H0 = H;

      if (BeamSpotConstraint) {
	H += BeamSpotCovI ;
        cterm += BeamSpotCovI * BeamSpotPos ;
        DW1D  += BeamSpotCovI ;
      }

      // update vertex position
      TMatrixDSym H1 = RegInv3(H);
      x = H1*cterm;
      
      // Update vertex covariance
      covX = DW1D.Similarity(H1);

      // Update phases and chi^2
      Chi2 = 0.0;
      for (Int_t i = 0; i < Ntr; i++)
	{
	  TVectorD lambda = (*Di[i])*(*x0i[i] - x);
	  TMatrixDSym Wm1 = *Winvi[i];
	  Double_t addChi2 = Wm1.Similarity(lambda);;
	  //Chi2 += Wm1.Similarity(lambda);
	  Chi2 += addChi2;
	  final_chi2[i] = addChi2;
	  TVectorD a = *ai[i];
	  TVectorD b = (*Wi[i])*(x - *x0i[i]);
	  for (Int_t j = 0; j < 3; j++)fi[i] += a(j)*b(j) / a2i[i];
	  final_phases[i] = fi[i];
	  
	  TMatrixD ta(TMatrixD::kTransposed, *Ai[i]);
	  TMatrixDSym kk(5);
	  kk = *Covi[i];
	  final_delta_alpha[i] =  kk * ta * lambda;  // that's minus delta_alpha
	}
      //

      TVectorD dx = x - x0;
      x0 = x;
      // update vertex stability
      TMatrixDSym Hess = RegInv3(covX);

      epsi = Hess.Similarity(dx);
      Ntry++;
      //if ( Ntry >= TryMax) std::cout << " ... in VertexFitterSimple, Ntry >= TryMax " << std::endl;

      if (BeamSpotConstraint) {
        
        // add the following term to the chi2 :
        TVectorD dx_beamspot = x - BeamSpotPos ;
        Double_t chi2_bsc = BeamSpotCovI.Similarity( dx_beamspot );
        //Chi2 += chi2_bsc -3;
        Chi2 += chi2_bsc ;

      }



      //
      // Cleanup
      //
      for (Int_t i = 0; i < Ntr; i++)
	{
	  x0i[i]->Clear();
	  Winvi[i]->Clear();
	  Wi[i]->Clear();
	  ai[i]->Clear();
	  Di[i]->Clear();
	  Ai[i]->Clear();
	  Covi[i]->Clear();
	  
	  delete x0i[i];
	  delete Winvi[i];
	  delete Wi[i];
	  delete ai[i];
	  delete Di[i];
	  delete Ai[i];
	  delete Covi[i];
	}
    }
  //
  delete[] fi;		// Phases 
  delete[] x0i;		// Track expansion point
  delete[] ai;		// dx/dphi
  delete[] a2i;		// a'Wa
  delete[] Di;		// W-WBW
  delete[] Wi;		// (ACA')^-1
  delete[] Winvi;	// ACA'
  delete[] Ai ;           // A
  delete[] Covi;          // Cov
  
  //
  //return Chi2;
  
  // store the results in an edm4hep::VertexData object
  // go back from meters to millimeters for the units 
  float conv = 1e3;
  std::array<float,6> covMatrix;	// covMat in edm4hep is a LOWER-triangle matrix.
  covMatrix[0] = covX(0,0) * pow(conv,2);
  covMatrix[1] = covX(1,0) * pow(conv,2);
  covMatrix[2] = covX(1,1) * pow(conv,2);
  covMatrix[3] = covX(2,0) * pow(conv,2);
  covMatrix[4] = covX(2,1) * pow(conv,2);
  covMatrix[5] = covX(2,2) * pow(conv,2);
  
  float Ndof = 2.0 * Ntr - 3.0; ;
  
  result.primary = Primary;
  result.chi2 = Chi2 /Ndof ;      // I store the normalised chi2 here
  result.position = edm4hep::Vector3f( x(0)*conv, x(1)*conv, x(2)*conv ) ;  // store the  vertex in mm
  result.covMatrix = covMatrix;
  result.algorithmType = 1;
  
  // Need to fill the associations ...
  
  double scale0 = 1e-3;   //convert mm to m
  double scale1 = 1;
  double scale2 = 0.5*1e3;  // C = rho/2, convert from mm-1 to m-1
  double scale3 = 1e-3 ;  //convert mm to m
  double scale4 = 1.;
  
  scale2 = -scale2 ;   // sign of omega (sign convention)
  
  for (Int_t i = 0; i < Ntr; i++) {
    
    edm4hep::TrackState t = tracks[i] ;
    TVectorD par = VertexingUtils::get_trackParam( t ) ;
    
    // initial momentum :
    //TVector3 ptrack_ini = ParToP( par );
    //std::cout << "----- Track # " << i << " initial track momentum : " << std::endl;
    //ptrack_ini.Print();
    
    // uncomment below to get the post-fit track parameters :
    par -= final_delta_alpha[i] ;
    
    //std::cout << " Track i = " << i << " --- delta_alpha : " << std::endl;
    //final_delta_alpha[i].Print();
    
    // ( px, py, pz) of the track
    TVector3 ptrack = ParToP( par );
    //std::cout << "         updates track param :" << std::endl;
    //ptrack.Print();
    
    // and (px, py) at the vertex instead of the dca :
    double phi0 = par(1);
    double phi = final_phases[i]  ;
    double px_at_vertex = ptrack.Pt() * TMath::Cos( phi0 + phi );
    double py_at_vertex = ptrack.Pt() * TMath::Sin( phi0 + phi );
    TVector3 ptrack_at_vertex( px_at_vertex, py_at_vertex, ptrack.Pz() );
    //std::cout << "         momentum at the vertex : " << std::endl;
    //std::cout << " phi0 at dca = " << phi0 << " phi at vertex = " << phi0+phi << " C = " << par(2) << " phase " << phi << std::endl;
    //ptrack_at_vertex.Print();
    
    updated_track_momentum_at_vertex.push_back( ptrack_at_vertex );
    
    // back to EDM4HEP units...
    par[0] = par[0] / scale0 ;
    par[1] = par[1] / scale1 ;
    par[2] = par[2] / scale2 ;
    par[3] = par[3] / scale3 ;
    par[4] = par[4] / scale4 ;
    updated_track_parameters.push_back( par );
    
    reco_chi2.push_back( final_chi2[i] );
    final_track_phases.push_back( final_phases[i] );
    
  }
  
  TheVertex.vertex = result;
  TheVertex.reco_chi2 = reco_chi2;
  TheVertex.reco_ind = reco_ind;
  TheVertex.updated_track_parameters = updated_track_parameters ;
  TheVertex.updated_track_momentum_at_vertex = updated_track_momentum_at_vertex;
  TheVertex.final_track_phases = final_track_phases;
  
  //std::cout << " end of VertexFitter " << std::endl;
  /*
    for ( Int_t i = 0; i < Ntr; i++) {
    std::cout << " Track #" << i << " chi2 = " << reco_chi2[i] << std::endl;
    std::cout << "        Initial parameters: " << std::endl;
    VertexingUtils::get_trackParam( tracks[i] ).Print();
    std::cout << "        Updated parameters : " << std::endl;
    updated_track_parameters[i].Print();
    }
  */
  
  delete final_chi2;
  delete final_phases;
 
  if(debug) std::cout << "Finished VertexFitter_Tk!" << std::endl;
  return TheVertex;
}


////////////////////////////////////////////////////



ROOT::VecOps::RVec<edm4hep::TrackState>   VertexFitterSimple::get_PrimaryTracks( VertexingUtils::FCCAnalysesVertex  initialVertex,
                                                                        ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
                                                                        bool BeamSpotConstraint,
                                                                        double bsc_sigmax, double bsc_sigmay, double bsc_sigmaz,
                                                                        double bsc_x, double bsc_y, double bsc_z,
                                                                        int ipass  )  {


// iterative procedure to determine the primary vertex - and the primary tracks
// Start from a vertex reconstructed from all tracks, remove the one with the highest chi2, fit again etc

// tracks = the collection of tracks that was used in the first step


if(debug) std::cout << "Starting get_PrimaryTracks!" << std::endl;
float CHI2MAX = 25  ;

if (debug) {
        if (ipass == 0) std::cout << " \n --------------------------------------------------------\n" << std::endl;
        std::cout << " ... enter in VertexFitterSimple::get_PrimaryTracks   ipass = " << ipass <<  std::endl;
        if (ipass == 0) std::cout  << "    initial number of tracks =  " << tracks.size() <<  std::endl;
}

ROOT::VecOps::RVec<edm4hep::TrackState> seltracks = tracks;
ROOT::VecOps::RVec<float> reco_chi2 = initialVertex.reco_chi2;

if ( seltracks.size() <= 1 ) return seltracks;

int isPrimaryVertex = initialVertex.vertex.primary  ;

int maxElementIndex = std::max_element(reco_chi2.begin(),reco_chi2.end()) - reco_chi2.begin();
auto minmax = std::minmax_element(reco_chi2.begin(), reco_chi2.end());
float chi2max = *minmax.second ;

if ( chi2max < CHI2MAX ) {
        if (debug) {
            std::cout << " --- DONE, all tracks have chi2 < CHI2MAX " << std::endl;
            std::cout  << "     number of primary tracks selected = " << seltracks.size() << std::endl;
        }
        return seltracks ;
}

if (debug) std::cout << " remove a track that has chi2 = " << chi2max << std::endl;

seltracks.erase( seltracks.begin() + maxElementIndex );
ipass ++;

 VertexingUtils::FCCAnalysesVertex vtx = VertexFitterSimple::VertexFitter_Tk(  isPrimaryVertex,
                                                                                seltracks,
                                                                         BeamSpotConstraint,
                                                                         bsc_sigmax, bsc_sigmay, bsc_sigmaz,
                                                                         bsc_x, bsc_y, bsc_z )  ;

 if(debug) std::cout << "Finished get_PrimaryTracks!" << std::endl;
 return VertexFitterSimple::get_PrimaryTracks( vtx, seltracks, BeamSpotConstraint, bsc_sigmax, bsc_sigmay, bsc_sigmaz,
                                                bsc_x,  bsc_y, bsc_z, ipass ) ;



}


ROOT::VecOps::RVec<edm4hep::TrackState>   VertexFitterSimple::get_NonPrimaryTracks( ROOT::VecOps::RVec<edm4hep::TrackState> allTracks,
                                                                                    ROOT::VecOps::RVec<edm4hep::TrackState> primaryTracks ) {
  if(debug) std::cout << "Starting get_NonPrimaryTracks!" << std::endl;
  ROOT::VecOps::RVec<edm4hep::TrackState> result;
  for (auto & track: allTracks) {
     bool isInPrimary = false;
     for ( auto &  primary:  primaryTracks) {
        if ( track.D0 == primary.D0 && track.Z0 == primary.Z0 &&  track.phi == primary.phi &&  track.omega == primary.omega && track.tanLambda == primary.tanLambda ) {
                isInPrimary = true;
                break;
        }
     }
     if ( !isInPrimary) result.push_back( track );
  }
 if(debug) std::cout << "Finished get_NonPrimaryTracks!" << std::endl;
 return result;
}


ROOT::VecOps::RVec<bool> VertexFitterSimple::IsPrimary_forTracks( ROOT::VecOps::RVec<edm4hep::TrackState> allTracks,
                                                                  ROOT::VecOps::RVec<edm4hep::TrackState> primaryTracks ) {

  if(debug) std::cout << "Starting IsPrimary_forTracks!" << std::endl;
  ROOT::VecOps::RVec<bool> result;
  for (auto & track: allTracks) {
     bool isInPrimary = false;
     for ( auto &  primary:  primaryTracks) {
        if ( track.D0 == primary.D0 && track.Z0 == primary.Z0 &&  track.phi == primary.phi &&  track.omega == primary.omega && track.tanLambda == primary.tanLambda ) {
                isInPrimary = true;
                break;
        }
     }
     result.push_back( isInPrimary );
  }
 if(debug) std::cout << "Finished IsPrimary_forTracks!" << std::endl;
 return result;
}



///////////////////////////
//** SV Finder (LCFI+) **//
///////////////////////////

VertexingUtils::FCCAnalysesSV VertexFitterSimple::get_SV_jets(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
							      ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
							      VertexingUtils::FCCAnalysesVertex PV,
							      ROOT::VecOps::RVec<bool> isInPrimary,
							      ROOT::VecOps::RVec<fastjet::PseudoJet> jets,
							      std::vector<std::vector<int>> jet_consti,
							      bool V0_rej,
							      double chi2_cut, double invM_cut, double chi2Tr_cut) {

  // find SVs using LCFI+ (clustering first)
  // change to vec of vec (RVec of RVec breaking) to separate SV from diff jets, currently don't separate SVs by jet
  
  VertexingUtils::FCCAnalysesSV SV;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
  SV.vtx = result;

  // find SV inside the jet loop (only from non-primary tracks)
  // first separate reco particles by jet then get the associated tracks
  int nJet = jets.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks;

  ROOT::VecOps::RVec<edm4hep::TrackState> tracks   = ReconstructedParticle2Track::getRP2TRK( recoparticles, thetracks );
  ROOT::VecOps::RVec<int> reco_ind_tracks     = ReconstructedParticle2Track::get_recoindTRK( recoparticles, thetracks );
  if(tracks.size() != reco_ind_tracks.size()) std::cout<<"ERROR: reco index vector not the same size as no of tracks"<<std::endl;

  if(tracks.size() != isInPrimary.size()) std::cout<<"ERROR: isInPrimary vector size not the same as no of tracks"<<std::endl;

  if(debug) std::cout<<"tracks extracted from the reco particles"<<std::endl;

  //
  for (unsigned int j=0; j<nJet; j++) {

    std::vector<int> i_jetconsti = jet_consti[j];
    for (int ctr=0; ctr<tracks.size(); ctr++) {
      if(isInPrimary[ctr]) continue; // remove primary tracks
      if(std::find(i_jetconsti.begin(), i_jetconsti.end(), reco_ind_tracks[ctr]) == i_jetconsti.end()) {
	np_tracks.push_back(tracks[ctr]); // separate tracks by jet
      }
    }
    
    if(debug) std::cout<<"primary tracks removed; there are "<<np_tracks.size()<<" non-primary tracks in jet#"<<j+1<<std::endl;
    
    // V0 rejection (tight)
    ROOT::VecOps::RVec<edm4hep::TrackState> tracks_fin;
    if(V0_rej) {
      bool tight = true;
      ROOT::VecOps::RVec<bool> isInV0 = isV0(np_tracks, PV, tight);
      for(unsigned int i=0; i<isInV0.size(); i++) {
	if (!isInV0[i]) tracks_fin.push_back(np_tracks[i]);
      }
    }
    else tracks_fin = np_tracks;
    
    if(debug) {
      std::cout<<np_tracks.size()-tracks_fin.size()<<" V0 tracks removed"<<std::endl;
      std::cout<<"now starting to find secondary vertices..."<<std::endl;
    }
    
    while(tracks_fin.size() > 1) {
      // find vertex seed
      ROOT::VecOps::RVec<int> vtx_seed = VertexSeed_best(tracks_fin, PV, chi2_cut, invM_cut);
      // constraint thresholds can be chosen by user, here using default cuts
      if(vtx_seed.size() == 0) break;
      
      // add tracks to the seed
      // check if a track is added; if not break loop
      ROOT::VecOps::RVec<int> vtx_fin = vtx_seed;
      int vtx_fin_size = 0; // to start the loop
      while(vtx_fin_size != vtx_fin.size()) {
	vtx_fin_size = vtx_fin.size();
	vtx_fin = addTrack_best(tracks_fin, vtx_fin, PV, chi2_cut, invM_cut, chi2Tr_cut);
	// constraint thresholds can be chosen by user, here using default cuts
      }
      
      // fit tracks to SV and remove from tracks_fin
      ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx_fin;
      for(int i_tr : vtx_fin) tr_vtx_fin.push_back(tracks_fin[i_tr]);
      VertexingUtils::FCCAnalysesVertex sec_vtx = VertexFitter_Tk(0, tr_vtx_fin);

      //fill the indices of the tracks
      sec_vtx.reco_ind = get_reco_ind(recoparticles,thetracks);

      result.push_back(sec_vtx);
      //
      ROOT::VecOps::RVec<edm4hep::TrackState> temp = tracks_fin;
      tracks_fin.clear();
      for(unsigned int t=0; t<temp.size(); t++) {
	if(std::find(vtx_fin.begin(), vtx_fin.end(), t) == vtx_fin.end()) tracks_fin.push_back(temp[t]);
      }
      // all this cause don't know how to remove multiple elements at once
      tr_vtx_fin.clear();
      
      if(debug) std::cout<<result.size()<<" SV found"<<std::endl;
    }
    
    // clean-up
    np_tracks.clear();
    tracks_fin.clear();
  }

  if(debug) std::cout<<"no more SVs can be reconstructed"<<std::endl;
  
  // currently don't know which SV is from which jet (FIX SOON)
  SV.vtx = result;
  //
  return SV;
}

VertexingUtils::FCCAnalysesSV VertexFitterSimple::get_SV_event(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
							       ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
							       VertexingUtils::FCCAnalysesVertex PV,
							       ROOT::VecOps::RVec<bool> isInPrimary,
							       bool V0_rej,
							       double chi2_cut, double invM_cut, double chi2Tr_cut) {


  if(debug) std::cout << "Starting SV finding!" << std::endl;

  // find SVs using LCFI+ (w/o clustering)
  // still need to think a little about jet clustering using SVs & pseudo-vertices as seeds
  
  VertexingUtils::FCCAnalysesSV SV;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
  SV.vtx = result;

  // retrieve the tracks associated to the recoparticles
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks = ReconstructedParticle2Track::getRP2TRK( recoparticles, thetracks );

  if(debug) std::cout<<"tracks extracted from the reco particles"<<std::endl;

  if(tracks.size() != isInPrimary.size()) std::cout<<"ISSUE: track vector and primary-nonprimary vector of diff sizes"<<std::endl;

  // find SV from non-primary tracks
  ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks;
  for(unsigned int i=0; i<tracks.size(); i++) {
    if (!isInPrimary[i]) np_tracks.push_back(tracks[i]);
  }

  if(debug) std::cout<<"primary tracks removed; there are "<<np_tracks.size()<<" non-primary tracks in the event"<<std::endl;

  // V0 rejection (tight)
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks_fin;
  if(V0_rej) {
    bool tight = true;
    ROOT::VecOps::RVec<bool> isInV0 = isV0(np_tracks, PV, tight);
    for(unsigned int i=0; i<isInV0.size(); i++) {
      if (!isInV0[i]) tracks_fin.push_back(np_tracks[i]);
    }
  }
  else tracks_fin = np_tracks;

  if(debug) {
    std::cout<<np_tracks.size()-tracks_fin.size()<<" V0 tracks removed"<<std::endl;
    std::cout<<"now starting to find secondary vertices..."<<std::endl;
  }
  

  if(debug) std::cout << "tracks_fin.size() = " << tracks_fin.size() << std::endl;
  while(tracks_fin.size() > 1) {
    // find vertex seed
    ROOT::VecOps::RVec<int> vtx_seed = VertexSeed_best(tracks_fin, PV, chi2_cut, invM_cut);
    
    if(debug){
        std::cout << "tracks_fin.size(): " << tracks_fin.size() << std::endl;
        for(int i=0; i<vtx_seed.size();i++)
	  std::cout << "vtx_seed: " << vtx_seed[i] << std::endl;
    }
    if(vtx_seed.size() == 0) break;
    
    // add tracks to the seed
    // check if a track is added; if not break loop
    ROOT::VecOps::RVec<int> vtx_fin = vtx_seed;
    int vtx_fin_size = 0; // to start the loop
    while(vtx_fin_size != vtx_fin.size()) {
      vtx_fin_size = vtx_fin.size();
      vtx_fin = addTrack_best(tracks_fin, vtx_fin, PV, chi2_cut, invM_cut, chi2Tr_cut);
    }
    
    // fit tracks to SV and remove from tracks_fin
    ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx_fin;
    for(int i_tr : vtx_fin){
      tr_vtx_fin.push_back(tracks_fin[i_tr]);
      if(debug) std::cout << "Pushing back tracks_fin[i_tr]" << std::endl;
    }
    VertexingUtils::FCCAnalysesVertex sec_vtx = VertexFitter_Tk(0, tr_vtx_fin);

    //fill the indices of the tracks
    sec_vtx.reco_ind = get_reco_ind(recoparticles,thetracks);

    result.push_back(sec_vtx);
    //
    ROOT::VecOps::RVec<edm4hep::TrackState> temp = tracks_fin;
    tracks_fin.clear();
    for(unsigned int t=0; t<temp.size(); t++) {
      if(std::find(vtx_fin.begin(), vtx_fin.end(), t) == vtx_fin.end()) tracks_fin.push_back(temp[t]);
    }
    // all this cause don't know how to remove multiple elements at once

    if(debug) std::cout<<result.size()<<" SV found"<<std::endl;
  }

  if(debug) std::cout<<"no more SVs can be reconstructed"<<std::endl;
  
  SV.vtx = result;
  //
  return SV;
}

VertexingUtils::FCCAnalysesSV VertexFitterSimple::get_SV_event(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
							       ROOT::VecOps::RVec<edm4hep::TrackState> thetracks,
							       ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
							       VertexingUtils::FCCAnalysesVertex PV,
							       bool V0_rej,
							       double chi2_cut, double invM_cut, double chi2Tr_cut) {
  
  // find SVs from non-primary tracks using LCFI+ (w/o clustering)
  // still need to think a little about jet clustering using SVs & pseudo-vertices as seeds
  // primary - non-primary separation done externally
  
  VertexingUtils::FCCAnalysesSV SV;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
  SV.vtx = result;

  // V0 rejection (tight)
  ROOT::VecOps::RVec<edm4hep::TrackState> tracks_fin;
  if(V0_rej) {
    bool tight = true;
    ROOT::VecOps::RVec<bool> isInV0 = isV0(np_tracks, PV, tight);
    for(unsigned int i=0; i<isInV0.size(); i++) {
      if (!isInV0[i]) tracks_fin.push_back(np_tracks[i]);
    }
  }
  else tracks_fin = np_tracks;

  if(debug) {
    std::cout<<np_tracks.size()-tracks_fin.size()<<" V0 tracks removed"<<std::endl;
    std::cout<<"now starting to find secondary vertices..."<<std::endl;
  }
  
  while(tracks_fin.size() > 1) {
    // find vertex seed
    ROOT::VecOps::RVec<int> vtx_seed = VertexSeed_best(tracks_fin, PV, chi2_cut, invM_cut);
    if(vtx_seed.size() == 0) break;
    
    // add tracks to the seed
    // check if a track is added; if not break loop
    ROOT::VecOps::RVec<int> vtx_fin = vtx_seed;
    int vtx_fin_size = 0; // to start the loop
    while(vtx_fin_size != vtx_fin.size()) {
      vtx_fin_size = vtx_fin.size();
      vtx_fin = addTrack_best(tracks_fin, vtx_fin, PV, chi2_cut, invM_cut, chi2Tr_cut);
    }
    
    // fit tracks to SV and remove from tracks_fin
    ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx_fin;
    for(int i_tr : vtx_fin) tr_vtx_fin.push_back(tracks_fin[i_tr]);
    VertexingUtils::FCCAnalysesVertex sec_vtx = VertexFitter_Tk(0, tr_vtx_fin);

    //
    ROOT::VecOps::RVec<edm4hep::TrackState> temp = tracks_fin;
    tracks_fin.clear();
    for(unsigned int t=0; t<temp.size(); t++) {
      if(std::find(vtx_fin.begin(), vtx_fin.end(), t) == vtx_fin.end()) tracks_fin.push_back(temp[t]);
    }
    // all this cause don't know how to remove multiple elements at once

    //fill the indices of the tracks
    sec_vtx.reco_ind = get_reco_ind(recoparticles,thetracks);
    
    result.push_back(sec_vtx);    

    if(debug) std::cout<<result.size()<<" SV found"<<std::endl;
  }

  if(debug) std::cout<<"no more SVs can be reconstructed"<<std::endl;
  
  SV.vtx = result;
  //
  return SV;
}


ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> VertexFitterSimple::get_all_vertices(VertexingUtils::FCCAnalysesVertex PV,
											   VertexingUtils::FCCAnalysesSV SV) {
  // Returns a vector of all vertices (PV and SVs)
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> result;
  result.push_back(PV);
  for (auto &p:SV.vtx){
    result.push_back(p);
  }
  return result;  
}



ROOT::VecOps::RVec<int> VertexFitterSimple::VertexSeed_best(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							    VertexingUtils::FCCAnalysesVertex PV,
							    double chi2_cut, double invM_cut) {

  // gives indices of the best pair of tracks
  // maybe also update and write one to get first pair to pass constraints

  ROOT::VecOps::RVec<int> result;
  int isel = 0;
  int jsel = 1;
  
  int nTr = tracks.size();
  // push empty tracks to make a size=2 vector
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  tr_pair.push_back(tracks[0]);
  tr_pair.push_back(tracks[1]);
  VertexingUtils::FCCAnalysesVertex vtx_seed;
  double chi2_min = 99;
  
  for(unsigned int i=0; i<nTr-1; i++) {
    tr_pair[0] = tracks[i];
    
    for(unsigned int j=i+1; j<nTr; j++) {
      tr_pair[1] = tracks[j];
      
      // V0 rejection (loose)
      ROOT::VecOps::RVec<bool> isInV0 = isV0(tr_pair, PV, false);
      if(isInV0[0] && isInV0[1]) continue;
      
      vtx_seed = VertexFitter_Tk(0, tr_pair);
      
      // Constraints
      // chi2 < cut (9)
      double chi2_seed = vtx_seed.vertex.chi2; // normalised
      if(chi2_seed >= chi2_cut) continue; // nDOF for 2 track vtx = 1
      //
      // invM < cut (10GeV)
      double invM_seed = VertexingUtils::get_invM(vtx_seed);
      if(invM_seed >= invM_cut) continue;
      //
      // invM < sum of energy
      double E_pair = 0.;
      for(edm4hep::TrackState tr_e : tr_pair) E_pair += VertexingUtils::get_trackE(tr_e);
      if(invM_seed >= E_pair) continue;
      //
      // momenta sum & vtx r on same side
      double angle = VertexingUtils::get_PV2vtx_angle(tr_pair, vtx_seed, PV);
      if(angle<0) continue;

      // if a pair passes all constraints compare chi2, store lowest chi2
      if(chi2_seed < chi2_min) {
	isel = i; jsel =j;
	chi2_min = chi2_seed;
      }
    }
  }

  result.push_back(isel); result.push_back(jsel);
  return result;
}

ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> VertexFitterSimple::VertexSeed_all(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
								 VertexingUtils::FCCAnalysesVertex PV,
								 double chi2_cut, double invM_cut) {

  // gives indices of the all pairs of tracks which pass the constraints

  ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>> result;
  ROOT::VecOps::RVec<int> ij_sel;
  
  int nTr = tracks.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  // push empty tracks to make a size=2 vector
  edm4hep::TrackState tr_i;
  edm4hep::TrackState tr_j;
  tr_pair.push_back(tr_i);
  tr_pair.push_back(tr_j);
  VertexingUtils::FCCAnalysesVertex vtx_seed;
  
  for(unsigned int i=0; i<nTr-1; i++) {
    if(i!=0) tr_pair[0] = tracks[i];

    for(unsigned int j=i+1; j<nTr; j++) {
      if(j!=1) tr_pair[1] = tracks[j];

      // V0 rejection (loose)
      ROOT::VecOps::RVec<bool> isInV0 = isV0(tr_pair, PV, false);
      if(isInV0[0] && isInV0[1]) continue;
      
      vtx_seed = VertexFitter_Tk(0, tr_pair);

      // Constraints
      // chi2 < cut (9)
      double chi2_seed = vtx_seed.vertex.chi2; // normalised
      if(chi2_seed >= chi2_cut) continue; // nDOF for 2 track vtx = 1
      //
      // invM < cut (10GeV)
      double invM_seed = VertexingUtils::get_invM(vtx_seed);
      if(invM_seed >= invM_cut) continue;
      //
      // invM < sum of energy
      double E_pair = 0.;
      for(edm4hep::TrackState tr_e : tr_pair) E_pair += VertexingUtils::get_trackE(tr_e);
      if(invM_seed >= E_pair) continue;
      //
      // momenta sum & vtx r on same side
      double angle = VertexingUtils::get_PV2vtx_angle(tr_pair, vtx_seed, PV);
      if(angle<0) continue;

      // if a pair passes all constraints, store indices
      ij_sel.push_back(i); ij_sel.push_back(j);
      result.push_back(ij_sel);
      ij_sel.clear();
    }
  }

  return result;
}

ROOT::VecOps::RVec<int> VertexFitterSimple::addTrack_best(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							  ROOT::VecOps::RVec<int> vtx_tr,
							  VertexingUtils::FCCAnalysesVertex PV,
							  double chi2_cut, double invM_cut, double chi2Tr_cut) {
  // adds index of the best track to the (seed) vtx
  
  ROOT::VecOps::RVec<int> result = vtx_tr;
  if(tracks.size() == vtx_tr.size()) return result;
  
  int isel = -1;

  int nTr = tracks.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx;
  VertexingUtils::FCCAnalysesVertex vtx;
  double chi2_min = 99;

  // add tracks of the previously formed vtx to a vector
  for(int tr : vtx_tr) {
    if(debug) std::cout << "Track integer: " << tr << std::endl;
    if(debug) std::cout <<  "Track value: " << tracks[tr] << std::endl;
    tr_vtx.push_back(tracks[tr]);
  }
  int iTr = tr_vtx.size();
  // add an empty track to increase vector size by 1
  tr_vtx.push_back(tracks[0]);

  // find best track to add to the vtx
  for(unsigned int i=0; i<nTr; i++) {
    if(std::find(vtx_tr.begin(), vtx_tr.end(), i) != vtx_tr.end()) continue;
    tr_vtx[iTr] = tracks[i];
    
    vtx = VertexFitter_Tk(0, tr_vtx);

    // Constraints
    // chi2_contribution(track) < threshold
    ROOT::VecOps::RVec<float> chi2_tr = vtx.reco_chi2;
    if(chi2_tr[iTr] >= chi2Tr_cut) continue; // threshold = 5 ok?
    //
    // chi2 < cut (9)
    double chi2_vtx = vtx.vertex.chi2; // normalised
    double nDOF = 2*(iTr+1) - 3; // nDOF = 2*nTr - 3
    chi2_vtx = chi2_vtx * nDOF;
    if(chi2_vtx >= chi2_cut) continue;
    //
    // invM < cut (10GeV)
    double invM_vtx = VertexingUtils::get_invM(vtx);
    if(invM_vtx >= invM_cut) continue;
    //
    // invM < sum of energy (should it be or not?)
    double E_vtx = 0.;
    for(edm4hep::TrackState tr_e : tr_vtx) E_vtx += VertexingUtils::get_trackE(tr_e);
    if(invM_vtx >= E_vtx) continue;
    //
    // momenta sum & vtx r on same side
    double angle = VertexingUtils::get_PV2vtx_angle(tr_vtx, vtx, PV);
    if(angle<0) continue;
    
    // if a track passes all constraints compare chi2, store lowest chi2
    if(chi2_vtx < chi2_min) {
      isel = i;
      chi2_min = chi2_vtx;
    }    
  }

  if(isel>=0) result.push_back(isel);
  return result;
}

ROOT::VecOps::RVec<int> VertexFitterSimple::addTrack_multi(ROOT::VecOps::RVec<edm4hep::TrackState> tracks,
							   ROOT::VecOps::RVec<int> vtx_tr,
							   VertexingUtils::FCCAnalysesVertex PV,
							   double chi2_cut, double invM_cut, double chi2Tr_cut) {
  // adds indices of all tracks passing constraints to the (seed) vtx
  
  ROOT::VecOps::RVec<int> result = vtx_tr;
  if(tracks.size() == vtx_tr.size()) return result;

  int nTr = tracks.size();
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_vtx;
  VertexingUtils::FCCAnalysesVertex vtx;

  // tracks from the previously formed vtx
  for(int tr : vtx_tr) {
    tr_vtx.push_back(tracks[tr]);
  }
  int iTr = tr_vtx.size();
  
  // find best track to add to the vtx
  for(unsigned int i=0; i<nTr; i++) {
    if(std::find(vtx_tr.begin(), vtx_tr.end(), i) != vtx_tr.end()) continue;

    if(iTr != tr_vtx.size()) tr_vtx[iTr] = tracks[i];
    else tr_vtx.push_back(tracks[i]);
    
    vtx = VertexFitter_Tk(0, tr_vtx);

    // Constraints
    // chi2_contribution < threshold.
    ROOT::VecOps::RVec<float> chi2_tr = vtx.reco_chi2;
    if(chi2_tr[iTr] >= chi2Tr_cut) continue; // threshold = 5 ok?
    //
    // chi2 < cut (9)
    double chi2_vtx = vtx.vertex.chi2; // normalised
    double nDOF = 2*(iTr+1) - 3; // nDOF = 2*nTr - 3
    chi2_vtx = chi2_vtx * nDOF;
    if(chi2_vtx >= chi2_cut) continue;
    //
    // invM < cut (10GeV)
    double invM_vtx = VertexingUtils::get_invM(vtx);
    if(invM_vtx >= invM_cut) continue;
    //
    // invM < sum of energy (should it be or not?)
    double E_vtx = 0.;
    for(edm4hep::TrackState tr_e : tr_vtx) E_vtx += VertexingUtils::get_trackE(tr_e);
    if(invM_vtx >= E_vtx) continue;
    //
    // momenta sum & vtx r on same side
    double angle = VertexingUtils::get_PV2vtx_angle(tr_vtx, vtx, PV);
    if(angle<0) continue;
    
    // if a pair passes all constraints add to the vtx
    result.push_back(i);
    iTr++;
  }

  return result;
}

ROOT::VecOps::RVec<bool> VertexFitterSimple::isV0(ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
						  VertexingUtils::FCCAnalysesVertex PV,
						  bool tight) {
  // V0 rejection
  // fn can be updated to reconstruct V0 and output its momentum, PID, etc instead

  // take all non-primary tracks & assign "true" to pairs that form V0
  // if(tight)  -> tight constraints
  // if(!tight) -> loose constraints

  int nTr = np_tracks.size();

  ROOT::VecOps::RVec<bool> result(nTr, false);
  // true -> forms a V0, false -> doesn't form a V0
  if(nTr<2) return result;
  
  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm  
  
  ROOT::VecOps::RVec<edm4hep::TrackState> t_pair;
  // push empty tracks to make a size=2 vector
  edm4hep::TrackState tr_i;
  edm4hep::TrackState tr_j;
  t_pair.push_back(tr_i);
  t_pair.push_back(tr_j);
  VertexingUtils::FCCAnalysesVertex V0;
  //
  const double m_pi = 0.13957039; // pi+- mass [GeV]
  const double m_p  = 0.93827208; // p+- mass
  const double m_e  = 0.00051099; // e+- mass
  //
  for(unsigned int i=0; i<nTr-1; i++) {
    if(result[i] == true) continue;
    t_pair[0] = np_tracks[i];

    for(unsigned int j=i+1; j<nTr; j++) {
      if(result[j] == true) continue;
      t_pair[1] = np_tracks[j];

      V0 = VertexFitter_Tk(0, t_pair);

      // invariant masses for V0 candidates
      double invM_Ks      = VertexingUtils::get_invM_pairs(V0, m_pi, m_pi);
      double invM_Lambda1 = VertexingUtils::get_invM_pairs(V0, m_pi, m_p);
      double invM_Lambda2 = VertexingUtils::get_invM_pairs(V0, m_p, m_pi);
      double invM_Gamma   = VertexingUtils::get_invM_pairs(V0, m_e, m_e);

      // V0 candidate distance from PV
      edm4hep::Vector3f r_V0 = V0.vertex.position; // in mm
      // does Vector3f class has similar functions as root vectors?
      TVector3 r_V0_PV(r_V0[0] - r_PV[0], r_V0[1] - r_PV[1], r_V0[2] - r_PV[2]);
      double r = r_V0_PV.Mag(); // in mm

      // angle b/n V0 candidate momentum & PV-V0 displacement vector
      double p_r = VertexingUtils::get_PV2V0angle(V0, PV);

      if(tight) {
	// Ks
	if(invM_Ks>0.493 && invM_Ks<0.503 && r>0.5 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	  break;
	}

	// Lambda0
	else if(invM_Lambda1>1.111 && invM_Lambda1<1.121 && r>0.5 && p_r>0.99995) {
	  result[i] = true;
	  result[j] = true;
	  break;
	}
	else if(invM_Lambda2>1.111 && invM_Lambda2<1.121 && r>0.5 && p_r>0.99995) {
	  result[i] = true;
	  result[j] = true;
	  break;
	}

	// photon conversion
	else if(invM_Gamma<0.005 && r>9 && p_r>0.99995) {
	  result[i] = true;
	  result[j] = true;
	  break;
	}	
      }

      else {
	// Ks
	if(invM_Ks>0.488 && invM_Ks<0.508 && r>0.3 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	}
	
	// Lambda0
	else if(invM_Lambda1>1.106 && invM_Lambda1<1.126 && r>0.3 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	}
	else if(invM_Lambda2>1.106 && invM_Lambda2<1.126 && r>0.3 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	}
	
	// photon conversion
	else if(invM_Gamma<0.01 && r>9 && p_r>0.999) {
	  result[i] = true;
	  result[j] = true;
	}	
      }
      
    }
  }

  return result;
}



ROOT::VecOps::RVec<int> VertexFitterSimple::get_reco_ind(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recoparticles,
							 ROOT::VecOps::RVec<edm4hep::TrackState> tracks){

 //fill the indices of the tracks
  ROOT::VecOps::RVec<int> reco_ind;
  int Ntr = tracks.size();
  for (auto & p: recoparticles) {
    //std::cout << " in VertexFitter:  a recoparticle with charge = " << p.charge << std::endl;
    if ( p.tracks_begin >=0 && p.tracks_begin<tracks.size()) {
      reco_ind.push_back( p.tracks_begin );
    }
  }
  if ( reco_ind.size() != Ntr ) std::cout << " ... problem in VertexFitterSimple::get_reco_ind, size of reco_ind != Ntr " << std::endl;
  return reco_ind;
}

///////////////////////////
//** V0 Reconstruction **//
///////////////////////////

VertexingUtils::FCCAnalysesV0 VertexFitterSimple::get_V0s(ROOT::VecOps::RVec<edm4hep::TrackState> np_tracks,
							  VertexingUtils::FCCAnalysesVertex PV,
							  bool tight,
							  double chi2_cut) {
  // V0 reconstruction
  // if(tight)  -> tight constraints
  // if(!tight) -> loose constraints

  // also look into how to reconstruct pi0 soon

  // make it stand-alone (removing primary tracks etc)

  VertexingUtils::FCCAnalysesV0 result;
  ROOT::VecOps::RVec<VertexingUtils::FCCAnalysesVertex> vtx; // FCCAnalyses vertex object
  ROOT::VecOps::RVec<int> pdgAbs;                            // absolute PDG ID
  ROOT::VecOps::RVec<double> invM;                           // invariant mass
  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;

  VertexingUtils::FCCAnalysesVertex V0_vtx;
  
  int nTr = np_tracks.size();
  if(nTr<2) return result;
  ROOT::VecOps::RVec<bool> isInV0(nTr, false);

  edm4hep::Vector3f r_PV = PV.vertex.position; // in mm  
  
  ROOT::VecOps::RVec<edm4hep::TrackState> tr_pair;
  // push empty tracks to make a size=2 vector
  edm4hep::TrackState tr_i;
  edm4hep::TrackState tr_j;
  tr_pair.push_back(tr_i);
  tr_pair.push_back(tr_j);
  //
  const double m_pi = 0.13957039; // pi+- mass [GeV]
  const double m_p  = 0.93827208; // p+- mass
  const double m_e  = 0.00051099; // e+- mass
  //
  for(unsigned int i=0; i<nTr-1; i++) {
    if(isInV0[i] == true) continue; // don't pair a track if it already forms a V0
    tr_pair[0] = np_tracks[i];

    for(unsigned int j=i+1; j<nTr; j++) {
      if(isInV0[j] == true) continue; // don't pair a track if it already forms a V0
      tr_pair[1] = np_tracks[j];

      V0_vtx = VertexFitter_Tk(0, tr_pair);

      // constraint on chi2: chi2 < cut (9)
      double chi2_V0 = V0_vtx.vertex.chi2; // normalised but DOF=1
      if(chi2_V0 >= chi2_cut) continue;

      // invariant masses for V0 candidates
      double invM_Ks      = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_pi);
      double invM_Lambda1 = VertexingUtils::get_invM_pairs(V0_vtx, m_pi, m_p);
      double invM_Lambda2 = VertexingUtils::get_invM_pairs(V0_vtx, m_p, m_pi);
      double invM_Gamma   = VertexingUtils::get_invM_pairs(V0_vtx, m_e, m_e);

      // V0 candidate distance from PV
      edm4hep::Vector3f r_V0 = V0_vtx.vertex.position; // in mm
      // does Vector3f class has similar functions as root vectors?
      TVector3 r_V0_PV(r_V0[0] - r_PV[0], r_V0[1] - r_PV[1], r_V0[2] - r_PV[2]);
      double r = r_V0_PV.Mag(); // in mm

      // angle b/n V0 candidate momentum & PV-V0 displacement vector
      double p_r = VertexingUtils::get_PV2V0angle(V0_vtx, PV);

      if(tight) {
	// Ks
	if(invM_Ks>0.493 && invM_Ks<0.503 && r>0.5 && p_r>0.999) {
	  if(debug) std::cout<<"Found a Ks"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(310);
	  invM.push_back(invM_Ks);
	  break;
	}
      
	// Lambda0
	else if(invM_Lambda1>1.111 && invM_Lambda1<1.121 && r>0.5 && p_r>0.99995) {
	  if(debug) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(invM_Lambda1);
	  break;
	}
	else if(invM_Lambda2>1.111 && invM_Lambda2<1.121 && r>0.5 && p_r>0.99995) {
	  if(debug) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(invM_Lambda2);
	  break;
	}
	
	// photon conversion
	else if(invM_Gamma<0.005 && r>9 && p_r>0.99995) {
	  if(debug) std::cout<<"Found a Photon coversion"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(22);
	  invM.push_back(invM_Gamma);
	  break;
	}
      }

      else {
	// Ks
	if(invM_Ks>0.488 && invM_Ks<0.508 && r>0.3 && p_r>0.999) {
	  if(debug) std::cout<<"Found a Ks"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(310);
	  invM.push_back(invM_Ks);
	  break;
	}
      
	// Lambda0
	else if(invM_Lambda1>1.106 && invM_Lambda1<1.126 && r>0.3 && p_r>0.999) {
	  if(debug) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(invM_Lambda1);
	  break;
	}
	else if(invM_Lambda2>1.106 && invM_Lambda2<1.126 && r>0.3 && p_r>0.999) {
	  if(debug) std::cout<<"Found a Lambda0"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(3122);
	  invM.push_back(invM_Lambda2);
	  break;
	}
	
	// photon conversion
	else if(invM_Gamma<0.01 && r>9 && p_r>0.999) {
	  if(debug) std::cout<<"Found a Photon coversion"<<std::endl;
	  isInV0[i] = true;
	  isInV0[j] = true;
	  vtx.push_back(V0_vtx);
	  pdgAbs.push_back(22);
	  invM.push_back(invM_Gamma);
	  break;
	}
      }
      
    }
  }

  //std::cout<<"Found "<<vtx.size()<<" V0s"<<std::endl;
  
  result.vtx = vtx;
  result.pdgAbs = pdgAbs;
  result.invM = invM;
  //
  return result;
}
