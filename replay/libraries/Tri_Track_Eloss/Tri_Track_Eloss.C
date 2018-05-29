//--- Adapted from the Gmp equivalent written by Barak Schmookler & Thir Gautam  12-01-15
//--- by Reynier Cruz Torres for the Tritium experiments.
//--- Revision History: 
//    1) Barak Schmookler, August 2016: Change from loop 2 to loop 3. Loop
//       3 was used in Spring 2016.
//    2) Longwu Ou, October 2016: Now the target geometry is determined by a
//       tag "<fprefix>.config" in the run database.
//    3) Reynier Cruz Torres, May 2018: Adapted code for the specifics
//       of the Tritium gas targets.
//
//////////////////////////////////////////////////////////////////////////
//
// Tri_Track_Eloss
//
// Correct energy loss of a track. For the gas target the
// correction will be vertex and angle dependent. Also, all the material
// after the target and before the spectrometer vacuum is included.
//
// The user can set a target Z offset while instantiating this class.
// However, this offset has a different meaning in the case of the gas
// target and the other targets. We also allow the user to set the length
// of the target (default is 25 cm). Since the amount of air between the
// scattering chamber and the spectrometer entrance differs for the LHRS
// and the RHRS, this can be set by the user.
//
// It is important to note that we assume the track originates from 
// (x,y) = (0,0). Also, we currently only include the material between
// the target and the spectrometer vacuum when using the gas target.
//////////////////////////////////////////////////////////////////////////

#include "Tri_Track_Eloss.h"
#include "THaSpectrometer.h"
#include "THaTrack.h"
#include "THaTrackInfo.h"
#include "THaVertexModule.h"
#include "TMath.h"
#include "TVector3.h"
#include "Riostream.h"

using namespace std;

//_____________________________________________________________________________
Tri_Track_Eloss::Tri_Track_Eloss( const char* name, 
		const char* description,
		const char* input_tracks,
		const char* vertex,
		const Double_t t_length,
		const Double_t z_offset,
		const Double_t air_length,
		Double_t particle_mass,
		Int_t hadron_charge ) :
	THaElossCorrection(name,description,input_tracks,particle_mass,
			hadron_charge), 
	fTLength(t_length), 
	fAirLength(air_length),
	fTrackModule(NULL)
{
	// Normal constructor.

	// Sets fExtPathMode = kTRUE, which allows the vertex 
	// module to be initialized. Also sets zref and the
	// scale factor.
	Double_t scale = 1.0 ; //Set scale factor to 1
	SetPathlength(vertex,z_offset,scale);

	Clear();
}

//_____________________________________________________________________________
Tri_Track_Eloss::~Tri_Track_Eloss()
{
	// Destructor

	DefineVariables( kDelete );
}

//_____________________________________________________________________________
void Tri_Track_Eloss::CalcEloss( THaTrackInfo* trkifo )
{
	// Compute the expected energy loss for the track given in trkifo
	// as well as the vertex location.
	//

	TVector3 MomDirectionHCS;
	Double_t vx =fVertexModule->GetVertex().X();
	Double_t vy =fVertexModule->GetVertex().Y();
	Double_t vz =fVertexModule->GetVertex().Z(); // vertex position
	Double_t thP=trkifo->GetPvect().X() / trkifo->GetPvect().Z();// thP=dx/dz
	Double_t phP=trkifo->GetPvect().Y() / trkifo->GetPvect().Z();// phP=dy/dz
	Double_t px=trkifo->GetPvect().X();
	Double_t py=trkifo->GetPvect().Y();
	Double_t pz=trkifo->GetPvect().Z();  
	Double_t p0 = trkifo->GetP();
	Double_t beta = p0 / TMath::Sqrt(p0*p0 + fM*fM);

	MomDirectionHCS(0)=px;
	MomDirectionHCS(1)=py;
	MomDirectionHCS(2)=pz;

	//Set Material Parameter

	//Aluminum
	// It's actually Al 7075 (rho = 2.81 g/cc) for the target,
	// and Al 2024-T3 (rho = 2.78 g/cc) for the scattering chamber
	// exit. But, just use 2.80 g/cc here...
	Double_t Z_Al   = 13. ;
	Double_t A_Al   = 26.98 ;
	Double_t rho_Al = 2.80 ;

	//Air
	Double_t Z_Air   = 7.22000 ;
	Double_t A_Air   = 14.46343 ;
	Double_t rho_Air = 1.20480E-03 ;

	//Kapton
	Double_t Z_Kap   = 5.02564 ;
	Double_t A_Kap   = 9.80345 ;
	Double_t rho_Kap = 1.42000 ;

	// Check if target corresponds to one of the gas cells
	if(fZmed < 3. && fAmed < 4.){	

		Double_t eloss_gas(0), eloss_Al1(0), eloss_Al2(0), eloss_Air(0), eloss_Kap(0) ;
		Double_t l_gas(0), l_Al1(0), l_Al2(0),l_Air,l_Kap(0) ;

		//Set Target Geometry
		Double_t R; // diameters is 0.77 inch
		Double_t Tip_Al = 0.11*1.e-3; //Tip thickness(m)

		// ---------------------------------------------------------------------------------------------
		// Thickness of the target side walls
		Double_t T_Al = 0; // Aluminum Wall Thickness (m) 

		if     ( fZmed>0.09 && fZmed<1.01 && fAmed>0.09 && fAmed<1.01 ) T_Al=(0.374+0.503)*1E-3/2.; // m - Hydrogen
		else if( fZmed>0.09 && fZmed<1.01 && fAmed>1.09 && fAmed<2.01 ) T_Al=(0.447+0.371)*1E-3/2.; // m - Deuterium
		else if( fZmed>0.09 && fZmed<1.01 && fAmed>2.09 && fAmed<3.01 ) T_Al=(0.435+0.447)*1E-3/2.; // m - Tritium
		else if( fZmed>1.09 && fZmed<2.01 && fAmed>2.09 && fAmed<3.01 ) T_Al=(0.487+0.478)*1E-3/2.; // m - Helium-3
		else T_Al=(0.374+0.503+0.447+0.371+0.435+0.447+0.487+0.478)*1E-3/8.; // m - Average of the numbers above

		// These numbers come from the following document:
		// https://logbooks.jlab.org/files/2018/01/3514205/TGT-RPT-17-007.pdf
		// for each target I average the "mid left" and "mid right" values
		// ---------------------------------------------------------------------------------------------

		R = 0.77/2.*0.0254; // 0.77/2 inches to meters
		// (got this radius from doc: TGT-103-1000-0101.pdf)

		Double_t Z_cap = (fTLength/2. - R) + fZref;
		Double_t pi = 3.14159265358979323846;
		// Double_t Z_rel = Z - Z_cap ;

		Double_t tx=0,ty=0,tz=Z_cap;//center postion of the tip
		TVector3 v_tipcenter; //vector from vertex to tip center, the direction is important
		Double_t distance_tip=0;//distance between the vertex and the tip center
		Double_t angleR=0,sin_angleR,cos_angleR=0; 

		//angleR is the angle between the scattering momentum and the vertex-tipcenter line
		distance_tip = TMath::Sqrt(TMath::Power(vx-tx,2)+TMath::Power(vy-ty,2)+TMath::Power(vz-tz-fZref,2)); 
		v_tipcenter(0)=tx-vx;
		v_tipcenter(1)=ty-vy;
		v_tipcenter(2)=tz-vz;
		cos_angleR =  (v_tipcenter.Dot(MomDirectionHCS))/( v_tipcenter.Mag()*MomDirectionHCS.Mag());
		sin_angleR=TMath::Sqrt(1-cos_angleR*cos_angleR);
		angleR = TMath::ACos( cos_angleR );

		Double_t anglex=0,sin_anglex=0;//anglex corresponds to distance_tip
		Double_t D=0,angleD=0;         //angleD corresponds to the traveling length
		Double_t x,y,z;                //the calculated position at the cylindral
		Double_t m1,m2,a,b,c,result1,result2;
		Double_t sin_theta=0;

		a=(px*px/py/py+1);
		m1=vx/R;
		m2=px/py/R*vy;
		b=2*px/py*(m1-m2);
		c=(m1-m2)*(m1-m2)-1;
		result1 = (-b+TMath::Sqrt(b*b-4*a*c))/2/a; 
		result2=(-b-TMath::Sqrt(b*b-4*a*c))/2/a;

		sin_theta = (R*result1-vy)*py>=0?result1:sin_theta ;
		sin_theta = (R*result2-vy)*py>=0?result2:sin_theta ;

		y=R*sin_theta;
		x=px/py*y-px/py*vy+vx;
		z=pz/py*(y-vy)+vz;

		if( z< Z_cap )
		{
			D = TMath::Sqrt((x-vx)*(x-vx)+(y-vy)*(y-vy)+(z-vz)*(z-vz));
			l_gas = D;
			l_Al1 = T_Al /TMath::Sqrt(thP*thP+phP*phP)*TMath::Sqrt(thP*thP+phP*phP+1) ; 
		}
		else
		{
			sin_anglex=distance_tip/R*sin(angleR);
			anglex = TMath::ASin(sin_anglex);
			angleD = pi - angleR - anglex;
			D=R/sin(angleR)*sin(angleD);

			l_gas = D;
			l_Al1 = Tip_Al /cos(anglex) ; 
		}

		// Aluminum Scattering Chamber Exit Window
		l_Al2 = 4.064E-4 ; // (i.e. 0.016")

		// Air between Scattering Chamber and Spectrometer
		l_Air = fAirLength ;
		// 13.95" (RHRS), 10.62" (LHRS) for Spring 16
		// 14.79" (RHRS), 15.23" (LHRS) for Fall 16		

		// Kapton window at Spectrometer Entrance
		l_Kap = 3.048E-4 ; // (i.e. 0.012")

		// --------------------------------------------------------------
		// Calculate energy loss with the parameters determined above

		if( fElectronMode ) {
			eloss_gas = ElossElectron( beta, fZmed, fAmed, fDensity, l_gas ); // Gas Target
			eloss_Al1 = ElossElectron( beta, Z_Al , A_Al , rho_Al  , l_Al1 ); // Aluminum Target Wall 
			eloss_Al2 = ElossElectron( beta, Z_Al , A_Al , rho_Al  , l_Al2 ); // Aluminum Scattering Chamber Exit Window
			eloss_Air = ElossElectron( beta, Z_Air, A_Air, rho_Air , l_Air ); // Air between Scattering Chamber and HRS
			eloss_Kap = ElossElectron( beta, Z_Kap, A_Kap, rho_Kap , l_Kap ); // Kapton window at Spectrometer Entrance
		}
		else{
			eloss_gas = ElossHadron( fZ, beta, fZmed, fAmed, fDensity, l_gas ); // Gas Target
			eloss_Al1 = ElossHadron( fZ, beta, Z_Al , A_Al , rho_Al  , l_Al1 ); // Aluminum Target Wall 
			eloss_Al2 = ElossHadron( fZ, beta, Z_Al , A_Al , rho_Al  , l_Al2 ); // Aluminum Scattering Chamber Exit Window
			eloss_Air = ElossHadron( fZ, beta, Z_Air, A_Air, rho_Air , l_Air ); // Air between Scattering Chamber and HRS
			eloss_Kap = ElossHadron( fZ, beta, Z_Kap, A_Kap, rho_Kap , l_Kap ); // Kapton window at Spectrometer Entrance
		}
		// Calculate Total Eloss
		fPathlength = l_gas; // Set as pathlength through gas...not really important
		fEloss = eloss_gas + eloss_Al1 + eloss_Al2 + eloss_Air + eloss_Kap ;

	}
	// if this is not one of the gas cells
	else{
		if( fElectronMode ) {
			fEloss = ElossElectron( beta, fZmed, fAmed, fDensity, fPathlength );
		}
		else {
			fEloss = ElossHadron( fZ, beta, fZmed, fAmed, fDensity, fPathlength );
		}
	}

}

//_____________________________________________________________________________
void Tri_Track_Eloss::Clear( Option_t* opt )
{
	// Clear all event-by-event variables.

	THaElossCorrection::Clear(opt);
	TrkIfoClear();
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus Tri_Track_Eloss::Init( const TDatime& run_time )
{
	// Initialize the module.
	// Locate the input tracking module named in fInputName and save
	// pointer to it.

	static const char* const here = "Init()";

	// Find the input tracking module
	fTrackModule = dynamic_cast<THaTrackingModule*>
		( FindModule( fInputName.Data(), "THaTrackingModule"));
	if( !fTrackModule )
		return fStatus;

	// Get the parent spectrometer apparatus from the input module  
	THaSpectrometer* spectro = fTrackModule->GetTrackInfo()->GetSpectrometer();
	if( !spectro ) {
		Error( Here(here), "Oops. Input tracking module has no pointer "
				"to a spectrometer?!?" );
		return fStatus = kInitError;
	}
	// Needed for initialization of dependent modules in a chain
	fTrkIfo.SetSpectrometer( spectro );

	// Standard initialization. Calls this object's DefineVariables() and
	// reads meterial properties from the run database.
	THaElossCorrection::Init( run_time );

	return fStatus;
}

//_____________________________________________________________________________
Int_t Tri_Track_Eloss::DefineVariables( EMode mode )
{
	// Define/delete global variables.

	if( mode == kDefine && fIsSetup ) return kOK;
	THaElossCorrection::DefineVariables( mode );

	return DefineVarsFromList( THaTrackingModule::GetRVarDef(), mode );
}

//_____________________________________________________________________________
Int_t Tri_Track_Eloss::Process( const THaEvData& )
{
	// Calculate corrections and adjust the track parameters.

	if( !IsOK() ) return -1;

	THaTrackInfo* trkifo = fTrackModule->GetTrackInfo();
	if( !trkifo->IsOK() ) return 2;

	// Copy the input track info
	fTrkIfo = *trkifo;

	// Compute the correction
	Double_t p_out = fTrkIfo.GetP(); 
	if( p_out <= 0.0 ) return 4; //oops
	Double_t E_out = TMath::Sqrt(p_out*p_out + fM*fM);

	if( !fTestMode ) {
		// calculate pathlength for this event if we have a vertex module
		// For electron tracks, this is only used for solid target. It may
		// not be very useful because there is no angle dependence.
		// Also note that fZref here means the vertex value where the
		// pathlength is zero; while for the gas target, fZref means
		// how much the target is shifted along Z
		if( fExtPathMode ) {
			if( !fVertexModule->HasVertex() )
				return 1;
			fPathlength = 
				TMath::Abs(fVertexModule->GetVertex().Z() - fZref) * fScale;
		}
		//update fEloss
		CalcEloss(trkifo); 
	}
	Double_t p_in = TMath::Sqrt(p_out*p_out + fEloss*fEloss + 2.0*E_out*fEloss);

	// Apply the correction
	fTrkIfo.SetP(p_in);

	fDataValid = true;
	return 0;
}
