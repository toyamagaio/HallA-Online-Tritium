// =========================================================
// Coincidence Replay Script
// =========================================================
// Before changing anything in this code, please, notify:
// Florian Hauenstein:  hauenst@jlab.org
// Reynier Cruz-Torres: reynier@mit.edu

#include "def_tritium.h"
using namespace std;

void replay_coinc_new(Int_t runnumber=0,Int_t numevents=0,Int_t fstEvt=0,Bool_t QuietRun = kFALSE, Bool_t OnlineReplay =kFALSE, Bool_t bPlots = kFALSE, Bool_t autoreplay = kFALSE){

	char buf[300];
	Int_t nrun=0;
	if (runnumber<=0)
	{
		cout << "\nreplay: Please enter a Run Number (-1 to exit):";
		cin >> nrun;
		fgets(buf,300,stdin);//get the extra '\n' from stdin
		if( nrun<=0 ) return;
		runnumber = nrun;
	}
	// ---------------------------
	//Enable modules
	Bool_t bScaler   = kTRUE  ;
	Bool_t bHelicity = kFALSE ;
	Bool_t bBeam     = kTRUE  ;
	Bool_t bPhysics  = kTRUE  ;
	Bool_t bOldTrack = kFALSE ;
	// ---------------------------
	TString rootname;
	if(OnlineReplay) rootname = "%s/tritium_online_%d.root";
	else             rootname = "%s/tritium_%d.root"       ;
	// ---------------------------
	const char* RNAME=rootname.Data();
	TString ODEF = "coinc.odef";
	TString CUTS = "coinc.cuts";
	if(autoreplay)  ODEF=Form(REPLAY_DIR_PREFIX,"coinc_auto.odef");

	//=====================================================================================================================
	// Detectors
	//==================================
	// LHRS Detectors
	//==================================	
	//THaHRS *HRSL = new THaHRS("L","Left arm HRS"); //Add vdc,s2...uses s0 for track beta
	Tritium_HRS* HRSL = new Tritium_HRS("L","Left arm HRS");
	HRSL->AutoStandardDetectors(kFALSE);
	gHaApps->Add( HRSL );
	HRSL->AddDetector( new TriFadcXscin    ("s0"  , "s0 scintillator",kFALSE )); // boolean is about the fact that pmts A and B are flipped in the HRSs
	HRSL->AddDetector( new THaVDC          ("vdc" , "Vertical Drift Chamber" ));
	HRSL->AddDetector( new TriFadcCherenkov("cer" , "Gas Cherenkov counter"  ));
	HRSL->AddDetector( new TriFadcScin     ("s2"  , "S2 Scintillator"        ));
	HRSL->AddDetector( new THaShower       ("prl1", "Pre-shower pion rej."   ));
	HRSL->AddDetector( new THaShower       ("prl2", "Show pion rej."         )); 

	THaHRS* FbusHRSL = new THaHRS("FbusL", "Fastbus LHRS Readout");
	FbusHRSL->AutoStandardDetectors(kFALSE);
	gHaApps->Add(FbusHRSL);
	FbusHRSL->AddDetector( new THaCherenkov   ("cer", "Gas Cherenkov counter - Fastbus"  ));
	FbusHRSL->AddDetector( new THaScintillator("s2" , "S2 Scintillator - Fastbus"        ));
	FbusHRSL->AddDetector( new Tritium_Xscin  ("s0" , "S0 Scintillator - Fastbus", kFALSE));

	//==================================
	// RHRS Detectors
	//==================================
	// THaHRS* HRSR = new THaHRS("R","Right arm HRS");
	Tritium_HRS* HRSR = new Tritium_HRS("R","Right arm HRS");
	HRSR->AutoStandardDetectors(kFALSE);
	gHaApps->Add( HRSR );
	HRSR->AddDetector( new TriFadcXscin    ("s0" , "s0 scintillator",kTRUE        ));
	HRSR->AddDetector( new THaVDC          ("vdc", "Vertical Drift Chamber"       ));
	HRSR->AddDetector( new TriFadcCherenkov("cer", "Gas Cherenkov counter - FADC" ));
	HRSR->AddDetector( new TriFadcScin     ("s2" , "S2 Scintillator - FADC"       ));
	HRSR->AddDetector( new THaShower       ("ps" , "Pre-shower pion rej."         ));
	HRSR->AddDetector( new THaShower       ("sh" , "Show pion rej."               ));

	THaHRS* FbusHRSR = new THaHRS("FbusR", "Fastbus RHRS Readout");
	FbusHRSR->AutoStandardDetectors(kFALSE);
	gHaApps->Add(FbusHRSR);
	FbusHRSR->AddDetector( new THaCherenkov   ("cer", "Gas Cherenkov counter - Fastbus" ));
	//FbusHRSR->AddDetector( new THaScintillator("s2" , "S2 Scintillator - Fastbus"       ));
	FbusHRSR->AddDetector( new Tritium_Xscin  ("s0" , "S0 Scintillator - Fastbus", kTRUE));

	//==================================
	if(bOldTrack){
		THaApparatus* OldTrackL = new TriHRS("OldTrackL","old analyzer track");
		OldTrackL->AddDetector( new TriXscin("s0","s0 sintillator",kFALSE) );
		gHaApps->Add( OldTrackL );

		THaApparatus* OldTrackR = new TriHRS("OldTrackR","old analyzer track");
		OldTrackR->AddDetector( new TriXscin("s0","s0 sintillator",kTRUE) );
		gHaApps->Add( OldTrackR );
	}

	//=====================================================================================================================
	// Scalers
	//==================================
	if(bScaler){	
		//==================================
		// LHRS Scalers
		//==================================
		TriScalerEvtHandler* lscaler = new TriScalerEvtHandler("Left","HA scaler event type 140 on L-HRS");
		gHaEvtHandlers->Add(lscaler);

		Tritium_TSScaler* levscaler = new Tritium_TSScaler("evLeft","HA scaler event type 1-14 on L-HRS");
		gHaEvtHandlers->Add(levscaler);

		Tritium_THaScaler100EvtHandler* lEndscaler = new Tritium_THaScaler100EvtHandler("EndLeft","HA scaler event type 100");
		gHaEvtHandlers->Add(lEndscaler);

		// Marco - for F1 tdc:
		gHaEvtHandlers->Add (new TdcDataEvtHandler("LTDC","F1 TDCs lHRS")); // do not change the "LTDC" word
		// Evan - V1495 Clock Counter:
		gHaEvtHandlers->Add (new ClockCountEvtHandler("LV1495","V1495 LHRS"));
		
		//==================================
		// RHRS Scalers
		//==================================
		TriScalerEvtHandler* rscaler = new TriScalerEvtHandler("Right","HA scaler event type 140 on R-HRS");
		gHaEvtHandlers->Add(rscaler);

		Tritium_TSScaler* revscaler = new Tritium_TSScaler("evRight","HA scaler event type 1-14 on R-HRS");
		gHaEvtHandlers->Add(revscaler);

		Tritium_THaScaler100EvtHandler* rEndscaler = new Tritium_THaScaler100EvtHandler("EndRight","HA scaler event type 100");
		gHaEvtHandlers->Add(rEndscaler);
		// Marco - F1 and VETROC tdcs:
		gHaEvtHandlers->Add (new TdcDataEvtHandler("RTDC","F1 and VETROC TDCs rHRS")); // do not change the "RTDC" word
		// Evan - V1495 Clock Counter:
		gHaEvtHandlers->Add (new ClockCountEvtHandler("RV1495","V1495 RHRS"));
	}

	//=====================================================================================================================
	// Decoder Data
	//==================================
	// RHRS
	THaApparatus* decR = new THaDecData("DR","Misc. Decoder Data");
	gHaApps->Add( decR );
	//==================================
	// LHRS
	THaApparatus* decL = new THaDecData("DL","Misc. Decoder Data");
	gHaApps->Add( decL );

	//=====================================================================================================================
	// Beam
	//==================================
	if(bBeam){	

		THaIdealBeam* ib = new THaIdealBeam("ib","Ideal beam");
                gHaApps->Add(ib);

		//==================================
		// LHRS Beam
		//==================================

		TriFadcRasteredBeam* Lrb = new TriFadcRasteredBeam("Lrb", "Rastered beam to L-HRS");
		gHaApps->Add(Lrb);

		THaRasteredBeam* FbusLrb = new THaRasteredBeam("FbusLrb", "Fastbus Rastered beam to L-HRS");
		FbusLrb->AddDetector(new THaRaster("Raster2", "Downstream Raster"));
		FbusLrb->AddDetector(new THaBPM   ("BPMA", "First BPM"));
		FbusLrb->AddDetector(new THaBPM   ("BPMB", "Second BPM"));
		gHaApps->Add(FbusLrb);
		
		//==================================
		// RHRS Beam
		//==================================	

		TriFadcRasteredBeam* Rrb = new TriFadcRasteredBeam("Rrb", "Rastered beam to the R-HRS");
		gHaApps->Add(Rrb);

		THaRasteredBeam* FbusRrb = new THaRasteredBeam("FbusRrb", "Fastbus Rastered beam to R-HRS");
		FbusRrb->AddDetector(new THaRaster("Raster2", "Downstream Raster"));
		FbusRrb->AddDetector(new THaBPM   ("BPMA", "First BPM"));
		FbusRrb->AddDetector(new THaBPM   ("BPMB", "Second BPM"));
		gHaApps->Add(FbusRrb);
	}

	//=====================================================================================================================
	// Physics
	//==================================
	if(bPhysics){
		Double_t mass_el  = 0.000511      ; // electron mass
		Double_t amu      = 931.494*1.e-3 ; // amu to GeV
		Double_t mass_He3 = 3.0160293 *amu;
		Double_t mass_H2  = 2.01410178*amu;
		Double_t mass_H3  = 3.0160492 *amu;
		Double_t mass_tg  = 0.938         ; //default target 
		Double_t mass_prot= 0.938         ; //GeV proton mass

		string word[5],line;

		// File to look up the target enconder position to determine the target for this run
		TString filename = Form("/adaqfs/home/adaq/epics/runfiles_tritium_Coinc/Start_of_Run_%d.epics",runnumber);

		ifstream infile;
		infile.open(filename);
		double pos=0       ;
		double t2 =33106235;
		double d2 =29367355;
		double he3=21875520;

		while(getline(infile,line)){
			istringstream str(line);
			str>>word[0]>>word[1];
			if(word[0]=="Encoder" && word[1]=="Position"){
				str>>word[2]>>word[3];
				pos = atof(word[3].c_str()); 
				if     (abs(pos-t2 )<50) mass_tg=mass_H3 /3.0; // "proton mass"
				else if(abs(pos-d2 )<50) mass_tg=mass_H2 /2.0;
				else if(abs(pos-he3)<50) mass_tg=mass_He3/3.0;
				break;
			}
		}

		//==================================
		// LHRS
		//==================================
		THaPhysicsModule *Lgold   = new THaGoldenTrack  ("L.gold"   ,"HRS-L Golden Track"                     ,"L"                 );
		THaPhysicsModule *Lvdceff = new TriVDCeff       ("L.vdceff" ,"Left vdc efficiency"                                         );
		THaPhysicsModule *EKL     = new THaPrimaryKine  ("EKL"      ,"Electron kinematics in HRS-L"	      ,"L"   ,"ib" ,mass_tg);
		THaPhysicsModule *EKLc    = new THaPrimaryKine  ("EKLc"     ,"Corrected Electron kinematics in HRS-L" ,"L"   ,"Lrb",mass_tg);
		THaReactionPoint *rpl     = new THaReactionPoint("rpl"      ,"Reaction vertex for HRS-L"              ,"L"   ,"Lrb"        );
		THaExtTarCor     *exL     = new THaExtTarCor    ("exL"      ,"Corrected for extended target, HRS-L"   ,"L"   ,"rpl"        );
		THaPhysicsModule *EKLx    = new THaPrimaryKine  ("EKLx"     ,"Better Corrected Electron kinem in LHRS","exL" ,"Lrb",mass_tg);
		THaPhysicsModule *BCML    = new TriBCM          ("LeftBCM"  ,"Beam Current Monitors"                  ,"Left",""   ,0      );
		THaPhysicsModule *BCMevL  = new TriBCM          ("LeftBCMev","Beam Current Monitors"                  ,"Left","ev" ,0      );

		gHaPhysics->Add(Lgold  );
		gHaPhysics->Add(Lvdceff);
		gHaPhysics->Add(EKL    );
		gHaPhysics->Add(EKLc   );
		gHaPhysics->Add(rpl    );
		gHaPhysics->Add(exL    );
		gHaPhysics->Add(EKLx   );
		gHaPhysics->Add(BCML   );
		gHaPhysics->Add(BCMevL );

		//==================================
		// RHRS
		//==================================
		THaPhysicsModule *Rgold   = new THaGoldenTrack  ("R.gold"    ,"HRS-R Golden Track"                     ,"R"                     );
		THaPhysicsModule *Rvdceff = new TriVDCeff       ("R.vdceff"  ,"Reft vdc efficiency"                                             );
		THaPhysicsModule *EKR     = new THaSecondaryKine("EKR"       ,"Proton kinematics in HRS-R"             ,"R"    ,"EKL" ,mass_prot);
		THaPhysicsModule *EKRc    = new THaSecondaryKine("EKRc"      ,"Corrected Proton kinematics in HRS-R"   ,"R"    ,"EKLc",mass_prot);
		THaReactionPoint *rpr     = new THaReactionPoint("rpr"       ,"Reaction vertex for HRS-R"              ,"R"    ,"Lrb"           );
		THaExtTarCor     *exR     = new THaExtTarCor    ("exR"       ,"Corrected for extended target, HRS-R"   ,"R"    ,"rpr"           );
		THaPhysicsModule *EKRx    = new THaSecondaryKine("EKRx"      ,"Better Corrected Proton kinem in RHRS"  ,"R"    ,"EKLx",mass_prot);
		THaPhysicsModule *BCMR    = new TriBCM          ("RightBCM"  ,"Beam Current Monitors"                  ,"Right",""    ,0        );
		THaPhysicsModule *BCMevR  = new TriBCM          ("RightBCMev","Beam Current Monitors"                  ,"Right","ev"  ,0        );

		gHaPhysics->Add(Rgold  );
		gHaPhysics->Add(Rvdceff);
		gHaPhysics->Add(EKR    );
		gHaPhysics->Add(EKRc   );
		gHaPhysics->Add(rpr    );
		gHaPhysics->Add(exR    );
		gHaPhysics->Add(EKRx   );
		gHaPhysics->Add(BCMR   );
		gHaPhysics->Add(BCMevR );

	}

	//=====================================================================================================================
	// Set up Analyzer and replay data
	//=====================================
	ReplayCore(
			runnumber,      // run #
			numevents,      // -1=replay all;0=ask for a number
			50000,          // default replay event num
			RNAME,          // output file format
			ODEF.Data(),	// out define
			CUTS.Data(),	// empty cut define
			bScaler,	// replay scalar?
			bHelicity,	// repaly helicity
			fstEvt,		// First Event To Replay
			QuietRun	// whether ask user for inputs
		  );

	//=====================================================================================================================
	// Generate online plots
	//=====================================
	if(bPlots){

		const char* path_to_plots = "/chafs1/work1/tritium/Run_pdfs/";

		const char* GUI_DIR = Form(REPLAY_DIR_PREFIX,"onlineGUI64/");
		const char* SUM_DIR = Form(REPLAY_DIR_PREFIX,"summaryfiles/");

		//==================================
		// LHRS
		//==================================
 		const char * config_online_L = Form( REPLAY_DIR_PREFIX,"onlineGUI64/config_coinc_2018/LHRS_online.cfg");

		gSystem->Exec(Form("%sonline -P -f %s -r %d"                                   ,GUI_DIR      ,config_online_L ,runnumber              ));
                gSystem->Exec(Form("mv %stemp_%d.pdf %sleft_online_%d.pdf"                     ,SUM_DIR      ,runnumber       ,path_to_plots,runnumber));
                gSystem->Exec(Form("unlink %sleft_online_latest.pdf"                           ,SUM_DIR                                               ));
                gSystem->Exec(Form("ln -s %sleft_online_%d.pdf %sleft_online_latest.pdf"       ,path_to_plots,runnumber       ,SUM_DIR                ));

		//==================================
		// RHRS
		//==================================
		const char * config_online_R = Form( REPLAY_DIR_PREFIX,"onlineGUI64/config_coinc_2018/RHRS_online.cfg");

		gSystem->Exec(Form("%sonline -P -f %s -r %d"                                     ,GUI_DIR      ,config_online_R ,runnumber              ));
                gSystem->Exec(Form("mv %stemp_%d.pdf %sright_online_%d.pdf"                      ,SUM_DIR      ,runnumber       ,path_to_plots,runnumber));
                gSystem->Exec(Form("unlink %sright_online_latest.pdf"                            ,SUM_DIR                                               ));
                gSystem->Exec(Form("ln -s %sright_online_%d.pdf %sright_online_latest.pdf"       ,path_to_plots,runnumber       ,SUM_DIR                ));

		//==================================
		// Both arms combined
		//==================================
		const char * CONFIGCOINC     = Form( REPLAY_DIR_PREFIX,"onlineGUI64/config_coinc_2018/Coinc.cfg"        );
		const char * CONFIGCOINCPHYS = Form( REPLAY_DIR_PREFIX,"onlineGUI64/config_coinc_2018/Coinc_physics.cfg");

		gSystem->Exec(Form("%sonline -P -f %s -r %d"                                     ,GUI_DIR      ,CONFIGCOINC,runnumber                   ));
		gSystem->Exec(Form("mv %stemp_%d.pdf %scoinc_%d.pdf"                             ,SUM_DIR      ,runnumber       ,path_to_plots,runnumber));
		gSystem->Exec(Form("unlink %scoinc_latest.pdf"                                   ,SUM_DIR                                               ));
		gSystem->Exec(Form("ln -s %scoinc_%d.pdf %scoinc_latest.pdf"                     ,path_to_plots,runnumber       ,SUM_DIR                ));    
		gSystem->Exec(Form("ln -sf %scoinc_%d.pdf %scoinc_latest.pdf"                    ,path_to_plots,runnumber       ,path_to_plots          ));

		gSystem->Exec(Form("%sonline -P -f %s -r %d"                                     ,GUI_DIR      ,CONFIGCOINCPHYS ,runnumber            ));
                gSystem->Exec(Form("mv %stemp_%d.pdf %scoinc_physics_%d.pdf"                     ,SUM_DIR      ,runnumber       ,path_to_plots,runnumber));
                gSystem->Exec(Form("unlink %scoinc_physics_latest.pdf"                           ,SUM_DIR                                               ));
                gSystem->Exec(Form("ln -s %scoinc_physics_%d.pdf %scoinc_physics_latest.pdf"     ,path_to_plots,runnumber       ,SUM_DIR                ));    
                gSystem->Exec(Form("ln -sf %scoinc_physics_%d.pdf %scoinc_physics_latest.pdf"    ,path_to_plots,runnumber       ,path_to_plots          ));
	}

	exit(0);
}


