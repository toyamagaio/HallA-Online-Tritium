//E/p distribution
//Tritium (e,e'p) SRC experiment
//Reynier Cruz Torres, Apr 2018
void E_over_p(Int_t flag, TString drawoption){
	gStyle->SetOptStat(0);
	TTree *tree = (TTree*)gDirectory->Get("T");

	// ---------------------------------------------------------------------------------
	// LHRS
	if(flag ==1){
		TH1F *h = new TH1F("h","E/p;E/p",100,0,1.5);
		h->GetYaxis()->CenterTitle();	h->GetXaxis()->CenterTitle();
		tree->Draw("(L.prl1.e+L.prl2.e)/L.tr.p[0]/1000.>>h","L.tr.n==1",drawoption);
	}
	// ---------------------------------------------------------------------------------
	// RHRS
	else if(flag ==2){
		TH1F *h1 = new TH1F("h1","E/p;E/p",100,0,2);
		h1->GetYaxis()->CenterTitle();	h1->GetXaxis()->CenterTitle();
		tree->Draw("(R.ps.e+R.sh.e)/R.tr.p[0]/1000.>>h1","R.tr.n==1",drawoption);
	}
	// ---------------------------------------------------------------------------------
}
