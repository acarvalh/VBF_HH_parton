{
/////////////////////////////////////////////
// here we put the options for ploting
TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
  defaultStyle->SetOptStat(0000);
  defaultStyle->SetOptFit(000); 
  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.05);
  defaultStyle->SetPadBottomMargin(0.13);
  defaultStyle->SetPadLeftMargin(0.13);
  defaultStyle->SetPadRightMargin(0.02);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetCanvasDefH(600);
  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0); 
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.25,"Y");
  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  defaultStyle->SetNdivisions(505,"Y");
  defaultStyle->SetLegendBorderSize(0);  // For the axis titles:

    defaultStyle->SetTitleColor(1, "XYZ");
    defaultStyle->SetTitleFont(42, "XYZ");
    defaultStyle->SetTitleSize(0.06, "XYZ");
 
    // defaultStyle->SetTitleYSize(Float_t size = 0.02);
    defaultStyle->SetTitleXOffset(0.9);
    defaultStyle->SetTitleYOffset(1.05);
    // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    // For the axis labels:
    defaultStyle->SetLabelColor(1, "XYZ");
    defaultStyle->SetLabelFont(42, "XYZ");
    defaultStyle->SetLabelOffset(0.007, "XYZ");
    defaultStyle->SetLabelSize(0.04, "XYZ");

    // For the axis:
    defaultStyle->SetAxisColor(1, "XYZ");
    defaultStyle->SetStripDecimals(kTRUE);
    defaultStyle->SetTickLength(0.03, "XYZ");
    defaultStyle->SetNdivisions(510, "XYZ");
    defaultStyle->SetPadTickX(1);   // To get tick marks on the opposite side of the frame
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();
/////////////////////////////////////////////////////

int nmass = 24;
//int mass[nmass]= {300,500,600,700,800,900,1500,2500,3000};

const char* channel[nmass]={
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
//
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
//
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root",
"histos/Control_shower_260.root"
};

//const char* lege[nmass]={"300 GeV","500 GeV","600 GeV","700 GeV","800 GeV","900 GeV","1500 GeV","2500 GeV","3000 GeV"};

const char* lege[nmass]={"250 GeV","300 GeV","350 GeV","400 GeV","450 GeV","500 GeV","550 GeV","600 GeV","650 GeV","700 GeV","750 GeV","800 GeV","850 GeV","900 GeV","950 GeV","1000 GeV","1100 GeV","1200 GeV","1300 GeV","1400 GeV","1500 GeV","2000 GeV","2500 GeV","3000 GeV"};
//for (int k=0; k<nmass; k++) channel[k] = Form("Control_shower_%d.root", k);


int nplots =25;

 TLegend *leg = new TLegend(0.65,0.250,0.95,0.90);
   leg->SetTextSize(0.03146853);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

 leg->SetHeader("BKG to 8TeV");


   const char* namplots[nplots] = {
        "DiPhotonMass.png","DiPhotonE.png", "DiPhotonPt.png", "DiPhotonEta.png", "DiPhotonPhi.png", 
	"njets.png", //1
	"DijetMass.png","DijetE.png","DijetPt.png","DijetEta.png","DijetPhi.png", // 5
	"radMass.png","radE.png","radPt.png","radEta.png","radPhi.png", // 5
	"DRphph.png",  //DPhiphph, // 2
	"DRhh.png", //DPhihh, // 2
	"DRjj.png",// DPhijj, // 2
	"empty.png","empty2.png","lead.png","sublead.png","fattagreco.png","fattagrecon.png"
	//0,fabs(Higgs_ph.rapidity()-Dijet.rapidity())//ph1_pt/ph2_pt,jpt[nlead]/jpt[nsublead] // 2
    }; //25


TH1D* plots[nmass][nplots];//[file][plot]
TFile *file[nmass];
for(int i=0;i<nmass;i++){
 TFile *file[i] = TFile::Open(channel[i]);
 TH1D* plots[i][0] = (TH1D* ) file[i]->Get("PhotonsMass_ct4;1");//
 TH1D* plots[i][1] = (TH1D* ) file[i]->Get("dipho_E_ct4;1");
 TH1D* plots[i][2] = (TH1D* ) file[i]->Get("dipho_pt_ct4;1");
 TH1D* plots[i][3] = (TH1D* ) file[i]->Get("dipho_eta_ct4;1");
 TH1D* plots[i][4] = (TH1D* ) file[i]->Get("dipho_phi_ct4;1");
 TH1D* plots[i][5] = (TH1D* ) file[i]->Get("njets_passing_kLooseID_ct4;1");
 TH1D* plots[i][6] = (TH1D* ) file[i]->Get("JetsMass_ct4;1");
 TH1D* plots[i][7] = (TH1D* ) file[i]->Get("dijet_E_ct4;1");
 TH1D* plots[i][8] = (TH1D* ) file[i]->Get("dijet_pt_ct4;1");
 TH1D* plots[i][9] = (TH1D* ) file[i]->Get("dijet_eta_ct4;1");
 TH1D* plots[i][10] = (TH1D* ) file[i]->Get("dijet_phi_ct4;1");
 TH1D* plots[i][11] = (TH1D* ) file[i]->Get("RadMass_ct4;1");
 TH1D* plots[i][12] = (TH1D* ) file[i]->Get("radion_E_ct4;1");
 TH1D* plots[i][13] = (TH1D* ) file[i]->Get("radion_pt_ct4;1");
 TH1D* plots[i][14] = (TH1D* ) file[i]->Get("radion_eta_ct4;1");
 TH1D* plots[i][15] = (TH1D* ) file[i]->Get("radion_phi_ct4;1");
 TH1D* plots[i][16] = (TH1D* ) file[i]->Get("DRphph_ct4;1");
 TH1D* plots[i][17] = (TH1D* ) file[i]->Get("DRhh_ct4;1");
 TH1D* plots[i][18] = (TH1D* ) file[i]->Get("DRjj_ct4;1");
 TH1D* plots[i][19] = (TH1D* ) file[i]->Get("aaass_ct4;1");
 TH1D* plots[i][20] = (TH1D* ) file[i]->Get("jjass_ct4;1");
 TH1D* plots[i][21] = (TH1D* ) file[i]->Get("leadsel_ct4;1");
 TH1D* plots[i][22] = (TH1D* ) file[i]->Get("sublead_ct4;1");
 TH1D* plots[i][23] = (TH1D* ) file[i]->Get("fattagreco_ct4;1");
 TH1D* plots[i][24] = (TH1D* ) file[i]->Get("fattagrecon_ct4;1");
} 

const int sigcolor[nmass]={2,3,5,6,7,8,9,11,12,2,3,5,6,7,8,9,11,12,2,3,5,6,7,8};
for(int k=0;k<nplots;k++) for(int l=0;l<nmass;l++){
//plots[0][i]->SetFillColor(kCyan+1);
//plots[1][i]->SetFillColor(kCyan+2);// terra qcd_30_8TeV_ff.root
//int k = 0;
plots[l][k]->SetLineColor(sigcolor[l]);//black qcd_30_8TeV_pf
plots[l][k]->SetLineStyle(0);
plots[l][k]->SetLineWidth(3);

//
//cout<<"here"<<endl;
}

TCanvas* PT_HAT = new TCanvas();
PT_HAT->cd(); 

for(int i=0;i<nplots;i++) {
//if(i==16 || i==4 || i==5 || i==6 || i==7  || i==12) PT_HAT->SetLogy(1); else PT_HAT->SetLogy(0);

	for(int j=0;j<nmass;j++) leg->AddEntry(plots[j][i],lege[j],"l");
	//plots[0][i].Scale(plots[0][i].Integral());
	plots[0][i].SetMaximum(4*plots[0][i].GetMaximum());
	plots[0][i].Draw("Hist");
	leg->Draw("same");
	for(int j=1;j<nmass;j++) {
	//	plots[j][i].Scale(plots[j][i].Integral());
		plots[j][i].Draw("Hist,same");
	}
	
	PT_HAT->SaveAs(namplots[i]);
	PT_HAT->Clear();
	leg->Clear();
}

PT_HAT->Close(); 

}
