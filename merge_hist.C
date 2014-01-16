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

int nmass = 11;
//int mass[nmass]= {300,500,600,700,800,900,1500,2500,3000};

const char* channel[nmass]={
"histos/Control_shower_260.root",
"histos/Control_shower_0.root",
"histos/Control_shower_1.root",
"histos/Control_shower_137.root",
"histos/Control_shower_450.root",
"histos/Control_shower_500.root",
"histos/Control_shower_550.root",
"histos/Control_shower_600.root",
"histos/Control_shower_650.root",
"histos/Control_shower_700.root",
"histos/Control_shower_750.root",
};

//const char* lege[nmass]={"300 GeV","500 GeV","600 GeV","700 GeV","800 GeV","900 GeV","1500 GeV","2500 GeV","3000 GeV"};

//const char* lege[nmass]={"260 GeV","300 GeV","350 GeV","400 GeV","450 GeV","500 GeV","550 GeV","600 GeV","650 GeV","700 GeV","750 GeV"};
const char* lege[nmass]={"260 GeV","cg =0","cg=1","cg=0.0137","450 GeV","500 GeV","550 GeV","600 GeV","650 GeV","700 GeV","750 GeV"};
//for (int k=0; k<nmass; k++) channel[k] = Form("Control_shower_%d.root", k);

 TLegend *leg = new TLegend(0.75,0.550,0.99,0.99);
   leg->SetTextSize(0.03146853);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

 leg->SetHeader("BKG to 8TeV");

 int nplots =25;
 const char* namplots[nplots] = {
           "njets_passing_kLooseID_ct4.png",
           "H1hist.png",
           "H1histpt.png",
           "H1histeta.png",
           "H1histphi.png",
           "H2hist.png",
           "H2histpt.png",
           "H2histeta.png",
           "H2histphi.png",
           "RadMass_ct4.png",
           "radion_pt_ct4.png",
           "radion_eta_ct4.png",
           "radion_phi_ct4.png",  
           "fattag_ct4.png",
           "btag_ct4.png",
           "bmistag_ct4.png",
           "j1histpt.png",
           "j1histeta.png",
           "j2histpt.png",
           "j2histeta.png",
           "RadMass_ct4vbf.png",
           "Deta.png",
           "fattag_ct4vbf.png",
           "btag_ct4vbf.png",
           "bmistag_ct4vbf.png"
    }; //25

TH1D* plots[nmass][nplots];//[file][plot]
TFile *file[nmass];
for(int i=0;i<nmass;i++){
 TFile *file[i] = TFile::Open(channel[i]);
 TH1D* plots[i][0] = (TH1D* ) file[i]->Get("njets_passing_kLooseID_ct4;1"); 
 TH1D* plots[i][1] = (TH1D* ) file[i]->Get("H1hist;1"); 
 TH1D* plots[i][2] = (TH1D* ) file[i]->Get("H1histpt;1"); 
 TH1D* plots[i][3] = (TH1D* ) file[i]->Get("H1histeta;1"); 
 TH1D* plots[i][4] = (TH1D* ) file[i]->Get("H1histphi;1");
 TH1D* plots[i][5] = (TH1D* ) file[i]->Get("H2hist;1"); 
 TH1D* plots[i][6] = (TH1D* ) file[i]->Get("H2histpt;1"); 
 TH1D* plots[i][7] = (TH1D* ) file[i]->Get("H2histeta;1"); 
 TH1D* plots[i][8] = (TH1D* ) file[i]->Get("H2histphi;1"); 
 TH1D* plots[i][9] = (TH1D* ) file[i]->Get("RadMass_ct4;1"); 
 TH1D* plots[i][10] = (TH1D* ) file[i]->Get("radion_pt_ct4;1"); 
 TH1D* plots[i][11] = (TH1D* ) file[i]->Get("radion_eta_ct4;1"); 
 TH1D* plots[i][12] = (TH1D* ) file[i]->Get("radion_phi_ct4;1");
 TH1D* plots[i][13] = (TH1D* ) file[i]->Get("fattag_ct4;1"); 
 TH1D* plots[i][14] = (TH1D* ) file[i]->Get("btag_ct4;1"); 
 TH1D* plots[i][15] = (TH1D* ) file[i]->Get("bmistag_ct4;1"); 
 TH1D* plots[i][16] = (TH1D* ) file[i]->Get("j1histpt;1"); 
 TH1D* plots[i][17] = (TH1D* ) file[i]->Get("j1histeta;1"); 
 TH1D* plots[i][18] = (TH1D* ) file[i]->Get("j1histeta;1");//"j2histpt;1"); 
 TH1D* plots[i][19] = (TH1D* ) file[i]->Get("j1histeta;1");//"j2histeta;1"); 
 TH1D* plots[i][20] = (TH1D* ) file[i]->Get("RadMass_ct4vbf;1"); 
 TH1D* plots[i][21] = (TH1D* ) file[i]->Get("Deta;1"); 
 TH1D* plots[i][22] = (TH1D* ) file[i]->Get("fattag_ct4vbf;1"); 
 TH1D* plots[i][23] = (TH1D* ) file[i]->Get("btag_ct4vbf;1"); 
 TH1D* plots[i][24] = (TH1D* ) file[i]->Get("bmistag_ct4vbf;1");
} 

const int sigcolor[nmass]={
	2,3,5,6,7,
	8,9,11,12,
	50,227};//,5,6,7,
	//8,9,11,12,
	//2,3,5,6,5,8};
for(int k=0;k<nplots;k++) for(int l=0;l<nmass;l++){
//plots[0][i]->SetFillColor(kCyan+1);
//plots[1][i]->SetFillColor(kCyan+2);// terra qcd_30_8TeV_ff.root
//int k = 0;
plots[l][k]->SetLineColor(sigcolor[l]);//black qcd_30_8TeV_pf
plots[l][k]->SetLineStyle(0);
plots[l][k]->SetLineWidth(3);

//
cout<<"here "<<k<<" "<<l<<endl;
}
cout<<"here"<<endl;
TCanvas* PT_HAT = new TCanvas();
PT_HAT->cd(); 
int max=4;//nmass;
double nevents[4]={20000,50000,50000,50000};
for(int i=0;i<nplots;i++) {
//if(i==16 || i==4 || i==5 || i==6 || i==7  || i==12) PT_HAT->SetLogy(1); else PT_HAT->SetLogy(0);
	
	for(int j=0;j<max;j++) leg->AddEntry(plots[j][i],lege[j],"l");
	//plots[0][i].Scale(plots[0][i].Integral());
        //plots[0][i].Scale(nevents[j]);
	plots[0][i].SetMaximum(1.2*plots[0][i].GetMaximum());
	plots[0][i].Draw("Hist");
	leg->Draw("same");
	for(int j=1;j<max;j++) {
        	//plots[j][i].Scale(nevents[j]);
		//plots[j][i].Scale(plots[j][i].Integral());
		plots[j][i].Draw("Hist,same");
	
	}
	
	PT_HAT->SaveAs(namplots[i]);
	PT_HAT->Clear();
	leg->Clear();
}

PT_HAT->Close(); 

}
