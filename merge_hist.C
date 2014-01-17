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

int nmass = 9;
//int mass[nmass]= {300,500,600,700,800,900,1500,2500,3000};

const char* channel[nmass]={
"histosnonres/Control_shower_1001.root",
"histosnonres/Control_shower_1000.root",
"histosnonres/Control_shower_101010.root",
//
"histosnonres/Control_shower_51010.root",
"histosnonres/Control_shower_151010.root",
//
"histosnonres/Control_shower_100010.root",
"histosnonres/Control_shower_102010.root",
//
"histosnonres/Control_shower_101000.root",
"histosnonres/Control_shower_101020.root"
};
int maxtodo=3; 
int todo[maxtodo]={0,1,2};
double masses[maxtodo] = {0,1,2};//,650,700,750};

const char* lege[nmass]={
"SM",
"SM VBFcuts",
"CV=1 C2V=1 C3=1",
//
"CV=0.5 C2V=1 C3=1",
"CV=1.5 C2V=1 C3=p1",
//
"CV=1 C2V=2 C3=1",
"CV=1 C2V=0 C3=1",
//
"CV=1 C2V=1 C3=0",
"CV=1 C2V=1 C3=2"
};

 TLegend *leg = new TLegend(0.55,0.550,0.99,0.99);
   leg->SetTextSize(0.03146853);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

 TLegend *leg2 = new TLegend(0.85,0.550,0.99,0.99);
   leg2->SetTextSize(0.03146853);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);

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
	8,9,11,12};//,
	//50,227};//,5,6,7,
	//8,9,11,12,
	//2,3,5,6,5,8};

for(int k=0;k<nplots;k++) for(int l=0;l<nmass;l++){
plots[l][k]->SetLineColor(sigcolor[l]);
plots[l][k]->SetLineStyle(0);
plots[l][k]->SetLineWidth(3);
cout<<"here "<<k<<" "<<l<<endl;
}
cout<<"here"<<endl;
TCanvas* PT_HAT = new TCanvas();
PT_HAT->cd(); 
int max=nmass;
//Double_t nevents[4]={20000.,20000.,20000.,20000.};//20000.,20000.,20000.,20000.};
double high[nplots]={1.2,1.2,1.2,2.1,1.2,
		     1.2,1.2,1.7,1.2,1.2,
		     1.2,1.2,1.2,1.2,1.2,
		     1.2,1.2,1.7,1.7,1.7,
                     1.2,1.2,1.2,1.2,1.2}; 

  for(int i=0;i<nplots;i++) {
  //if(i==16 || i==4 || i==5 || i==6 || i==7  || i==12) PT_HAT->SetLogy(1); else PT_HAT->SetLogy(0);
	
	for(int j=0;j<maxtodo;j++) leg->AddEntry(plots[todo[j]][i],lege[todo[j]],"l");
	//plots[0][i].Scale(plots[0][i].Integral());
        //plots[0][i].Scale(nevents[j]);
	plots[todo[0]][i].SetMaximum(high[i]*plots[todo[0]][i].GetMaximum());
	plots[todo[0]][i].Draw("Hist");
	leg->Draw("same");
	for(int j=1;j<maxtodo;j++) {
        //	plots[j][i].Scale(nevents[j]);
		//plots[j][i].Scale(plots[j][i].Integral());
		plots[todo[j]][i].Draw("Hist,same");
	
	}
	
	PT_HAT->SaveAs(namplots[i]);
	PT_HAT->Clear();
	leg->Clear();
  }
//////////////////////////////////////////////////////////////////////
// njets plot
vector<double> nj3,nj4,nj5,nj6,nj7,nj8,nj9,nj10,njmore,ntot;
    for(int j=0;j<maxtodo;j++) {
       cout<<lege[todo[j]]<<endl;
       // get number of njets
       int nbins = plots[todo[j]][0].GetNbinsX();
       //cout<<nbins<<endl;
       for(int k=0;k<nbins;k++) {
	double njets = plots[todo[j]][0].GetBinContent(k); 
        if (k == 3)      {nj3.push_back(njets); cout<<"3 jets "<<nj3[j]<<endl;}
        else if (k == 4) {nj4.push_back(njets); cout<<"4 jets "<<nj4[j]<<endl;}
        else if (k == 5) {nj5.push_back(njets); cout<<"5 jets "<<nj5[j]<<endl;}
        else if (k == 6) {nj6.push_back(njets); cout<<"6 jets "<<nj6[j]<<endl;}
        else if (k == 7) {nj7.push_back(njets); cout<<"7 jets "<<nj7[j]<<endl;}
        else if (k == 8) {nj8.push_back(njets); cout<<"8 jets "<<nj8[j]<<endl;}
        else if (k == 9) {nj9.push_back(njets); cout<<"9 jets "<<nj9[j]<<endl;}
        else if (k == 10) {nj10.push_back(njets); cout<<"10 jets "<<nj10[j]<<endl;}
        else {njmore.push_back(njets); cout<<"more "<< k<<" "<<njmore[j]<<endl;}
       } // close for bins
       ntot.push_back(nj3[j]+nj4[j]+nj5[j]+nj6[j]+nj7[j]+nj8[j]+nj9[j]+nj10[j]+njmore[j]);
       cout<<"total "<<ntot[j]<<endl;
   } // close for masses
   TMultiGraph *mg = new TMultiGraph();
   int nevents=1;//20000;
   //
   TGraphErrors *TreeJets = new TGraphErrors(1);
   TreeJets->SetMarkerStyle(kFullDotLarge); 
   TreeJets->SetLineColor(kMagenta);
   TreeJets->SetMarkerColor(kMagenta);
   TreeJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) TreeJets->SetPoint(j, masses[j], nj3[j]/nevents);
   //mg->Add(TreeJets,"L,P");//->Draw("L,P");
   //leg->AddEntry(TreeJets, " 3 jets", "LP");
   //
   TGraphErrors *fourJets = new TGraphErrors(1);
   fourJets->SetMarkerStyle(kFullDotLarge); 
   fourJets->SetLineColor(8);
   fourJets->SetMarkerColor(8);
   fourJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) fourJets->SetPoint(j, masses[j], nj4[j]/nevents);
   mg->Add(fourJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(fourJets, " 4 jets", "LP");
   //
   TGraphErrors *fiveJets = new TGraphErrors(1);
   fiveJets->SetMarkerStyle(kFullDotLarge); 
   fiveJets->SetLineColor(kCyan);
   fiveJets->SetMarkerColor(kCyan);
   fiveJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) fiveJets->SetPoint(j, masses[j], nj5[j]/nevents);
   mg->Add(fiveJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(fiveJets, " 5 jets", "LP");
   //
   TGraphErrors *sixJets = new TGraphErrors(1);
   sixJets->SetMarkerStyle(kFullDotLarge); 
   sixJets->SetLineColor(kBlue);
   sixJets->SetMarkerColor(kBlue);
   sixJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) sixJets->SetPoint(j, masses[j], nj6[j]/nevents);
   mg->Add(sixJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(sixJets, " 6 jets", "LP");
   //
   TGraphErrors *sevenJets = new TGraphErrors(1);
   sevenJets->SetMarkerStyle(kFullDotLarge); 
   sevenJets->SetLineColor(kMagenta);
   sevenJets->SetMarkerColor(kMagenta);
   sevenJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) sevenJets->SetPoint(j, masses[j], nj7[j]/nevents);
   mg->Add(sevenJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(sevenJets, " 7 jets", "LP");
   //
   TGraphErrors *eightJets = new TGraphErrors(1);
   eightJets->SetMarkerStyle(kFullDotLarge); 
   eightJets->SetLineColor(kGreen);
   eightJets->SetMarkerColor(kGreen);
   eightJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) eightJets->SetPoint(j, masses[j], nj8[j]/nevents);
   mg->Add(eightJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(eightJets, " 8 jets", "LP");
   //
   TGraphErrors *nineJets = new TGraphErrors(1);
   nineJets->SetMarkerStyle(kFullDotLarge); 
   nineJets->SetLineColor(50);
   nineJets->SetMarkerColor(50);
   nineJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) nineJets->SetPoint(j, masses[j], nj9[j]/nevents);
   mg->Add(nineJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(nineJets, " 9 jets", "LP");
   //
   TGraphErrors *tenJets = new TGraphErrors(1);
   tenJets->SetMarkerStyle(kFullDotLarge); 
   tenJets->SetLineColor(kGray);
   tenJets->SetMarkerColor(kGray);
   tenJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) tenJets->SetPoint(j, masses[j], nj10[j]/nevents);
   mg->Add(tenJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(tenJets, " 10 jets", "LP");
   //
   TGraphErrors *moreJets = new TGraphErrors(1);
   moreJets->SetMarkerStyle(kFullDotLarge); 
   moreJets->SetLineColor(kPink);
   moreJets->SetMarkerColor(kPink);
   moreJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) moreJets->SetPoint(j, masses[j], njmore[j]/nevents);
   mg->Add(moreJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(moreJets, " > 6 jets", "LP");
   //
   TGraphErrors *totJets = new TGraphErrors(1);
   totJets->SetMarkerStyle(kFullDotLarge); 
   totJets->SetLineColor(kBlack);
   totJets->SetMarkerColor(kBlack);
   totJets->SetLineWidth(3);
   for(int j=0;j<maxtodo;j++) totJets->SetPoint(j, masses[j], ntot[j]/nevents);
   mg->Add(totJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(totJets, " total", "LP");
   //
  mg->Draw("AP");
  leg2->Draw();
  PT_HAT->SetLogy(0);
  PT_HAT->SaveAs("njets_flat.png");
//PT_HAT->Close(); 

}
