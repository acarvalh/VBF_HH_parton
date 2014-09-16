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
  defaultStyle->SetLabelSize(0.07,"XY");
  defaultStyle->SetLabelFont(46,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.25,"Y");
  defaultStyle->SetTitleSize(0.07,"Y");
  defaultStyle->SetTitleFont(44, "XYZ");
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
int nmass = 68;
//int mass[nmass]= {300,500,600,700,800,900,1500,2500,3000};

const char* channel[nmass]={
"nonresonant/Control_shower_1001.root", //1
"spin0/shower/Control_shower_260.root",//2
"spin0/shower/Control_shower_300.root",//3
"spin0/shower/Control_shower_350.root",//4
//"olivierMar14/MGraviton__500_HHtobbbb_NoPileUp.root",//2
//"olivierMar14/MGraviton__500_HHtobbbb_50PileUp.root",//3
//"olivierMar14/MGraviton__500_HHtobbbb_140PileUp.root",//4
"spin0/shower/Control_shower_400.root",//5
"spin0/shower/Control_shower_450.root",//6
"spin0/shower/Control_shower_500.root",//7
"spin0/shower/Control_shower_550.root",//8
"spin0/shower/Control_shower_600.root",//9
"spin0/shower/Control_shower_650.root",//10
"spin0/shower/Control_shower_700.root",//11
"spin0/shower/Control_shower_750.root",//12
"spin0/shower/Control_shower_800.root",//13
"spin0/shower/Control_shower_850.root",//14
"spin0/shower/Control_shower_900.root",//15
"spin0/shower/Control_shower_950.root",//16
"spin0/shower/Control_shower_1000.root",//17
"spin0/shower/Control_shower_1050.root",//18
"spin0/shower/Control_shower_1100.root",//19
"spin0/shower/Control_shower_1150.root",//20
"spin0/shower/Control_shower_1200.root",//21
"spin0/shower/Control_shower_1250.root",//22
"spin0/shower/Control_shower_1300.root",//23
"spin0/shower/Control_shower_1350.root",//24
"spin0/shower/Control_shower_1400.root",//25
"spin0/shower/Control_shower_1450.root",//26
"spin0/shower/Control_shower_1500.root",//27
// bkgs 4b
"4bsbkg/Control_shower_0.root",//28
"4bsbkg/Control_shower_1.root",//29
"4bsbkg/Control_shower_2.root",//30
"4bsbkg/Control_shower_3.root",//31
"4bsbkg/Control_shower_4.root",//32
"4bsbkg/Control_shower_5.root",//33
"4bsbkg/Control_shower_7.root",//34
"4bsbkg/Control_shower_8.root",//35
//
// WWbb
"nonresonant/Control_shower_1000.root",
"bulk_graviton_mad_WWbb/Control_shower_260.root",
"bulk_graviton_mad_WWbb/Control_shower_300.root",
"bulk_graviton_mad_WWbb/Control_shower_350.root",
"bulk_graviton_mad_WWbb/Control_shower_400.root",
"bulk_graviton_mad_WWbb/Control_shower_450.root",
"bulk_graviton_mad_WWbb/Control_shower_500.root",
"bulk_graviton_mad_WWbb/Control_shower_550.root",
"bulk_graviton_mad_WWbb/Control_shower_550.root",
"bulk_graviton_mad_WWbb/Control_shower_550.root",
"bulk_graviton_mad_WWbb/Control_shower_700.root",
"bulk_graviton_mad_WWbb/Control_shower_750.root",
"bulk_graviton_mad_WWbb/Control_shower_800.root",
"bulk_graviton_mad_WWbb/Control_shower_850.root",
"bulk_graviton_mad_WWbb/Control_shower_900.root",
"bulk_graviton_mad_WWbb/Control_shower_950.root",
"bulk_graviton_mad_WWbb/Control_shower_1000.root",
"bulk_graviton_mad_WWbb/Control_shower_1050.root",
"bulk_graviton_mad_WWbb/Control_shower_1100.root",
"bulk_graviton_mad_WWbb/Control_shower_1150.root",
"bulk_graviton_mad_WWbb/Control_shower_1200.root",
"bulk_graviton_mad_WWbb/Control_shower_1250.root",
"bulk_graviton_mad_WWbb/Control_shower_1300.root",
"bulk_graviton_mad_WWbb/Control_shower_1350.root",
"bulk_graviton_mad_WWbb/Control_shower_1400.root",
"bulk_graviton_mad_WWbb/Control_shower_1450.root",
"bulk_graviton_mad_WWbb/Control_shower_1500.root",
// WWbb BKG
"4bsbkg/Control_shower_9.root",
"4bsbkg/Control_shower_10.root",
//
"4bsbkg/Control_shower_11.root", // veronica
"4bsbkg/Control_shower_12.root",
"4bsbkg/Control_shower_13.root",
"4bsbkg/Control_shower_14.root"
//"histos/Madgraph0/Control_shower_260.root",
//"histos/Madgraph0_0137/Control_shower_260.root",
//"histos/Madgraph1/Control_shower_260.root",
//
//"histos/Madgraph0/Control_shower_500.root",
//"histos/Madgraph0_0137/Control_shower_500.root",
//"histos/Madgraph1/Control_shower_500.root"
};

//const char* lege[nmass]={"300 GeV","500 GeV","600 GeV","700 GeV","800 GeV","900 GeV","1500 GeV","2500 GeV","3000 GeV"};

const char* lege[nmass]={"SM H",
                         "260 GeV","300 GeV","350 GeV",
                         //"no pileup", "LHC Run 3","HL-LHC ",
                         "400 GeV","450 GeV","500 GeV",
                         "550 GeV","600 GeV","650 GeV","700 GeV","750 GeV",
			 "800 GeV","850 GeV","900 GeV","950 GeV",
                         "1000 GeV","1050 GeV","1100 GeV","1150 GeV",
                         "1200 GeV","1250 GeV","1300 GeV","1350 GeV",
                         "1400 GeV","1450 GeV", "1500 GeV",
			"z(bb) bb jj", // 10k
			"z(bb) z(bb) jj", // 10k
			"dd> 4bdd", // 10k
			"gg> 4bddbar", // 10k
			"gg> 4buubar", // 10k
			"ud> 4bud", // 10k
			"uu> 4buu", // 2.3k events
			"4b2j alpgen",
                        "SM HH","260 GeV","300 GeV","350 GeV","400 GeV","450 GeV","calchep 500 GeV",
                         "550 GeV","600 GeV","650 GeV","700 GeV","750 GeV",
			 "800 GeV","850 GeV","900 GeV","950 GeV",
                         "1000 GeV","1050 GeV","1100 GeV","1150 GeV",
                         "1200 GeV","1250 GeV","1300 GeV","1350 GeV",
                         "1400 GeV","1450 GeV", "1500 GeV",
			"t(wb)t(wb) jj", // 100kevents
			"wwbb jj", // 10k
			 " 400 GeV 14 TeV - veronica", "M_400_8tev",
"M_400_10tev",
"M_400_14tev"
			// "Mad cg=0 (500)","Mad cg=0.0137 (500)","Mad cg=1 (500)"
                        };

/*
int maxtodo=4;//6; 
int todo[maxtodo]={ 1,2,3,4};//,6,64,65,66,67};
double masses[maxtodo] = { 0,0,0,0};//,1, 260, 260, 260, 260 };//,
*/

// 4b
int maxtodo=6; 
int todo[maxtodo]={1,2,3,4,5,6};//3 ,2 ,1 };//27,28,34, // 0,21,16,13, 9 ,5 ,
double masses[maxtodo] = {1,2,3,4,5,6};//,

/////////////////
// 4b all
/*
int maxtodo=25;//25; 
int todo[maxtodo]={ 1, 2, 3, 4, 5,
                    7, 8, 9, 10, 11,
		    12,13,14,15,16,17,18,19,20,21,22,23,24,25};
double masses[maxtodo] = { 260, 300 ,350 ,400 ,450 ,
                           500 ,550 ,600, 650 ,700 ,
                           750,800 ,850  ,950 ,1000,
                          1050,1100,1150,1200,1250,
                          1300,1350,1400,1450,1500};//,
*/
/*
// wwbb
int maxtodo=12; 
int todo[maxtodo]={   36,37,40,42,44,47,50,53,56,61,62,63};
double masses[maxtodo] = { 0, 260, 300 ,350 ,400 ,450 ,500 ,550 ,600,650 ,700 ,750};//,
*/
 TLegend *leg = new TLegend(0.65,0.60,0.99,0.99);
   leg->SetTextSize(0.04146853);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);

 TLegend *leg2 = new TLegend(0.85,0.60,0.5,0.9);
   leg2->SetTextSize(0.03146853);
   leg2->SetLineColor(1);
   leg2->SetLineStyle(1);
   leg2->SetLineWidth(1);
   leg2->SetFillColor(0);

 TLegend *leg3 = new TLegend(0.20,0.60,0.5,0.9);
   leg3->SetTextSize(0.03146853);
   leg3->SetLineColor(1);
   leg3->SetLineStyle(1);
   leg3->SetLineWidth(1);
   leg3->SetFillColor(0);

 //leg->SetHeader("MX = 500 GeV");

 int nplots =42;
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
           "RadPt_ct4vbf.png",
           "Deta.png",
           "fattag_ct4vbf.png",
           "btag_ct4vbf.png",
           "bmistag_ct4vbf.png",
           "cat.png",
           "detahh.png",
	   "gen_higgs.png",
	   "btagselected.png",
"E1histpt.png",  
"E1histeta.png",  
"E2histpt.png",  
"E2histeta.png",  
"LLMass_ct4.png",  
"LL_pt_ct4.png",  
"LL_Deta.png",  
"MetMass_ct4.png",
"Detabb.png",
"DRll_test.png",
"DRbb_test.png",
"MetMass1_ct4.png"
    }; //25

 const char* namplotspdf[nplots] = {
           "njets_passing_kLooseID_ct4.pdf",
           "H1hist.pdf",
           "H1histpt.pdf",
           "H1histeta.pdf",
           "H1histphi.pdf",
           "H2hist.pdf",
           "H2histpt.pdf",
           "H2histeta.pdf",
           "H2histphi.pdf",
           "RadMass_ct4.pdf",
           "radion_pt_ct4.pdf",
           "radion_eta_ct4.pdf",
           "radion_phi_ct4.pdf",  
           "fattag_ct4.pdf",
           "btag_ct4.pdf",
           "bmistag_ct4.pdf",
           "j1histpt.pdf",
           "j1histeta.pdf",
           "j2histpt.pdf",
           "j2histeta.pdf",
           "RadMass_ct4vbf.pdf",
           "RadPt_ct4vbf.pdf",
           "Deta.pdf",
           "fattag_ct4vbf.pdf",
           "btag_ct4vbf.pdf",
           "bmistag_ct4vbf.pdf",
           "cat.pdf",
           "detahh.pdf",
	   "gen_higgs.pdf",
	   "btagselected.pdf",
"E1histpt.pdf",  
"E1histeta.pdf",  
"E2histpt.pdf",  
"E2histeta.pdf",  
"LLMass_ct4.pdf",  
"LL_pt_ct4.pdf",  
"LL_Deta.pdf",  
"MetMass_ct4.pdf",
"Detabb.pdf",
"DRll_test.pdf",
"DRbb_test.pdf",
"MetMass1_ct4.pdf"
    }; //25

TH1D* plots[maxtodo][nplots];//[file][plot]
TFile *file[maxtodo];
for(int i=0;i<maxtodo;i++){
 TFile *file[i] = TFile::Open(channel[todo[i]]);
 cout<<channel[todo[i]]<<endl;
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
 TH1D* plots[i][18] = (TH1D* ) file[i]->Get("j2histpt;1");//"j2histpt;1"); 
 TH1D* plots[i][19] = (TH1D* ) file[i]->Get("j2histeta;1");//"j2histeta;1"); 
 TH1D* plots[i][20] = (TH1D* ) file[i]->Get("RadMass_ct4vbf;1"); 
 TH1D* plots[i][21] = (TH1D* ) file[i]->Get("RadPt_ct4vbf;1"); 
 TH1D* plots[i][22] = (TH1D* ) file[i]->Get("Deta;1"); 
 TH1D* plots[i][23] = (TH1D* ) file[i]->Get("fattag_ct4vbf;1"); 
 TH1D* plots[i][24] = (TH1D* ) file[i]->Get("btag_ct4vbf;1"); 
 TH1D* plots[i][25] = (TH1D* ) file[i]->Get("bmistag_ct4vbf;1");
 TH1D* plots[i][26] = (TH1D* ) file[i]->Get("cat_ct4;1");
 TH1D* plots[i][27] = (TH1D* ) file[i]->Get("detahh_ct4;1");
 TH1D* plots[i][28] = (TH1D* ) file[i]->Get("gen_higgs;1");
 TH1D* plots[i][29] = (TH1D* ) file[i]->Get("btagselected;1");
 TH1D* plots[i][30] = (TH1D* ) file[i]->Get("E1histpt;1");   
 TH1D* plots[i][31] = (TH1D* ) file[i]->Get("E1histeta;1");   
 TH1D* plots[i][32] = (TH1D* ) file[i]->Get("E2histpt;1");   
 TH1D* plots[i][33] = (TH1D* ) file[i]->Get("E2histeta;1");   
 TH1D* plots[i][34] = (TH1D* ) file[i]->Get("LLMass_ct4;1");   
 TH1D* plots[i][35] = (TH1D* ) file[i]->Get("LL_pt_ct4;1");   
 TH1D* plots[i][36] = (TH1D* ) file[i]->Get("LL_Deta;1");   
 TH1D* plots[i][37] = (TH1D* ) file[i]->Get("MetMass_ct4;1"); //
 TH1D* plots[i][38] = (TH1D* ) file[i]->Get("Detabb;1");//Detabb
 TH1D* plots[i][39] = (TH1D* ) file[i]->Get("DRll_test;1");//Detabb
 TH1D* plots[i][40] = (TH1D* ) file[i]->Get("DRbb;1");//Detabb
 TH1D* plots[i][41] = (TH1D* ) file[i]->Get("MetMass1_ct4;1"); //
} 

const int sigcolor[nmass]={
	2 ,1  ,8 ,6 ,7 ,
	9 ,30  ,11,12,70,
	50,227,5, 6, 7,
	227, 9, 11, 12,8,
	2, 3, 5,  6, 5,
        8,9,15,8, 70, 
        5,  6, 4,50,3,
        5 ,6 ,7 ,
	9 ,30  ,11,12,70,
	50,227,5, 6, 7,
	227, 9, 11, 12,8,
	2, 3, 5,  6, 5,
        8,9,15,8, 70, 
        5,40,2,3,5
};

for(int k=0;k<nplots;k++) for(int l=0;l<maxtodo;l++){
plots[l][k]->SetLineColor(sigcolor[todo[l]]);
//plots[l][k]->SetLineStyle(0);
//if(l<3) plots[l][k]->SetLineStyle(2); else 
plots[l][k]->SetLineStyle(0);
//if(l==1) plots[l][k]->SetLineStyle(2); 
plots[l][k]->SetLineWidth(3);
cout<<"here "<<k<<" "<<l<<endl;
}
cout<<"here"<<endl;
TCanvas* PT_HAT = new TCanvas();
PT_HAT->cd(); 
int max=nmass;
//Double_t nevents[4]={20000.,20000.,20000.,20000.};//20000.,20000.,20000.,20000.};
double high[nplots]={1,1,1.2,1.1,1.2,
		     1,1.2,1.7,1.7,1.7,
		     1.2,1.2,1.2,1.2,1.2,
		     1.2,1.2,1.7,2.7,2.7,
                     1.2,1.2,2.5,1.2,1.2,
                     1.2,1.9,1.7,2.5,2.7,
                     1.7,1.7,1.7,1.7,1.7,
		     1.7,1.5,1.2,1.7,1.7,1.7}; 
  for(int i=0;i<nplots;i++) {
  //if(i==16 || i==4 || i==5 || i==6 || i==7  || i==12) PT_HAT->SetLogy(1); else PT_HAT->SetLogy(0);
	for(int j=0;j<maxtodo;j++) {
        leg->AddEntry(plots[j][i],lege[todo[j]],"l");
        cout<<"here "<<j<<" "<<i<<endl;
        }
	plots[0][i].Scale(1./plots[0][i].Integral());
        //plots[0][i].Scale(nevents[j]);
	plots[0][i].SetMaximum(high[i]*plots[0][i].GetMaximum());
	plots[0][i].Draw("Hist");
	leg->Draw("same");
	for(int j=1;j<maxtodo;j++) {
        //	plots[j][i].Scale(nevents[j]);
		//plots[j][i].Scale(plots[j][i].Integral());
		if(i!=26)plots[j][i].Scale(1./plots[j][i].Integral());
	//plots[j][i].SetMaximum(high[i]*plots[j][i].GetMaximum());
		plots[j][i].Draw("Hist,same");
	
	}
        if(i==0) {TLine li(6,0.00001,6,0.25); li->Draw("same");}
        if(i==1 || i==5) {TLine li3(125,0.0002,125,0.09); li3->Draw("same");}
        if(i==29) {TLine li2(4,0.00001,4,0.16); li2->Draw("same");}
        if(i==9) {TLine li2(500,0.00001,500,0.040); li2->Draw("same");}
	if( i==1 || i==5) {
           PT_HAT->SetLogy(1); 
           //TLine li(500,0.00001,500,0.001);
           //li->Draw("same");
        } else PT_HAT->SetLogy(0); 	
	PT_HAT->SaveAs(namplots[i]);
	if(i==28 || i==22 || i==9 || i==1 || i==0 || i==29) PT_HAT->SaveAs(namplotspdf[i]);
	PT_HAT->Clear();
	leg->Clear();
  }
//////////////////////////////////////////////////////////////////////
// njets plot
vector<double> nj3,nj4,nj5,nj6,nj7,nj8,nj9,nj10,njmore,ntot;
    for(int j=0;j<maxtodo;j++) {
       cout<<lege[todo[j]]<<endl;
       // get number of njets
       int nbins = plots[j][0].GetNbinsX();
       //cout<<nbins<<endl;
       for(int k=0;k<nbins;k++) {
	double njets = plots[j][0].GetBinContent(k); 
        if (k == 4)      {nj3.push_back(njets); cout<<"3 jets "<<nj3[j]<<endl;}
        else if (k == 5) {nj4.push_back(njets); cout<<"4 jets "<<nj4[j]<<endl;}
        else if (k == 6) {nj5.push_back(njets); cout<<"5 jets "<<nj5[j]<<endl;}
        else if (k == 7) {nj6.push_back(njets); cout<<"6 jets "<<nj6[j]<<endl;}
        else if (k == 8) {nj7.push_back(njets); cout<<"7 jets "<<nj7[j]<<endl;}
        else if (k == 9) {nj8.push_back(njets); cout<<"8 jets "<<nj8[j]<<endl;}
        else if (k == 10) {nj9.push_back(njets); cout<<"9 jets "<<nj9[j]<<endl;}
        else if (k == 11) {nj10.push_back(njets); cout<<"10 jets "<<nj10[j]<<endl;}
        else if(k > 11){njmore.push_back(njets); cout<<"more "<< k<<" "<<njmore[j]<<endl;}
       } // close for bins
       ntot.push_back(nj3[j]+nj4[j]+nj5[j]+nj6[j]+nj7[j]+nj8[j]+nj9[j]+nj10[j]+njmore[j]);
       cout<<"total "<<ntot[j]<<endl;
   } // close for masses
   TMultiGraph *mg = new TMultiGraph();
   mg->SetMaximum(2);
   int nevents=1.;//20000;
   //
   TGraphErrors *TreeJets = new TGraphErrors(1);
   TreeJets->SetMarkerStyle(kFullDotLarge); 
   TreeJets->SetLineColor(kYellow+2);
   TreeJets->SetMarkerColor(kYellow+2);
   TreeJets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) TreeJets->SetPoint(j, masses[j], nj3[j]/nevents);
   mg->Add(TreeJets,"L,P");//->Draw("L,P");
   leg3->AddEntry(TreeJets, " 3 jets", "LP");
   //
   TGraphErrors *fourJets = new TGraphErrors(1);
   fourJets->SetMarkerStyle(kFullDotLarge); 
   fourJets->SetLineColor(8);
   fourJets->SetMarkerColor(8);
   fourJets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) fourJets->SetPoint(j, masses[j], nj4[j]/nevents);
   mg->Add(fourJets,"L,P");//->Draw("L,P");
   leg3->AddEntry(fourJets, " 4 jets", "LP");
   //
   TGraphErrors *fiveJets = new TGraphErrors(1);
   fiveJets->SetMarkerStyle(kFullDotLarge); 
   fiveJets->SetLineColor(kCyan);
   fiveJets->SetMarkerColor(kCyan);
   fiveJets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) fiveJets->SetPoint(j, masses[j], nj5[j]/nevents);
   mg->Add(fiveJets,"L,P");//->Draw("L,P");
   leg3->AddEntry(fiveJets, " 5 jets", "LP");
   //
   TGraphErrors *sixJets = new TGraphErrors(1);
   sixJets->SetMarkerStyle(kFullDotLarge); 
   sixJets->SetLineColor(kBlue);
   sixJets->SetMarkerColor(kBlue);
   sixJets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) sixJets->SetPoint(j, masses[j], nj6[j]/nevents);
   mg->Add(sixJets,"L,P");//->Draw("L,P");
   leg3->AddEntry(sixJets, " 6 jets", "LP");
   //
   TGraphErrors *sevenJets = new TGraphErrors(1);
   sevenJets->SetMarkerStyle(kFullDotLarge); 
   sevenJets->SetLineColor(kMagenta);
   sevenJets->SetMarkerColor(kMagenta);
   sevenJets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) sevenJets->SetPoint(j, masses[j], nj7[j]/nevents);
   mg->Add(sevenJets,"L,P");//->Draw("L,P");
   leg3->AddEntry(sevenJets, " 7 jets", "LP");
   //
   TGraphErrors *eightJets = new TGraphErrors(1);
   eightJets->SetMarkerStyle(kFullDotLarge); 
   eightJets->SetLineColor(kGreen);
   eightJets->SetMarkerColor(kGreen);
   eightJets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) eightJets->SetPoint(j, masses[j], nj8[j]/nevents);
   mg->Add(eightJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(eightJets, " 8 jets", "LP");
   //
   TGraphErrors *nineJets = new TGraphErrors(1);
   nineJets->SetMarkerStyle(kFullDotLarge); 
   nineJets->SetLineColor(50);
   nineJets->SetMarkerColor(50);
   nineJets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) nineJets->SetPoint(j, masses[j], nj9[j]/nevents);
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
   for(int j=0;j<maxtodo-1;j++) moreJets->SetPoint(j, masses[j], njmore[j]/nevents);
   mg->Add(moreJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(moreJets, " > 10 jets", "LP");
   //
   TGraphErrors *totJets = new TGraphErrors(1);
   totJets->SetMarkerStyle(kFullDotLarge); 
   totJets->SetLineColor(kBlack);
   totJets->SetMarkerColor(kBlack);
   totJets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) totJets->SetPoint(j, masses[j], ntot[j]/nevents);
   mg->Add(totJets,"L,P");//->Draw("L,P");
   leg2->AddEntry(totJets, " total", "LP");
   //
  mg->Draw("AP");
  leg2->Draw();
  leg3->Draw();
  PT_HAT->SetLogy(0);
  PT_HAT->SaveAs("njets_flat.png");
  PT_HAT->SaveAs("njets_flat.pdf");
  PT_HAT->SaveAs("njets_flat.root");
  PT_HAT->Clear();
  leg2->Clear();
  leg3->Clear();
//////////////////////////////////////////////////////////////////////
// cat plot
   TMultiGraph *mg1 = new TMultiGraph();
   //mg1->SetMaximum(0.6);
  //mg->GetXaxis()->SetRangeUser(490,520);

   int nevents=1;//20000;
   vector<double> cat0,cat1,cat2,ctot;
    for(int j=0;j<maxtodo;j++) {
       cout<<lege[todo[j]]<<endl;
       // get number of njets
       int nbins = plots[j][26].GetNbinsX();
       //cout<<nbins<<endl;
       for(int k=0;k<nbins;k++) {
	double njets = plots[j][26].GetBinContent(k); 
        if (k == 2)      {cat0.push_back(njets); cout<<"0 tag "<<cat0[j]<<endl;}
        else if (k == 3) {cat1.push_back(njets); cout<<"1 tag "<<cat1[j]<<endl;}
        else if (k == 4) {cat2.push_back(njets); cout<<"2 tag "<<cat2[j]<<endl;}
       } // close for bins
       ctot.push_back(cat0[j]+cat1[j]+cat2[j]);
       cout<<"total "<<ctot[j]<<endl;
   } // close for masses
   for(int j=0;j<maxtodo;j++) cout<<lege[todo[j]]<<" "<<ctot[j]<<" "<<plots[j][26].Integral()<<endl;
   TMultiGraph *mg = new TMultiGraph();
   //
   TGraphErrors *cat0Jets = new TGraphErrors(1);
   cat0Jets->SetMarkerStyle(kFullDotLarge); 
   cat0Jets->SetLineColor(kMagenta);
   cat0Jets->SetMarkerColor(kMagenta);
   cat0Jets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) cat0Jets->SetPoint(j, masses[j], cat0[j]/nevents);
   mg1->Add(cat0Jets,"L,P");//->Draw("L,P");
   mg1.SetMaximum(1.1);
   leg2->AddEntry(cat0Jets, " 0 tag", "LP");
   //
   TGraphErrors *cat1Jets = new TGraphErrors(1);
   cat1Jets->SetMarkerStyle(kFullDotLarge); 
   cat1Jets->SetLineColor(8);
   cat1Jets->SetMarkerColor(8);
   cat1Jets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) cat1Jets->SetPoint(j, masses[j], cat1[j]/nevents);
   mg1->Add(cat1Jets,"L,P");//->Draw("L,P");
   leg2->AddEntry(cat1Jets, " 1 tag", "LP");
   //
   TGraphErrors *cat2Jets = new TGraphErrors(1);
   cat2Jets->SetMarkerStyle(kFullDotLarge); 
   cat2Jets->SetLineColor(kCyan);
   cat2Jets->SetMarkerColor(kCyan);
   cat2Jets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) cat2Jets->SetPoint(j, masses[j], cat2[j]/nevents);
   mg1->Add(cat2Jets,"L,P");//->Draw("L,P");
   leg3->AddEntry(cat2Jets, " 2 tags", "LP");
   //
   TGraphErrors *totJets = new TGraphErrors(1);
   totJets->SetMarkerStyle(kFullDotLarge); 
   totJets->SetLineColor(1);
   totJets->SetMarkerColor(1);
   totJets->SetLineWidth(3);
   for(int j=0;j<maxtodo-1;j++) totJets->SetPoint(j, masses[j], ctot[j]/nevents);
   mg1->Add(totJets,"L,P");//->Draw("L,P");
   leg3->AddEntry(totJets, " total", "LP");
   //
  mg1->Draw("AP");
  mg1->GetXaxis()->SetTitle("M_{Gr} GeV");
  mg1->GetYaxis()->SetTitle("Efficiency");
  leg2->Draw();
  leg3->Draw();
  PT_HAT->SetLogy(0);
  PT_HAT->SaveAs("fatcat.png");
  PT_HAT->SaveAs("fatcat.pdf");
  PT_HAT->SaveAs("fatcat.root");
//PT_HAT->Close(); 

}
