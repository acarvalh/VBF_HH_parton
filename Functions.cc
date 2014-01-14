#include <fstream>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "deltaphi.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TArray.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TArray.h>
#include <TVector.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TObject.h>
#include <TMath.h>
#include "root_plot.h" //declare the histos
using namespace fastjet;
using namespace std;

void hello(){
  cout<<"\n\n\n HELLO!!!! \n\n"<<endl;
}
/////////////////////////////////////////////////////////////////
int fill_hist(
        // vector that keep the higgses constituents
        vector<PseudoJet> H1,
        vector<PseudoJet> H2,
        // bvf jets
        PseudoJet j1, PseudoJet j2
        ){
        PseudoJet Higgs_1 = H1.at(0)+H1.at(1);
        PseudoJet Higgs_2 = H2.at(0)+H2.at(1);
        PseudoJet Xres= Higgs_1 + Higgs_2;
        //
        double DRhh = Higgs_1.delta_R(Higgs_2); // in between Higgs
        //double DRjj=-10, deltarap=-10;
        //if(!fattag) {DRjj= j1.delta_R(j2); deltarap= abs(higgs_b.eta()-Higgs_ph.eta());}// in between the jets
        //
        const int numvar=25;
        double monovar[numvar] = {
          Higgs_1.m(), Higgs_1.e(), Higgs_1.pt(), Higgs_1.eta(), Higgs_1.phi(), //5
//	  0,0,0,0,0, //5
	  0, //1
	  Higgs_2.m(),Higgs_2.e(),Higgs_2.pt(),Higgs_2.eta(),Higgs_2.phi(), // 5
	  Xres.m(),Xres.e(),Xres.pt(),Xres.eta(),Xres.phi(), // 5
//	  0,0,0,0,0, //5
//	  0,0,0,0,0, //5
	  0,  //DPhiphph, // 2
	  DRhh, //DPhihh, // 2
//	  0,
	  0,// DPhijj, // 2
	  0,0,0,0,0
    }; //25
  for (int j=0; j<numvar; j++) basic[j]->Fill(monovar[j],1);
  return 1;
}
/////////////////////////////////////////////////////////////////////////
// save the histos
int save_hist(int nmass){
const char* Mass = Form("histos/Control_shower_%d.root",nmass);
TFile f1(Mass, "recreate");
//cout<<basic.size()<<endl;
f1.cd();
//cout<<nmass<<endl;
Njets_passing_kLooseID->Write();
//Sel->Write();
//cout<<nmass<<endl;
for (int i=0; i<basic.size(); i++) {basic[i]->Write(); }
f1.Close();
for (int i=0; i<24; i++) basic[i]->Reset();

return 1;
}
///////////////////////////////////////////////////////////////////////////
// declare the histos
int decla(int mass){
   const char* label="kk graviton ";

	Njets_passing_kLooseID = new TH1F("njets_passing_kLooseID_ct4",  
		label, 
		11, -1.5, 9.5);
	Njets_passing_kLooseID->GetYaxis()->SetTitle("");
	Njets_passing_kLooseID->GetXaxis()->SetTitle("Njets after showering");

	TH1F *PPhotonsMass = new TH1F("PhotonsMass_ct4",  
		label, 
		1000, 95, 185);
	PPhotonsMass->GetYaxis()->SetTitle("Events/ 2 GeV");
	PPhotonsMass->GetXaxis()->SetTitle("M_{#gamma #gamma} (GeV)");
	basic.push_back (PPhotonsMass); //69
	//
	TH1F *Dipho_E = new TH1F("dipho_E_ct4",  
		label, 
		40, 90, 800);
	Dipho_E->GetYaxis()->SetTitle("Events");
	Dipho_E->GetXaxis()->SetTitle("E_{#gamma #gamma} (GeV)");
	basic.push_back (Dipho_E); //70
	//
	TH1F *Dipho_pt = new TH1F("dipho_pt_ct4",  
		label, 
		25, -1, 510);
	Dipho_pt->GetYaxis()->SetTitle("Events/ 2 GeV");
	Dipho_pt->GetXaxis()->SetTitle("Pt_{#gamma #gamma} (GeV)");
	basic.push_back (Dipho_pt); //71
	//
	TH1F *Dipho_eta = new TH1F("dipho_eta_ct4",  
		label, 
		35, -9, 9);
	Dipho_eta->GetYaxis()->SetTitle("Events/ 2 GeV");
	Dipho_eta->GetXaxis()->SetTitle("#eta_{#gamma #gamma}");
	basic.push_back (Dipho_eta); //72
	//
	TH1F *Dipho_phi = new TH1F("dipho_phi_ct4",  
		label, 
		25, -4.14, 7.14);
	Dipho_phi->GetYaxis()->SetTitle("Events/ 2 GeV");
	Dipho_phi->GetXaxis()->SetTitle("#phi_{#gamma #gamma}");
	basic.push_back (Dipho_phi); //73
	//
	// dijet variables
	//
	// "JetsMass","dijet_E","dijet_Pt","dijet_Eta","dijet_Phi",
	TH1F *Njets_passing_kLooseID2 = new TH1F("njets_passing_kLooseID_ct42",  
		label, 
		11, -1.5, 9.5);
	Njets_passing_kLooseID2->GetYaxis()->SetTitle("");
	Njets_passing_kLooseID2->GetXaxis()->SetTitle("Njets after reconstruction");
	basic.push_back (Njets_passing_kLooseID2);
	//
	TH1F *JJetsMass = new TH1F("JetsMass_ct4",  
		label, 
		50, 0, 400);
	JJetsMass->GetXaxis()->SetTitle("M_{jj} (GeV)");
	basic.push_back (JJetsMass);
	//
	TH1F *Dijet_E = new TH1F("dijet_E_ct4",  
		label, 
		40, 0, 600);
	Dijet_E->GetYaxis()->SetTitle("Events/ 2 GeV");
	Dijet_E->GetXaxis()->SetTitle("E_{jj} (GeV)");
	basic.push_back (Dijet_E);
	//
	TH1F *Dijet_pt = new TH1F("dijet_pt_ct4",  
		label, 
		60, 0, 700);
	Dijet_pt->GetXaxis()->SetTitle("Pt_{jj} (GeV)");
	basic.push_back (Dijet_pt);
	//
	TH1F *Dijet_eta = new TH1F("dijet_eta_ct4",  
		label, 
		35, -9, 9);
	Dijet_eta->GetYaxis()->SetTitle("Events/ 2 GeV");
	Dijet_eta->GetXaxis()->SetTitle("#eta_{jj}");
	basic.push_back (Dijet_eta);
	//
	TH1F *Dijet_phi = new TH1F("dijet_phi_ct4",  
		label, 
		25, -4.14, 7.14);
	Dijet_phi->GetYaxis()->SetTitle("");
	Dijet_phi->GetXaxis()->SetTitle("#phi_{jj}");
	basic.push_back (Dijet_phi);
	//
	// radion variables
	//
	//  "RadMass", "radion_E", "radion_Pt", "radion_Eta","radion_Phi",
	TH1F *RRadMass = new TH1F("RadMass_ct4",  
		label, 
		50, 50, 2500);
	RRadMass->GetXaxis()->SetTitle("M_{#gamma #gamma j j } (GeV)");
	basic.push_back (RRadMass);
	//
	TH1F *Radion_E = new TH1F("radion_E_ct4",  
		label, 
		40, 150, 2000);
	Radion_E->GetXaxis()->SetTitle("E_{##gamma #gamma j j } (GeV)");
	basic.push_back (Radion_E);
	//
	TH1F *Radion_pt = new TH1F("radion_pt_ct4",  
		label, 
		40, 0, 700);
	Radion_pt->GetYaxis()->SetTitle("");
	Radion_pt->GetXaxis()->SetTitle("Pt_{#gamma #gamma j j } (GeV)");
	basic.push_back (Radion_pt);
	//
	TH1F *Radion_eta = new TH1F("radion_eta_ct4",  
		label, 
		35, -9 , 9);
	Radion_eta->GetYaxis()->SetTitle("");
	Radion_eta->GetXaxis()->SetTitle("#eta_{#gamma #gamma j j }");
	basic.push_back (Radion_eta);
	//
	TH1F *Radion_phi = new TH1F("radion_phi_ct4",  
		label, 
		25, -4.14, 7.14);
	Radion_phi->GetYaxis()->SetTitle("");
	Radion_phi->GetXaxis()->SetTitle("#phi_{#gamma #gamma j j }");
	basic.push_back (Radion_phi);
	//
	// distances no var on tree
	//
	TH1F *DRphph = new TH1F("DRphph_ct4",  
		label, 
		30, -2, 9);
	DRphph->GetYaxis()->SetTitle("");
	DRphph->GetXaxis()->SetTitle("#Delta R( #gamma #gamma )");
	basic.push_back (DRphph);
	//
	TH1F *DRhh = new TH1F("DRhh_ct4",  
		label, 
		50, -2, 9);
	DRhh->GetYaxis()->SetTitle("");
	DRhh->GetXaxis()->SetTitle("#Delta R( H H )");
	basic.push_back (DRhh);
	//
	TH1F *DRjj = new TH1F("DRjj_ct4",  
		label, 
		50, -2, 9);
	DRjj->GetYaxis()->SetTitle("");
	DRjj->GetXaxis()->SetTitle("#Delta R( jj )");
	basic.push_back (DRjj);
	//
	TH1F *AAass = new TH1F("aaass_ct4",  
		label, 
		20, 0, 6);
	AAass->GetYaxis()->SetTitle("");
	AAass->GetXaxis()->SetTitle("delta pseudorapidity #gamma");
	basic.push_back (AAass);
	//
	TH1F *JJass = new TH1F("jjass_ct4",  
		label, 
		4, -0.5, 3.5);
	JJass->GetYaxis()->SetTitle("");
	JJass->GetXaxis()->SetTitle("generated jets selected");
	basic.push_back (JJass);
	//
	TH1F *leadsel = new TH1F("leadsel_ct4",  
		label, 
		10, -0.5, 9.5);
	leadsel->GetYaxis()->SetTitle("");
	leadsel->GetXaxis()->SetTitle("leading photon selected resolved");
	basic.push_back (leadsel);
	//
	TH1F *subleadsel = new TH1F("sublead_ct4",  
		label, 
		10, -0.5, 9.5);
	subleadsel->GetYaxis()->SetTitle("");
	subleadsel->GetXaxis()->SetTitle("subleading photon selected resolved");
	basic.push_back (subleadsel);
	//
	TH1F *fattagreco = new TH1F("fattagreco_ct4",  
		label, 
		5, -0.5, 4.5);
	fattagreco->GetYaxis()->SetTitle("");
	fattagreco->GetXaxis()->SetTitle("fat jet selected");
	basic.push_back (fattagreco);
	//
	TH1F *fattagrecon = new TH1F("fattagrecon_ct4",  
		label, 
		5, -0.5, 4.5);
	fattagrecon->GetYaxis()->SetTitle("");
	fattagrecon->GetXaxis()->SetTitle("# of fat jet reconstructed");
	basic.push_back (fattagrecon);
        //basic[0][0][0][2].Rebin(200, mass-100, mass+100);
return 1;
}


