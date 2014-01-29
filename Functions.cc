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
#include "cuts.h" // basic cuts 
using namespace fastjet;
using namespace std;

void hello(){cout<<"\n\n\n HELLO!!!! \n\n"<<endl;}
/////////////////////////////////////////////////////////////////
bool analyse_4b(
	vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag,vector<int> vbftag){
  // pair the jets that was not VBF tagged
  bool found = false;
  // now we separate analysis
  // number of fattags
  int nfat=0, nbtag=0, nmistag=0;
  for(int i=0;i<jets.size();i++) 
	{nfat = nfat + fattag[i]; nbtag = nbtag + btag[i]; nmistag = nmistag + bmistag[i];}
  PseudoJet H1,H2;
  float Hmin = higgs_mass*(1-tolerance);
  float Hmax = higgs_mass*(1+tolerance);
  if(jets.size() > 3 && nfat >1) { // close if 1 tag
    //std::cout<<"2 tag! "<<std::endl;
    //if(nfat ==2){
      H1=jets.at(fattag[0]);
      H2=jets.at(fattag[1]);
    //} //else found= false; // if more than 2 thinhk 
    // quality requirements   
    double massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
    double rapDiff = abs(H1.eta() - H2.eta());
    if( 
       massDiff < tolerance && //rapDiff < deltaEtaHH &&
       (H1.m() > Hmin && H1.m() < Hmax) &&
       (H2.m() > Hmin && H2.m() < Hmax) 
       && (H1+H2).m() >MHH
       && rapDiff<DetaHH
      ){  found=true; Cat->Fill(2.,weight);}   
  } // close 1 tag
   else if(nfat >0 && jets.size() > 4 &&  found==false) { 
   //std::cout<<"1 tag! "<<std::endl;
   if(nfat>1) cout<<"ops! "<<endl;
   H1=jets.at(fattag[0]);
   // pair H2 the jets by the minimum invariant mass difference with H1
   std::vector<double> a1; 
   std::vector< int > jetn1, jetn2; // to keep the pairs
   double invmassB =  H1.m();
   for(int nj1=0; nj1< jets.size(); nj1++)  if(nj1 != vbftag[0] && nj1 != vbftag[1]) 
     for(int nj2=nj1+1; nj2< jets.size(); nj2++) if(nj2 != vbftag[0] && nj2 != vbftag[1]) { 
	   //std::cout<<nj1<<nj2<<" "<<nj3<<nj4<<std::endl;
	   double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
	   a1.push_back((invmassA-invmassB)*(invmassA-invmassB)); 
	   jetn1.push_back(nj1);jetn2.push_back(nj2); // we also what to keep the nj...           
    } // loop on jets
    int minM;
    //Find the minumum value of the vector (iterator version)
    minM = TMath::LocMin(a1.size(), &a1[0]);
    //std::cout<<"hi, the jets pairs are !!!! "<<jetn1[minM]<<jetn2[minM]<<" "
	//<<fattag[0]<<std::endl;
    H2=jets.at(jetn1[minM])+jets.at(jetn2[minM]);
    // quality requirements   
    double massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
    double rapDiff = abs(H1.eta() - H2.eta());
    if( 
       massDiff < tolerance && //rapDiff < deltaEtaHH &&
       (H1.m() > Hmin && H1.m() < Hmax) &&
       (H2.m() > Hmin && H2.m() < Hmax)
       && (H1+H2).m() >MHH
       && rapDiff<DetaHH
       && abs(jets.at(jetn1[minM]).eta()+jets.at(jetn2[minM]).eta())<DetaH
      ){ found=true; Cat->Fill(1.,weight);} 
  } // close if 1 tag
  else if(jets.size() > 5 && found==false) { // resolved
   // pair the jets by the minimum invariant mass difference
   std::vector<double> a1,a2; 
   //std::cout<<"resolved! "<<std::endl;
   // pair by minimum distance of higgs of one pair:test
   std::vector< int > jetn1, jetn2,jetn3, jetn4; // to keep the pairs
   for(int nj1=0; nj1< jets.size(); nj1++)   
     for(int nj2=0; nj2< jets.size(); nj2++) 
       for(int nj3=0; nj3< jets.size(); nj3++)     
	 for(int nj4=0; nj4< jets.size(); nj4++)   
           { 
           if( 1>0
           && nj1 != vbftag[0] && nj1 != vbftag[1]
	   && nj2 != vbftag[0] && nj2 != vbftag[1] 
           && nj3 != vbftag[0] && nj3 != vbftag[1]
           && nj4 != vbftag[0] && nj4 != vbftag[1] 
           && nj1 !=nj2 && nj1 !=nj3 && nj1 !=nj4 
           && nj2 !=nj3 && nj2 !=nj4 
           && nj1 !=nj4  
              ){
	      double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
	      double invmassB =  (jets.at(nj3)+jets.at(nj4)).m();
	      //if(invmassA> 100 && invmassB> 100) {
	      //std::cout<<nj1<<nj2<<" "<<nj3<<nj4<<std::endl;
              a1.push_back((invmassA-invmassB)*(invmassA-invmassB)); 
              //cout << " higgs masses "<<invmassA<<" "<<invmassB <<endl;
              //cout << " higgs masses wrong "<<(jets.at(nj1)+jets.at(nj3)).m()<<" "<<(jets.at(nj4)+jets.at(nj2)).m() <<endl;
	      jetn1.push_back(nj1);jetn2.push_back(nj2); // we also what to keep the nj...
	      jetn3.push_back(nj3);jetn4.push_back(nj4);
           }// close if not btagged
    } // loop on jets
    int minM;
    //Find the minumum value of the vector (iterator version)
    minM = TMath::LocMin(a1.size(), &a1[0]);
    //std::cout<<"hi, the jets pairs are !!!! "<<jetn1[minM]<<jetn2[minM]<<" "
	//<<jetn3[minM]<<jetn4[minM]<<std::endl;
    H1=jets.at(jetn1[minM])+jets.at(jetn2[minM]);
    H2=jets.at(jetn3[minM])+jets.at(jetn4[minM]);  
    // quality requirements   
    double massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
    double rapDiff = abs(H1.eta() - H2.eta());
    if( 
       massDiff < tolerance && //rapDiff < deltaEtaHH &&
       (H1.m() > Hmin && H1.m() < Hmax) &&
       (H2.m() > Hmin && H2.m() < Hmax)
       && (H1+H2).m() >MHH
       && rapDiff<DetaHH
       && abs(jets.at(jetn1[minM]).eta()+jets.at(jetn2[minM]).eta())<DetaH
       && abs(jets.at(jetn3[minM]).eta()+jets.at(jetn4[minM]).eta())<DetaH
      ){ //std::cout<<"getting there"<<std::endl; 
	found=true; Cat->Fill(0.,weight);} 
    //else if(massDiff > tolerance) 
    //cout << "2 tag failed , btags = "<<btag[jetn1[minM]]<<" "<<btag[jetn2[minM]]<<" "
//			<<btag[jetn3[minM]]<<" "<<btag[jetn4[minM]]<<" " <<endl;
//    cout << "               higgs masses "<<H1.m()<<" "<<H2.m() <<endl;
  } //else cout<<"njets"<<jets.size()<<" nfat "<< nfat<<endl; // close if resolved
  //////////////////////////////////////
  if(found){
       ///////////////////////////
       // fill the histos
       ///////////////////////////
       //cout<<"fat tag = " <<nfat <<" number of plots "<<basicHiggses.size()<<endl;
       const int numvar=16;
       PseudoJet Xres = H1 +H2; 
       double monovar[numvar] = {
        H1.m(),H1.pt(),H1.eta(),H1.phi(), //4
        H2.m(),H2.pt(),H2.eta(),H2.phi(), // 4
        Xres.m(),Xres.pt(),Xres.eta(),Xres.phi(), // 4
	nfat,nbtag,nmistag,abs(H1.eta() - H2.eta())
       }; //25
       for (int j=0; j<basicHiggses.size(); j++) basicHiggses[j]->Fill(monovar[j],weight); // weight
  } // close if fill 
return found; // close if 2 tags
} // close 4b analysis
/////////////////////////////////////////////////////////////////
bool findVBF(vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag, vector<int> & vbftag){
  // highest invariant mass among the non-tagged
  vector<int> listNonTag; // list the non tagged
  for(int nj1=0; nj1< jets.size(); nj1++) { 
  //cout<<"fat tagged = " <<fattag[nj1]<<" b-quarks mistagged = " <<bmistag[nj1] <<" b-quark = " <<btag[nj1] <<endl;
	if(fattag[nj1]==0 && btag[nj1] ==0 && bmistag[nj1]==0
          ) {listNonTag.push_back(nj1);}
  }
  if(listNonTag.size()>1){   // find the hightest inv mass pair 
    int nfat=0, nbtag=0, nmistag=0;
    for(int i=0;i<listNonTag.size();i++) 
	{nfat = nfat + fattag[listNonTag[i]]; 
         nbtag = nbtag + btag[listNonTag[i]]; nmistag = nmistag + bmistag[listNonTag[i]];}
    std::vector<double> a1; 
    std::vector< int > jetn1, jetn2; // to keep the pairs
    for(int nj1=0; nj1< listNonTag.size(); nj1++) 
	for(int nj2=nj1+1; nj2< listNonTag.size(); nj2++) { // we also what to keep the nj...
	  double invmass =  (jets.at(listNonTag[nj1])+jets.at(listNonTag[nj2])).m();
	  a1.push_back(invmass); jetn1.push_back(listNonTag[nj1]);jetn2.push_back(listNonTag[nj2]);
    } // loop on jets  
    //
    int i1; i1 = TMath::LocMax(a1.size(), &a1[0]); // max inv mass
    vbftag.push_back(jetn1[i1]); vbftag.push_back(jetn2[i1]); // save the pair number
    double etaVBF = abs(jets.at(vbftag[0]).eta()-jets.at(vbftag[1]).eta());
    if( 1>0 && a1[i1] > HTVBF && etaVBF > DeltayVBF 
        && jets.at(vbftag[0]).pt()>jet_ptminvbf
        && jets.at(vbftag[1]).pt()>jet_ptminvbf){ // apply the VBF cuts
       // std::cout<<"hi VBF jets really are !!!! "<<vbftag[0]<<" "<<vbftag[1]<<std::endl;
       //    cout << " dijet mass "<<(jets.at(vbftag[0])+jets.at(vbftag[1])).m() <<endl;
       //  cout << " dijet btag "<<btag[vbftag[0]]<<btag[vbftag[1]] <<endl;
       ///////////////////////////
       // fill the histos
       ///////////////////////////
       const int numvar=9;
       PseudoJet vbfmass = jets.at(vbftag[0]) +jets.at(vbftag[1]); 
       double Deta = abs(jets.at(vbftag[0]).eta() -jets.at(vbftag[1]).eta());
       //cout<<"plots vbf "<<basicvbf.size()<<endl;
       double monovarvbf[numvar] = {
        jets.at(vbftag[0]).pt(),jets.at(vbftag[0]).eta(), //2
        jets.at(vbftag[1]).pt(),jets.at(vbftag[1]).eta(), //2
        vbfmass.m(),Deta, // 2
	nfat,nbtag,nmistag
       }; //25
       for (int j=0; j<basicvbf.size(); j++) basicvbf[j]->Fill(monovarvbf[j],weight); // weight
    return true;
  } else return false; // close VBF cuts
 } else{ //cout<<"not enought vbf ! "<<listNonTag.size()<<endl; 
         return false;}// close if listnontag>1 
} // close VBF selection
/////////////////////////////////////////////////////////////////
void isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag);
/////////////////////////////
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets_akt, vector<int> & btag, vector<int> & bmistag, vector<int> & fattag){
  JetDefinition akt(antikt_algorithm, 0.5);
  ClusterSequence cs_akt(particles, akt);
  //vector<PseudoJet> jets_akt;
  Selector jet_selector = SelectorPtMin(jet_ptmin) && SelectorAbsRapMax(rapmax);
  if(shower){
    // first we do akt jets from particles
    jets_akt = sorted_by_pt(jet_selector(cs_akt.inclusive_jets()));
  }
  else{
    double const ptmin=jet_ptmin;
    Selector jet_selector_parton = SelectorPtMin(ptmin);
    jets_akt = sorted_by_pt(jet_selector_parton(cs_akt.inclusive_jets()));
/*
//////////////////
for (int i=0; i<jets_akt.size(); i++){
  vector<PseudoJet> constitu=jets_akt.at(i).constituents();
  cout<<"constituents size "<<constitu.size()<<endl;
  for (int j=0; j<constitu.size(); j++) {
    //int test=constitu.at(j).user_index();
    cout<<"constituents flavour "<< constitu.at(j).user_index()<<endl;
  }
}
//////////////////
*/
//  if(jets_akt.size()==5) {
//   cout<<"------------------------"<<endl;
//   for (int i=0; i<jets_akt.size(); i++) if(jets_akt.at(i).constituents().size() ==2 && jets_akt.at(i).m()<120) cout<< jets_akt.at(i).m()<<" "<<" "<<endl;
//  }
  }
  int njets = jets_akt.size();
  Njets_passing_kLooseID->Fill(njets,weight);
  isbtagged(jets_akt, btag, bmistag); // check wheather the b(c)jet is b--(mis)tagable
  /////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////// check tag
  JetDefinition CA10(cambridge_algorithm, Rsb);
  // Filter definition to improve mass resolution
  Filter filter(JetDefinition(cambridge_algorithm, Rfilt), SelectorNHardest(n_subjet));
  PseudoJet tagged_jet;
  for (int i = 0; i < jets_akt.size(); i++) {
    int see=0;
    // first recluster with some large CA (needed for mass-drop)
    ClusterSequence cs_tmp(jets_akt[i].constituents(), CA10);
    // next get hardest jet
    PseudoJet ca_jet = sorted_by_pt(cs_tmp.inclusive_jets())[0]; // find the cores
    // now run mass drop tagger
    MassDropTagger md_tagger(mu, ycut); // define the cut on mass drop
	// mu: ratio in between mass of cores, symetric splitting
    tagged_jet = md_tagger(ca_jet);
    if(tagged_jet.m()>0) {
      PseudoJet filtered_jet = filter(jets_akt.at(i)); // filter to tag
      if(filtered_jet.m()>Mfat) see = 1; else see = 0; 
    } else see = 0; 
    fattag.push_back(see);
  } // close find mass drop
  ////////////////////////////////////////////////////////////////////////////////////////////////
  //jets = jets_akt;
  return njets;
} // close cluster jets
/////////
void isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag){ 
  for (int i=0; i<jets.size(); i++) { // check wheter jet have inside a b's are taggable 
  int see=0,see2=0;
     vector<PseudoJet> constitu=jets.at(i).constituents();
     for (int j=0; j<constitu.size(); j++) {
       //cout<<"constituents flavour "<<constitu.at(j).user_index()<<endl;
       if((constitu.at(j).user_index() == 5 || constitu.at(j).user_index() == -5 ) // work !!
            && constitu.at(j).pt() > 10
	    && constitu.at(j).eta() < etab
          ) {see++;}//btag.push_back(1);} else btag.push_back(0); 
       if( abs(constitu.at(j).user_index()) == 4  // work !! 
	    && constitu.at(j).pt() > 10
         ) {see2++;}// bmistag.push_back(1);} bmistag.push_back(0);
     } // close constituents
     //bmistag.push_back(see2);
     if(see>0) btag.push_back(1); else btag.push_back(0); // count only one tag/jet
     if(see2>0) bmistag.push_back(1); else bmistag.push_back(0); // count only one tag/jet
     //cout<<"b-quarks mistagged = " <<bmistag[i] <<" b-quark = " <<btag[i] <<endl;
  } // close for each jet


} // close isbtagged
//////////////////////////////////////////////////////////////////////////////////////////////
void istagged(vector<PseudoJet> jets, vector<int> & fattag){
//  for (int i=0; i<jets.size(); i++) { // check wheter b's are taggable 
//     int see=0;
//     if( jets.at(i).m() > 100) see = 1; else see = 0; 
//     fattag.push_back(see);
//  } // close for each jets
  // do not filter yet
//  if(jets.size()==5) {
//   cout<<"------------------------"<<endl;
//   for (int i=0; i<jets.size(); i++) if(fattag[i]==0) cout<< jets.at(i).m()<<" "<<endl;
//  }
} // close istagged
/////////////////////////////////////////////////////////////////////////
// save the histos
int save_hist(int nmass, bool resonant,bool bkg){
  const char* Mass;
  if(resonant && !bkg) {
    Mass = Form("histos/Control_shower_%d.root",nmass);
    //Mass = Form("histos/Madgraph0_0137/Control_shower_%d.root",nmass);
  }
  else if(!bkg) Mass = Form("histosnonres/Control_shower_%d.root",nmass);
  else Mass = Form("4bsbkg/Control_shower_%d.root",nmass);
  TFile f1(Mass, "recreate");
  f1.cd();
  Njets_passing_kLooseID->Write();
  Cat->Write();
  for (int j=0; j<basicHiggses.size(); j++){basicHiggses[j]->Write();}
  for (int i=0; i<basicvbf.size(); i++) {basicvbf[i]->Write();}
  f1.Close();
  //
  Njets_passing_kLooseID->Reset();
  Cat->Reset();
  for (int i=0; i<basicHiggses.size(); i++) basicHiggses[i]->Reset();
  for (int i=0; i<basicvbf.size(); i++) basicvbf[i]->Reset();
  return 1;
}
///////////////////////////////////////////////////////////////////////////
// declare the histos
int decla(int mass){
   const char* label="kk graviton ";

	Njets_passing_kLooseID = new TH1F("njets_passing_kLooseID_ct4",  
		label, 
		13, -0.5, 12.5);
	Njets_passing_kLooseID->GetYaxis()->SetTitle("");
	Njets_passing_kLooseID->GetXaxis()->SetTitle("Njets after showering");

	TH1F *H1hist = new TH1F("H1hist",  
		label, 
		90, 95, 185);
	H1hist->GetYaxis()->SetTitle("Events/ 2 GeV");
	H1hist->GetXaxis()->SetTitle("M_{H1} (GeV)");
	basicHiggses.push_back (H1hist); 

	TH1F *H1histpt = new TH1F("H1histpt",  
		label, 
		120, 0, 700);
	H1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	H1histpt->GetXaxis()->SetTitle("H1 P_T (GeV)");
	basicHiggses.push_back (H1histpt); 

	TH1F *H1histeta = new TH1F("H1histeta",  
		label, 
		30, -6, 6);
	H1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	H1histeta->GetXaxis()->SetTitle("#eta_{H1} (GeV)");
	basicHiggses.push_back (H1histeta); 

	TH1F *H1histphi = new TH1F("H1histphi",  
		label, 
		30, 0, 5);
	H1histphi->GetYaxis()->SetTitle("Events/ 2 GeV");
	H1histphi->GetXaxis()->SetTitle("#phi_{H1} (GeV)");
	basicHiggses.push_back (H1histphi); 

	//

	TH1F *H2hist = new TH1F("H2hist",  
		label, 
		90, 95, 185);
	H2hist->GetYaxis()->SetTitle("Events/ 2 GeV");
	H2hist->GetXaxis()->SetTitle("M_{H2} (GeV)");
	basicHiggses.push_back (H2hist); 

	TH1F *H2histpt = new TH1F("H2histpt",  
		label, 
		120, 0, 700);
	H2histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	H2histpt->GetXaxis()->SetTitle("H2 P_T (GeV)");
	basicHiggses.push_back (H2histpt); 

	TH1F *H2histeta = new TH1F("H2histeta",  
		label, 
		30, -6, 6);
	H2histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	H2histeta->GetXaxis()->SetTitle("#eta_{H2} (GeV)");
	basicHiggses.push_back (H2histeta); 

	TH1F *H2histphi = new TH1F("H2histphi",  
		label, 
		30, 0, 5);
	H2histphi->GetYaxis()->SetTitle("Events/ 2 GeV");
	H2histphi->GetXaxis()->SetTitle("#phi_{H2} (GeV)");
	basicHiggses.push_back (H2histphi); 

	//

	TH1F *RRadMass = new TH1F("RadMass_ct4",  
		label, 
		90, 50, 5500);
	RRadMass->GetXaxis()->SetTitle("M_{ HH } (GeV)");
	basicHiggses.push_back (RRadMass);
	
	TH1F *Radion_pt = new TH1F("radion_pt_ct4",  
		label, 
		50, 0, 700);
	Radion_pt->GetYaxis()->SetTitle("");
	Radion_pt->GetXaxis()->SetTitle("Pt_{HH} (GeV)");
	basicHiggses.push_back (Radion_pt);
	
	TH1F *Radion_eta = new TH1F("radion_eta_ct4",  
		label, 
		35, -9 , 9);
	Radion_eta->GetYaxis()->SetTitle("");
	Radion_eta->GetXaxis()->SetTitle("#eta_{HH }");
	basicHiggses.push_back (Radion_eta);
	//
	TH1F *Radion_phi = new TH1F("radion_phi_ct4",  
		label, 
		25, 0, 7.14);
	Radion_phi->GetYaxis()->SetTitle("");
	Radion_phi->GetXaxis()->SetTitle("#phi_{HH}");
	basicHiggses.push_back (Radion_phi);

	//

	TH1F *Nfattag = new TH1F("fattag_ct4",  
		label, 
		6, -1.5, 4.5);
	Nfattag->GetYaxis()->SetTitle("");
	Nfattag->GetXaxis()->SetTitle("Number of fat tags");
	basicHiggses.push_back (Nfattag);

	TH1F *Nbtag = new TH1F("btag_ct4",  
		label, 
		9, -1.5, 7.5);
	Nbtag->GetYaxis()->SetTitle("");
	Nbtag->GetXaxis()->SetTitle("Number b-jets");
	basicHiggses.push_back (Nbtag);

	TH1F *Nbmistag = new TH1F("bmistag_ct4",  
		label, 
		7, -1.5, 5.5);
	Nbmistag->GetYaxis()->SetTitle("");
	Nbmistag->GetXaxis()->SetTitle("Number mistagged b jets's");
	basicHiggses.push_back (Nbmistag);

	TH1F *DetaHH = new TH1F("detahh_ct4",  
		label, 
		30, -1.5, 10.5);
	DetaHH->GetYaxis()->SetTitle("");
	DetaHH->GetXaxis()->SetTitle("#Delta #eta_{HH}");
	basicHiggses.push_back (DetaHH);

	Cat = new TH1F("cat_ct4",  
		label, 
		7, -1.5, 5.5);
	Cat->GetYaxis()->SetTitle("");
	Cat->GetXaxis()->SetTitle("Tag category");
	//basicHiggses.push_back (Cat);
	/////////////////////////////////////////////////////////////////////

	TH1F *j1histpt = new TH1F("j1histpt",  
		label, 
		60, 0, 300);
	j1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	j1histpt->GetXaxis()->SetTitle("vbf j1 P_T (GeV)");
	basicvbf.push_back (j1histpt); 

	TH1F *j1histeta = new TH1F("j1histeta",  
		label, 
		32, -6, 6);
	j1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	j1histeta->GetXaxis()->SetTitle("vbf #eta_{j1} (GeV)");
	basicvbf.push_back (j1histeta); 

	//

	TH1F *j2histpt = new TH1F("j2histpt",  
		label, 
		60, 0, 300);
	j2histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	j2histpt->GetXaxis()->SetTitle("vbf j2 P_T (GeV)");
	basicvbf.push_back (j2histpt); 

	TH1F *j2histeta = new TH1F("j2histeta",  
		label, 
		32, -6, 6);
	j2histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	j2histeta->GetXaxis()->SetTitle("vbf #eta_{j2} (GeV)");
	basicvbf.push_back (H2histeta); 

	//

	TH1F *RRadMassvbf = new TH1F("RadMass_ct4vbf",  
		label, 
		90, 50, 5000);
	RRadMassvbf->GetXaxis()->SetTitle("vbf M_{j j} (GeV)");
	basicvbf.push_back (RRadMassvbf);

	//

	TH1F *Detavbf = new TH1F("Deta",  
		label, 
		50, 0, 25);
	Detavbf->GetXaxis()->SetTitle("vbf #Delta #eta");
	basicvbf.push_back (Detavbf);
	
	//

	TH1F *Nfattagvbf = new TH1F("fattag_ct4vbf",  
		label, 
		6, -1.5, 4.5);
	Nfattagvbf->GetYaxis()->SetTitle("");
	Nfattagvbf->GetXaxis()->SetTitle("Number of fat tags vbf");
	basicvbf.push_back (Nfattagvbf);

	TH1F *Nbtagvbf = new TH1F("btag_ct4vbf",  
		label, 
		7, -1.5, 5.5);
	Nbtagvbf->GetYaxis()->SetTitle("");
	Nbtagvbf->GetXaxis()->SetTitle("Number b-jets vbf");
	basicvbf.push_back (Nbtagvbf);

	TH1F *Nbmistagvbf = new TH1F("bmistag_ct4vbf",  
		label, 
		7, -1.5, 5.5);
	Nbmistagvbf->GetYaxis()->SetTitle("");
	Nbmistagvbf->GetXaxis()->SetTitle("Number mistagged b jets's vbf");
	basicvbf.push_back (Nbmistagvbf);

return 1;
}


