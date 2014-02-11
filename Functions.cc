#include <fstream>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "deltaphi.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TArray.h>
#include <TH1D.h>
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
void genhiggs(int counterh, vector<PseudoJet> higgses){
double hhmass = (higgses.at(0)+higgses.at(1)).m(); 
gen_higgs->Fill(hhmass);
}
/////////////////////////////////////////////////////////////////
bool analyse_2b2w(
	vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag,vector<int> vbftag,
	vector<PseudoJet> leptons, vector<PseudoJet> neutrinos){
  PseudoJet H1,H2; 
  float Hmin = higgs_mass*(1-tolerance);
  float Hmax = higgs_mass*(1+tolerance);
  bool found = false;
  PseudoJet Xres;
  // we first find the leptonic higgs -- do not consider more than 2 leptons
  // do a isolation vector
  vector<double> LepIso;
  int minM;
  std::vector< int > jetn1, jetn2; // to keep the pair
  LepIso.push_back(leptons.at(0).delta_R(leptons.at(1)));
  for(int i = 0;i<jets.size();i++) {
	LepIso.push_back(leptons.at(0).delta_R(jets.at(i)));
	LepIso.push_back(leptons.at(1).delta_R(jets.at(i)));
  } // close iso vector
  int MinDRLep = TMath::LocMin(LepIso.size(), &LepIso[0]);
  //if(MinDRLep<0) for(int i = 0;i<LepIso.size();i++)
  //cout <<LepIso[MinDRLep]<< endl;   
  double MET=(neutrinos.at(0)+neutrinos.at(1)).m();
  if(1>0
   && LepIso[MinDRLep] > lepiso  
   && leptons.at(0).pt()> ptlepton 
   && leptons.at(1).pt()> ptlepton
   && leptons.at(0).eta()< etab 
   && leptons.at(1).eta()< etab
   && MET> 0
   //
   ){     
    H1 = leptons.at(0)+leptons.at(1)+neutrinos.at(0)+neutrinos.at(1);
    // deal with the jets -- separate analysis by number of fattags
    int nfat=0, nbtag=0, nmistag=0;
    for(int i=0;i<jets.size();i++) // count
	{ nbtag = nbtag + btag[i]; nmistag = nmistag + bmistag[i];} 
    if(jets.size() > 2 && fattag.size() >0) {
       H2=jets.at(fattag[0]); 
       // quality requirements   
       double massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
       double rapDiff = abs(H1.eta() - H2.eta());
       Xres = H1 +H2;
cout << " category "<<"1"<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl; 
       if( 
//         massDiff < tolerance && //rapDiff < deltaEtaHH &&
//         (H1.m() > Hmin && H1.m() < Hmax) &&
         (H2.m() > Hmin && H2.m() < Hmax)
         && (H1+H2).m() >MHH
         //&& rapDiff<DetaHH
	 && MET < MnunuMax
	 && (leptons.at(0)+leptons.at(1)).m() < MeeMax
         //&& abs(jets.at(jetn1[minM]).eta()+jets.at(jetn2[minM]).eta())<DetaH
         ){ found=true; Cat->Fill(0.,weight);
           //if(Xres.m()<250) 
           
         } // close quality 
    } // close 1 tag
    if(jets.size() > 3 && !found){
      // find the Hbb candidate
      // pair H2 the jets by the minimum invariant mass difference with H1
      std::vector<double> a1; 
      double invmassB = H1.m(); 
      for(int nj1=0; nj1< jets.size(); nj1++)  if(nj1 != vbftag[0] && nj1 != vbftag[1]) 
         for(int nj2=nj1+1; nj2< jets.size(); nj2++) if(nj2 != vbftag[0] && nj2 != vbftag[1]) { 
	   //std::cout<<nj1<<nj2<<" "<<nj3<<nj4<<std::endl;
	   double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
	   a1.push_back((invmassA-invmassB)*(invmassA-invmassB)); 
	   jetn1.push_back(nj1);jetn2.push_back(nj2); // we also what to keep the nj...           
       } // loop on jets
       //Find the minumum value of the vector (iterator version)
       minM = TMath::LocMin(a1.size(), &a1[0]);
       //std::cout<<"hi, the jets pairs are !!!! "<<jetn1[minM]<<jetn2[minM]<<" "
	  //<<fattag[0]<<std::endl;
       H2=jets.at(jetn1[minM])+jets.at(jetn2[minM]);
       Xres = H2+leptons.at(0)+leptons.at(1)+neutrinos.at(0)+neutrinos.at(1); 
       // quality requirements   
       //double massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
       //double rapDiff = abs(H1.eta() - H2.eta());
       if( 
         //massDiff < tolerance && //rapDiff < deltaEtaHH &&
         //(H1.m() > Hmin && H1.m() < Hmax) &&
         (H2.m() > Hmin && H2.m() < Hmax)
         && (H1+H2).m() >MHH
         && MET < MnunuMax
	 && (leptons.at(0)+leptons.at(1)).m() < MeeMax 
         //&& rapDiff<DetaHH
         && abs(jets.at(jetn1[minM]).eta()+jets.at(jetn2[minM]).eta())<DetaH
         ){ found=true; Cat->Fill(0.,weight);
           //if(Xres.m()<250) 
          //cout << " category "<<"0"<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
       } // close quality 
    } 
    else return found; // close if 2 tags// close if 0 tag
  //////////////////////////////////////
  if(found){
       ///////////////////////////
       // fill the histos
       ///////////////////////////
       //cout<<"fat tag = " <<nfat <<" number of plots "<<basicHiggses.size()<<endl;
       const int numvar1=16;
       if(nbtag >cat){
         double monovar[numvar1] = {
          H1.m(),H1.pt(),H1.eta(),H1.phi(), //4
          H2.m(),H2.pt(),H2.eta(),H2.phi(), // 4
          Xres.m(),Xres.pt(),Xres.eta(),Xres.phi(), // 4
	  nfat,nbtag,nmistag,abs(H1.eta() - H2.eta())
         }; //25
         for (int j=0; j<basicHiggses.size(); j++) basicHiggses[j]->Fill(monovar[j],weight); // weight
         ///////////////////////////
         // fill the jet histos
         ///////////////////////////
         const int numvar2=10;
         PseudoJet vbfmass = jets.at(vbftag[0]) +jets.at(vbftag[1]); 
         double Deta = abs(jets.at(vbftag[0]).eta() -jets.at(vbftag[1]).eta());
         //cout<<"plots vbf "<<basicvbf.size()<<endl;
         double monovarvbf[numvar2] = {
          jets.at(vbftag[0]).pt(),jets.at(vbftag[0]).eta(), //2
          jets.at(vbftag[1]).pt(),jets.at(vbftag[1]).eta(), //2
          vbfmass.m(),vbfmass.pt(),Deta, // 2
	  nfat,nbtag,nmistag
         }; //25
         for (int j=0; j<basicvbf.size(); j++) basicvbf[j]->Fill(monovarvbf[j],weight); // weight
         /////////////////////////////////
         // fill leptons histos
         ////////////////////////////////
         //basicLeptons
         const int numvar3=12;
         double detalep; detalep = abs((leptons.at(0).eta() - leptons.at(1).eta()));
         double drlep; drlep = LepIso[0];//jets.at(jetn1[minM]).delta_R(jets.at(jetn2[minM]));
	 //cout<<" "<<basicLeptons.size()<<" "<<drlep<<" "<<detalep<<" "<<(leptons.at(0)+leptons.at(1)).pt()<<endl;
         double detabb =-1 ,drbb=-1; 
         if (jets.at(jetn2[minM]).m()>0) {
            drbb=jets.at(jetn1[minM]).delta_R(jets.at(jetn2[minM]));//
            detabb=abs(jets.at(jetn1[minM]).eta()-jets.at(jetn2[minM]).eta());
         }
         //cout<<"plots vbf "<<basicvbf.size()<<endl;
         double monovarlep[numvar3] = {
          leptons.at(0).pt(),leptons.at(0).eta(), //2
          leptons.at(1).pt(),leptons.at(1).eta(), //2
          (leptons.at(0)+leptons.at(1)).m(),(leptons.at(0)+leptons.at(1)).pt(),detalep , // 2
	  (H2+leptons.at(0)+leptons.at(1)).m(),
          detabb,drlep,drbb,MET
         }; //25
         for (int j=0; j<basicLeptons.size(); j++) basicLeptons[j]->Fill(monovarlep[j],weight); // weight
       } // close if correct btag
     btagselected->Fill(nbtag,weight); 
  }  // close if fill 
     // histos -- met and Mee Ptee 
  } // close if the leptons pass the cuts
return found; // close if 2 tags
}//close WWbb
/////////////////////////////////////////////////////////////////
bool analyse_4b(
	vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag,vector<int> vbftag, int Xmass){
  // pair the jets that was not VBF tagged
  bool found = false;
  PseudoJet Xres;
  // now we separate analysis
  // number of fattags
  int nfat=0, nbtag=0, nmistag=0;
  double massDiff; 
  for(int i=0;i<jets.size();i++) 
	{ nbtag = nbtag + btag[i]; nmistag = nmistag + bmistag[i];}
  //for(int i=0;i<fattag.size();i++) 
  nfat = fattag.size();// nfat + fattag[i];
  PseudoJet H1,H2;
  float Hmin = higgs_mass*(1-tolerance);
  float Hmax = higgs_mass*(1+tolerance);
  double drbb =-1, drbb2 =-1;
  std::vector<double> a1,a2; 
  std::vector< int > jetn1, jetn2,jetn3, jetn4; // to keep the pairs resolved
  std::vector<double> a3; 
  std::vector< int > jetn11, jetn21; // to keep the pairs 1 tag
  int minM,minM2,cate=-2;
  if(jets.size() > 3 && nfat >1) { // if 2 tag
    //std::cout<<"2 tag! "<<std::endl;
      H1=jets.at(fattag[0]);
      H2=jets.at(fattag[1]);
    // quality requirements   
    massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
    double rapDiff = abs(H1.eta() - H2.eta());
    Xres = H1 +H2; 
    if( //massDiff1 < tolerance 
       //&& massDiff2 < tolerance //rapDiff < deltaEtaHH &&
       H1.m() > Hmin && H1.m() < Hmax
       && H2.m() > Hmin && H2.m() < Hmax
       && (H1+H2).m() >MHH
       && Xres.m() < Xmass*(1.+toleranceX)
       && Xres.m() > Xmass*(1.-toleranceX) 
       //&& rapDiff<DetaHH
      ){  found=true;
	cate=2; 
        //if(Xres.m()<250) cout << " category "<<2<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
      }   
  } // close 2 tag
   else if(nfat >0 && jets.size() > 4 &&  found==false) { // if 1 tag
   //std::cout<<"1 tag! "<<std::endl;
   if(nfat>1) {cout<<"ops! "<<endl;}
   //int nj; for(nj=0; nj< jets.size(); nj++) if(fattag[nj]==1) {
   H1=jets.at(fattag[0]);//}
   // pair H2 the jets by the minimum invariant mass difference with H1
   double invmassB =  H1.m();
   for(int nj1=0; nj1< jets.size(); nj1++)  if(nj1 != vbftag[0] && nj1 != vbftag[1] && nj1!=fattag[0]) 
     for(int nj2=nj1+1; nj2< jets.size(); nj2++) if(nj2 != vbftag[0] && nj2 != vbftag[1] && nj2!=fattag[0]) { 
	   //std::cout<<nj1<<nj2<<" "<<nj3<<nj4<<std::endl;
	   double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
	   a3.push_back((invmassA-invmassB)*(invmassA-invmassB)); 
	   jetn11.push_back(nj1);jetn21.push_back(nj2); // we also what to keep the nj...           
    } // loop on jets
    //int minM;
    //Find the minumum value of the vector (iterator version)
    minM = TMath::LocMin(a3.size(), &a3[0]);
    //std::cout<<"hi, the jets pairs are !!!! "<<jetn1[minM]<<jetn2[minM]<<" "
	//<<fattag[0]<<std::endl;
    H2=jets.at(jetn11[minM])+jets.at(jetn21[minM]);
    // quality requirements   
    massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
    Xres = H1 +H2; 
    double rapDiff = abs(H1.eta() - H2.eta());
    if( //massDiff < tolerance && //rapDiff < deltaEtaHH &&
       H1.m() > Hmin && H1.m() < Hmax &&
       H2.m() > Hmin && H2.m() < Hmax
       && (H1+H2).m() >MHH
       && Xres.m()< Xmass*(1+toleranceX)
       && Xres.m()> Xmass*(1-toleranceX) 
       //&& rapDiff<DetaHH
       //&& abs(jets.at(jetn1[minM]).eta()+jets.at(jetn2[minM]).eta())<DetaH
      ){ found=true; 
        //if(Xres.m()<250)
	cate=1; 
      } //if(jets.size() < 6 && found==false && nbtag==4) cout << " category "<<1<<" mass "<<Xres.m()<<" "<< massDiff<<endl;
  } // close if 1 tag
  else if(jets.size() > 5 && found==false) { // resolved
   // pair the jets by the minimum invariant mass difference
   //std::cout<<"resolved! "<<std::endl;
   // pair by minimum distance of higgs of one pair:test
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
              a2.push_back((invmassA-invmassB)*(invmassA-invmassB)); 
              //cout << " higgs masses "<<invmassA<<" "<<invmassB <<endl;
              //cout << " higgs masses wrong "<<(jets.at(nj1)+jets.at(nj3)).m()<<" "<<(jets.at(nj4)+jets.at(nj2)).m() <<endl;
	      jetn1.push_back(nj1);jetn2.push_back(nj2); // we also what to keep the nj...
	      jetn3.push_back(nj3);jetn4.push_back(nj4);
           }// close if not btagged
    } // loop on jets
    //Find the minumum value of the vector (iterator version)
    minM2 = TMath::LocMin(a2.size(), &a2[0]);
    //std::cout<<"hi, the jets pairs are !!!! "<<jetn1[minM]<<jetn2[minM]<<" "
	//<<jetn3[minM]<<jetn4[minM]<<std::endl;
    H1=jets.at(jetn1[minM2])+jets.at(jetn2[minM2]);
    H2=jets.at(jetn3[minM2])+jets.at(jetn4[minM2]);  
    // quality requirements   
    massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
    double rapDiff = abs(H1.eta() - H2.eta());
    Xres = H1 +H2;
    if( 
       //massDiff < tolerance && //rapDiff < deltaEtaHH &&
       H1.m() > Hmin && H1.m() < Hmax 
       && H2.m() > Hmin && H2.m() < Hmax
       && (H1+H2).m() >MHH
       && Xres.m()< Xmass*(1.+toleranceX)
       && Xres.m()> Xmass*(1.-toleranceX) 
       //&& abs(jets.at(jetn1[minM2]).eta()+jets.at(jetn2[minM2]).eta())<DetaH
       //&& abs(jets.at(jetn3[minM2]).eta()+jets.at(jetn4[minM2]).eta())<DetaH
      ){ //std::cout<<"getting there"<<std::endl; 
	found=true;
	cate=0;  
        //if(Xres.m()<250) cout << " category "<<0<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
      } 
//    
  } //else cout<<"njets"<<jets.size()<<" nfat "<< nfat<<endl; // close if resolved
  //////////////////////////////////////
  //if(massDiff > 0.1 && nbtag>3) cout << "failed , btags = "<<H1.m()<<" "<<H2.m()<<" "<<cate<<endl;
  //cout << "Mass hyp = "<<Xmass<<endl;
  if(found){
       ///////////////////////////
       // fill the histos
       ///////////////////////////
       //cout<<"fat tag = " <<nfat <<" number of plots "<<basicHiggses.size()<<endl;
       const int numvar1=16;
       if(nbtag >cat){
         double monovar[numvar1] = {
          H1.m(),H1.pt(),H1.eta(),H1.phi(), //4
          H2.m(),H2.pt(),H2.eta(),H2.phi(), // 4
          Xres.m(),Xres.pt(),Xres.eta(),Xres.phi(), // 4
	  nfat,nbtag,nmistag,abs(H1.eta() - H2.eta())
         }; //25
         for (int j=0; j<basicHiggses.size(); j++) basicHiggses[j]->Fill(monovar[j],weight); // weight
         ///////////////////////////
         // fill the jet histos
         ///////////////////////////
         const int numvar2=10;
         PseudoJet vbfmass = jets.at(vbftag[0]) +jets.at(vbftag[1]); 
         double Deta = abs(jets.at(vbftag[0]).eta() -jets.at(vbftag[1]).eta());
         //cout<<"plots vbf "<<basicvbf.size()<<endl;
         double monovarvbf[numvar2] = {
          jets.at(vbftag[0]).pt(),jets.at(vbftag[0]).eta(), //2
          jets.at(vbftag[1]).pt(),jets.at(vbftag[1]).eta(), //2
          vbfmass.m(),vbfmass.pt(),Deta, // 2
	  nfat,nbtag,nmistag
         }; //25
         for (int j=0; j<basicvbf.size(); j++) basicvbf[j]->Fill(monovarvbf[j],weight); // weight
        /////////////////////////////////
         // fill leptons histos with 0
         ////////////////////////////////
         //basicLeptons
         const int numvar3=12;
         if(minM2 >0) {
            drbb=jets.at(jetn1[minM2]).delta_R(jets.at(jetn2[minM2]));//
            drbb2=jets.at(jetn3[minM2]).delta_R(jets.at(jetn4[minM2]));
         } else if(minM >0) {drbb=jets.at(jetn11[minM]).delta_R(jets.at(jetn21[minM])); drbb2=0;}
         //cout<<"plots vbf "<<basicvbf.size()<<endl;
         double monovarlep[numvar3] = {
          0,0, //2
          0,0, //2
          0,0,0, // 2
	  0,0,drbb,drbb2,0
         }; //25
         for (int j=0; j<basicLeptons.size(); j++) basicLeptons[j]->Fill(monovarlep[j],weight); // weight
         Cat->Fill(cate,weight);
       } // close if correct btag
     btagselected->Fill(nbtag,weight); 
  }  // close if fill 
  //else cout << " category "<<nfat<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
return found; // close if 2 tags
} // close 4b analysis
/////////////////////////////////////////////////////////////////
bool findVBF(vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag, vector<int> & vbftag){
  // highest invariant mass among the non-tagged
  vector<int> listNonTag; // list the non tagged
  for(int nj1=0; nj1< jets.size(); nj1++) { 
  //cout<<"fat tagged = " <<fattag[nj1]<<" b-quarks mistagged = " <<bmistag[nj1] <<" b-quark = " <<btag[nj1] <<endl;
	if(1>0
           && (1>0
          // && fattag[nj1]==0 
           && btag[nj1] ==0 
           //&& bmistag[nj1]==0
           ) ) {listNonTag.push_back(nj1);}
  }
  if(listNonTag.size()>1){   // find the hightest inv mass pair 
    int nfat=0, nbtag=0, nmistag=0;
    for(int i=0;i<listNonTag.size();i++) 
	{//nfat = nfat + fattag[listNonTag[i]]; 
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
    double MJJ = (jets.at(vbftag[0])+jets.at(vbftag[1])).m();
    if( 1>0 
        && MJJ > Mjj && etaVBF > DeltayVBF 
        && jets.at(vbftag[0]).pt()>jet1_ptminvbf
        && jets.at(vbftag[1]).pt()>jet2_ptminvbf
        && jets.at(vbftag[0]).eta() < etaj
        && jets.at(vbftag[1]).eta() < etaj
      ){ // apply the VBF cuts
       // std::cout<<"hi VBF jets really are !!!! "<<vbftag[0]<<" "<<vbftag[1]<<std::endl;
       //    cout << " dijet mass "<<(jets.at(vbftag[0])+jets.at(vbftag[1])).m() <<endl;
       //  cout << " dijet btag "<<btag[vbftag[0]]<<btag[vbftag[1]] <<endl;
    return true;
  } else return false; // close VBF cuts
 } else  { //cout<<"not enought vbf ! "<<listNonTag.size()<<endl; 
         return false;}// close if listnontag>1 
} // close VBF selection
/////////////////////////////////////////////////////////////////
void isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag);
/////////////////////////////
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets_akt, vector<int> & btag, vector<int> & bmistag, vector<int> & fattag){
  JetDefinition akt(antikt_algorithm, RR);
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
    if(tagged_jet.m()>10) {
      PseudoJet filtered_jet = filter(jets_akt.at(i)); // filter to tag
      if(filtered_jet.m()>Mfat) fattag.push_back(i); //see = 1; else see = 0; // no fat tag 
    } //else see = 0; 
    
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
            && constitu.at(j).pt() > bjetpt
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
int save_hist(int nmass, bool resonant,bool bkg, bool fourb){
  const char* Mass;
  if(resonant && !bkg) {
    if(fourb) Mass = Form("bulk_graviton_mad/Control_shower_%d.root",nmass);
    if(!fourb) Mass = Form("bulk_graviton_mad_WWbb/Control_shower_%d.root",nmass);
    //Mass = Form("histos/Madgraph0_0137/Control_shower_%d.root",nmass);
  }
  else if(!bkg && fourb) Mass = Form("histosnonres/Control_shower_%d.root",nmass);
  else if(!bkg && !fourb) Mass = Form("nonresWWbb/Control_shower_%d.root",nmass);
  else Mass = Form("4bsbkg/Control_shower_%d.root",nmass);
  TFile f1(Mass, "recreate");
  f1.cd();
  Njets_passing_kLooseID->Write();
  Cat->Write();
  btagselected->Write();
  gen_higgs->Write();
  for (int j=0; j<basicLeptons.size(); j++){basicLeptons[j]->Write();} //
  for (int j=0; j<basicHiggses.size(); j++){basicHiggses[j]->Write();} //
  for (int i=0; i<basicvbf.size(); i++) {basicvbf[i]->Write();}
  f1.Close();
  //
  Njets_passing_kLooseID->Reset();
  Cat->Reset();
  gen_higgs->Reset();
  btagselected->Reset();
  for (int i=0; i<basicLeptons.size(); i++) basicLeptons[i]->Reset();
  for (int i=0; i<basicHiggses.size(); i++) basicHiggses[i]->Reset();
  for (int i=0; i<basicvbf.size(); i++) basicvbf[i]->Reset();
  return 1;
}
///////////////////////////////////////////////////////////////////////////
// declare the histos
int decla(int mass){
   const char* label="kk graviton ";

	Njets_passing_kLooseID = new TH1D("njets_passing_kLooseID_ct4",  
		label, 
		13, -0.5, 12.5);
	Njets_passing_kLooseID->GetYaxis()->SetTitle("");
	Njets_passing_kLooseID->GetXaxis()->SetTitle("Njets after showering");

	gen_higgs = new TH1D("gen_higgs",  
		label, 
		90, 50, 1550);
	gen_higgs->GetYaxis()->SetTitle("");
	gen_higgs->GetXaxis()->SetTitle("gen Mhh");

	btagselected = new TH1D("btagselected",  
		label, 
		13, -0.5, 12.5);
	btagselected->GetYaxis()->SetTitle("");
	btagselected->GetXaxis()->SetTitle("b-tagable b's on selected events");
        //

	TH1D *H1hist = new TH1D("H1hist",  
		label, 
		90, -10, 385);
	H1hist->GetYaxis()->SetTitle("Events/ 2 GeV");
	H1hist->GetXaxis()->SetTitle("M_{H1} (GeV)");
	basicHiggses.push_back (H1hist); 

	TH1D *H1histpt = new TH1D("H1histpt",  
		label, 
		20, 0, 300);
	H1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	H1histpt->GetXaxis()->SetTitle("H1 P_T (GeV)");
	basicHiggses.push_back (H1histpt); 

	TH1D *H1histeta = new TH1D("H1histeta",  
		label, 
		30, -6, 6);
	H1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	H1histeta->GetXaxis()->SetTitle("#eta_{H1} (GeV)");
	basicHiggses.push_back (H1histeta); 

	TH1D *H1histphi = new TH1D("H1histphi",  
		label, 
		30, 0, 5);
	H1histphi->GetYaxis()->SetTitle("Events/ 2 GeV");
	H1histphi->GetXaxis()->SetTitle("#phi_{H1} (GeV)");
	basicHiggses.push_back (H1histphi); 

	//

	TH1D *H2hist = new TH1D("H2hist",  
		label, 
		90, -10, 385);
	H2hist->GetYaxis()->SetTitle("Events/ 2 GeV");
	H2hist->GetXaxis()->SetTitle("M_{H2} (GeV)");
	basicHiggses.push_back (H2hist); 

	TH1D *H2histpt = new TH1D("H2histpt",  
		label, 
		20, 0, 300);
	H2histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	H2histpt->GetXaxis()->SetTitle("H2 P_T (GeV)");
	basicHiggses.push_back (H2histpt); 

	TH1D *H2histeta = new TH1D("H2histeta",  
		label, 
		30, -6, 6);
	H2histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	H2histeta->GetXaxis()->SetTitle("#eta_{H2} (GeV)");
	basicHiggses.push_back (H2histeta); 

	TH1D *H2histphi = new TH1D("H2histphi",  
		label, 
		30, 0, 5);
	H2histphi->GetYaxis()->SetTitle("Events/ 2 GeV");
	H2histphi->GetXaxis()->SetTitle("#phi_{H2} (GeV)");
	basicHiggses.push_back (H2histphi); 

	//

	TH1D *RRadMass = new TH1D("RadMass_ct4",  
		label, 
		80, 50, 1100);
	RRadMass->GetXaxis()->SetTitle("M_{ HH } (GeV)");
	basicHiggses.push_back (RRadMass);
	
	TH1D *Radion_pt = new TH1D("radion_pt_ct4",  
		label, 
		50, 0, 300);
	Radion_pt->GetYaxis()->SetTitle("");
	Radion_pt->GetXaxis()->SetTitle("Pt_{HH} (GeV)");
	basicHiggses.push_back (Radion_pt);
	
	TH1D *Radion_eta = new TH1D("radion_eta_ct4",  
		label, 
		35, -9 , 9);
	Radion_eta->GetYaxis()->SetTitle("");
	Radion_eta->GetXaxis()->SetTitle("#eta_{HH }");
	basicHiggses.push_back (Radion_eta);
	//
	TH1D *Radion_phi = new TH1D("radion_phi_ct4",  
		label, 
		25, 0, 7.14);
	Radion_phi->GetYaxis()->SetTitle("");
	Radion_phi->GetXaxis()->SetTitle("#phi_{HH}");
	basicHiggses.push_back (Radion_phi);

	//

	TH1D *Nfattag = new TH1D("fattag_ct4",  
		label, 
		6, -1.5, 4.5);
	Nfattag->GetYaxis()->SetTitle("");
	Nfattag->GetXaxis()->SetTitle("Number of fat tags");
	basicHiggses.push_back (Nfattag);

	TH1D *Nbtag = new TH1D("btag_ct4",  
		label, 
		9, -1.5, 7.5);
	Nbtag->GetYaxis()->SetTitle("");
	Nbtag->GetXaxis()->SetTitle("Number b-jets");
	basicHiggses.push_back (Nbtag);

	TH1D *Nbmistag = new TH1D("bmistag_ct4",  
		label, 
		7, -1.5, 5.5);
	Nbmistag->GetYaxis()->SetTitle("");
	Nbmistag->GetXaxis()->SetTitle("Number mistagged b jets's");
	basicHiggses.push_back (Nbmistag);

	TH1D *DetaHH = new TH1D("detahh_ct4",  
		label, 
		30, -1.5, 10.5);
	DetaHH->GetYaxis()->SetTitle("");
	DetaHH->GetXaxis()->SetTitle("#Delta #eta_{HH}");
	basicHiggses.push_back (DetaHH);

	Cat = new TH1D("cat_ct4",  
		label, 
		7, -1.5, 5.5);
	Cat->GetYaxis()->SetTitle("");
	Cat->GetXaxis()->SetTitle("Tag category");
	//basicHiggses.push_back (Cat);
	/////////////////////////////////////////////////////////////////////

	TH1D *j1histpt = new TH1D("j1histpt",  
		label, 
		30, 20, 300);
	j1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	j1histpt->GetXaxis()->SetTitle("vbf j1 P_T (GeV)");
	basicvbf.push_back (j1histpt); 

	TH1D *j1histeta = new TH1D("j1histeta",  
		label, 
		32, -6, 6);
	j1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	j1histeta->GetXaxis()->SetTitle("vbf #eta_{j1} (GeV)");
	basicvbf.push_back (j1histeta); 

	//

	TH1D *j2histpt = new TH1D("j2histpt",  
		label, 
		30, 20, 100);
	j2histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	j2histpt->GetXaxis()->SetTitle("vbf j2 P_T (GeV)");
	basicvbf.push_back (j2histpt); 

	TH1D *j2histeta = new TH1D("j2histeta",  
		label, 
		32, -6, 6);
	j2histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	j2histeta->GetXaxis()->SetTitle("vbf #eta_{j2} (GeV)");
	basicvbf.push_back (j2histeta); 

	//

	TH1D *RRadMassvbf = new TH1D("RadMass_ct4vbf",  
		label, 
		80, 50, 1200);
	RRadMassvbf->GetXaxis()->SetTitle("vbf M_{j j} (GeV)");
	basicvbf.push_back (RRadMassvbf);

	//

	TH1D *RRadPtvbf = new TH1D("RadPt_ct4vbf",  
		label, 
		60, 50, 250);
	RRadPtvbf->GetXaxis()->SetTitle("vbf pt_{j j} (GeV)");
	basicvbf.push_back (RRadPtvbf);

	//

	TH1D *Detavbf = new TH1D("Deta",  
		label, 
		30, 0, 10);
	Detavbf->GetXaxis()->SetTitle("vbf #Delta #eta");
	basicvbf.push_back (Detavbf);
	
	//

	TH1D *Nfattagvbf = new TH1D("fattag_ct4vbf",  
		label, 
		6, -1.5, 4.5);
	Nfattagvbf->GetYaxis()->SetTitle("");
	Nfattagvbf->GetXaxis()->SetTitle("Number of fat tags vbf");
	basicvbf.push_back (Nfattagvbf);

	TH1D *Nbtagvbf = new TH1D("btag_ct4vbf",  
		label, 
		7, -1.5, 5.5);
	Nbtagvbf->GetYaxis()->SetTitle("");
	Nbtagvbf->GetXaxis()->SetTitle("Number b-jets vbf");
	basicvbf.push_back (Nbtagvbf);

	TH1D *Nbmistagvbf = new TH1D("bmistag_ct4vbf",  
		label, 
		7, -1.5, 5.5);
	Nbmistagvbf->GetYaxis()->SetTitle("");
	Nbmistagvbf->GetXaxis()->SetTitle("Number mistagged b jets's vbf");
	basicvbf.push_back (Nbmistagvbf);

	/////////////////////////////////////////////////////////////////////////////
	// for leptons

	TH1D *E1histpt = new TH1D("E1histpt",  
		label, 
		50, 0, 200);
	E1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	E1histpt->GetXaxis()->SetTitle("lepton 1 P_T (GeV)");
	basicLeptons.push_back (E1histpt); 

	TH1D *E1histeta = new TH1D("E1histeta",  
		label, 
		30, -6, 6);
	E1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	E1histeta->GetXaxis()->SetTitle("#eta_{lepton1} (GeV)");
	basicLeptons.push_back (E1histeta); 

	//

	TH1D *E2histpt = new TH1D("E2histpt",  
		label, 
		50, 0, 200);
	E2histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	E2histpt->GetXaxis()->SetTitle("lepton 2 P_T (GeV)");
	basicLeptons.push_back (E2histpt); 

	TH1D *E2histeta = new TH1D("E2histeta",  
		label, 
		30, -6, 6);
	E2histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	E2histeta->GetXaxis()->SetTitle("#eta_{lepton2} (GeV)");
	basicLeptons.push_back (E2histeta); 

	TH1D *LLMass = new TH1D("LLMass_ct4",  
		label, 
		80, 0, 600);
	LLMass->GetXaxis()->SetTitle("M_{ ll } (GeV)");
	basicLeptons.push_back (LLMass);
	
	TH1D *LL_pt = new TH1D("LL_pt_ct4",  
		label, 
		50, 0, 1400);
	LL_pt->GetYaxis()->SetTitle("");
	LL_pt->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
	basicLeptons.push_back (LL_pt);

	TH1D *LL_Deta = new TH1D("LL_Deta",  
		label, 
		45, -10, 10);
	LL_Deta->GetYaxis()->SetTitle("");
	LL_Deta->GetXaxis()->SetTitle("#Delta#eta_{ll} (GeV)");
	basicLeptons.push_back (LL_Deta);

	TH1D *MetMass = new TH1D("MetMass_ct4",  
		label, 
		50, 0, 750);
	MetMass->GetXaxis()->SetTitle("M_{ bbee } (GeV)");
	basicLeptons.push_back (MetMass);

	TH1D *Detabb = new TH1D("Detabb",  
		label, 
		30, 0, 15);
	Detabb->GetXaxis()->SetTitle("#Delta #eta bb");
	basicLeptons.push_back (Detabb);

	TH1D *DRll = new TH1D("DRll_test",  
		label, 
		30, 0, 15);
	DRll->GetXaxis()->SetTitle("#Delta R ll");
	basicLeptons.push_back (DRll);

	TH1D *DRbb = new TH1D("DRbb",  
		label, 
		30, 0, 15);
	DRbb->GetXaxis()->SetTitle("#Delta R b");
	basicLeptons.push_back (DRbb);

	TH1D *MetMass1 = new TH1D("MetMass1_ct4",  
		label, 
		50, 0, 750);
	MetMass1->GetXaxis()->SetTitle("M_{ #nu#nu } (GeV)");
	basicLeptons.push_back (MetMass1);

return 1;
}


