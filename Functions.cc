#include <fstream>
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
//#include "fastjet/contrib/Nsubjettiness.hh"
//#include "fastjet/contrib/Njettiness.hh"
//#include "fastjet/contrib/NjettinessPlugin.hh"
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
//using namespace fastjet::contrib;
bool findVBFsimple(vector<PseudoJet> jets, vector<int> & vbftag, vector<int> & listNonTag);
/////////////////////////////////////////////////////////////////
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
  int nfat=0, nbtag=0, nmistag=0;
  float Hmin = higgs_mass*(1-toleranceH1);
  float Hmax = higgs_mass*(1+toleranceH1);
  bool found = false;
  PseudoJet Xres;
  int cate;
  // we first find the leptonic higgs -- do not consider more than 2 leptons
  // do a isolation vector
  vector<double> LepIso;
  unsigned int minM =0;
  unsigned int jsize = jets.size(), vbf1 = vbftag[0], vbf2 = vbftag[1];
  std::vector< int > jetn1, jetn2; // to keep the pair
  LepIso.push_back(leptons.at(0).delta_R(leptons.at(1)));
  for(unsigned int i = 0;i<jsize;i++) {
	LepIso.push_back(leptons.at(0).delta_R(jets.at(i)));
	LepIso.push_back(leptons.at(1).delta_R(jets.at(i)));
  } // close iso vector
  int MinDRLep = TMath::LocMin(LepIso.size(), &LepIso[0]);
  //if(MinDRLep<0) for(int i = 0;i<LepIso.size();i++)
  //cout <<LepIso[MinDRLep]<< endl;   
  double MTnunu=(neutrinos.at(0)+neutrinos.at(1)).mperp();
  if(1>0
   && LepIso[MinDRLep] > lepiso  
   && leptons.at(0).pt()> ptlepton 
   && leptons.at(1).pt()> ptlepton
   && leptons.at(0).eta()< etab 
   && leptons.at(1).eta()< etab
   && MTnunu > 0
   //
   ){     
    H1 = leptons.at(0)+leptons.at(1)+neutrinos.at(0)+neutrinos.at(1);
    // deal with the jets -- separate analysis by number of fattags
    for(unsigned int i=0;i<jsize;i++) { nbtag = nbtag + btag[i]; nmistag = nmistag + bmistag[i];} 
    if(jsize > 2 && fattag.size() >0) {
       H2=jets.at(fattag[0]); 
       // quality requirements   
       //double massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
       //double rapDiff = abs(H1.eta() - H2.eta());
       Xres = H1 +H2;
       // cout << " category "<<"1"<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl; 
       if( 
         (H2.m() > Hmin && H2.m() < Hmax)
         && (H1+H2).m() >MHH
	 && MTnunu < MnunuMax
	 && (leptons.at(0)+leptons.at(1)).m() < MeeMax
         //&& abs(jets.at(jetn1[minM]).eta()+jets.at(jetn2[minM]).eta())<DetaH
         ){ found=true; cate =1;         } // close quality 
           //cout << " category "<<cate<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
    } // close 1 tag
    if(jsize > 3 && !found){
      // pair H2 the jets by the minimum invariant mass difference with H1
      std::vector<double> a1; 
      double invmassB = H1.m(); 
      for(unsigned int nj1=0; nj1< jsize; nj1++)  if(nj1 != vbf1 && nj1 != vbf2) 
         for(unsigned int nj2=nj1+1; nj2< jsize; nj2++) if(nj2 != vbf1 && nj2 != vbf2) { 
	   double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
	   a1.push_back((invmassA-invmassB)*(invmassA-invmassB)); 
	   jetn1.push_back(nj1);jetn2.push_back(nj2); // we also what to keep the nj...           
       } // loop on jets
       //Find the minumum value of the vector (iterator version)
       minM = TMath::LocMin(a1.size(), &a1[0]);
       H2=jets.at(jetn1[minM])+jets.at(jetn2[minM]);
       Xres = H2+leptons.at(0)+leptons.at(1)+neutrinos.at(0)+neutrinos.at(1); 
       if( 
         (H2.m() > Hmin && H2.m() < Hmax)
         && (H1+H2).m() >MHH
         && MTnunu < MnunuMax
	 && (leptons.at(0)+leptons.at(1)).m() < MeeMax 
         //&& rapDiff<DetaHH
         && abs(jets.at(jetn1[minM]).eta()+jets.at(jetn2[minM]).eta())<DetaH
         ){ found=true; cate=0;//Cat->Fill(0.,weight);
           //if(Xres.m()<250) 
          //cout << " category "<<"0"<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
          } // close quality 
     } // close if 2 tags
    } else return found; // not found 2 leptons
  //////////////////////////////////////
  if(found==true){
       ///////////////////////////
       // fill the histos
       ///////////////////////////
       //cout<<"fat tag = " <<nfat <<" number of plots "<<basicHiggses.size()<<endl;
       const int numvar1=16;
       unsigned int numvar11=16;
      if( 1>0 
         && ((nbtag >cat && cate==0 ) ||
         ( nbtag >cat-1 && cate==1 )) 
         ){
         //cout << " category "<<"0"<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
         Cat->Fill(cate,weight);
         double monovar[numvar1] = {
          H1.m(),H1.pt(),H1.eta(),H1.phi(), //4
          H2.m(),H2.pt(),H2.eta(),H2.phi(), // 4
          Xres.m(),Xres.pt(),Xres.eta(),Xres.phi(), // 4
	  nfat,nbtag,nmistag,abs(H1.eta() - H2.eta())
         }; //25
         for (unsigned int j=0; j<numvar11; j++) basicHiggses[j]->Fill(monovar[j],weight); // weight
         ///////////////////////////
         // fill the jet histos
         ///////////////////////////
         const int numvar2=10;
         unsigned int numvar22=10;
         PseudoJet vbfmass = jets.at(vbf1) +jets.at(vbf2); 
         double Deta = abs(jets.at(vbf1).eta() -jets.at(vbf2).eta());
         //cout<<"plots vbf "<<basicvbf.size()<<endl;
         double monovarvbf[numvar2] = {
          jets.at(vbf1).pt(),jets.at(vbf1).eta(), //2
          jets.at(vbf2).pt(),jets.at(vbf2).eta(), //2
          vbfmass.m(),vbfmass.pt(),Deta, // 2
	  nfat,nbtag,nmistag
         }; //25
         for (unsigned int j=0; j<numvar22; j++) basicvbf[j]->Fill(monovarvbf[j],weight); // weight
         /////////////////////////////////
         // fill leptons histos
         ////////////////////////////////
         //basicLeptons
         const int numvar3=12;
         unsigned int numvar33=12;
         double detalep; detalep = abs((leptons.at(0).eta() - leptons.at(1).eta()));
         double drlep=-1; drlep = LepIso[0];
         double detabb =-1 ,drbb=-1; 
         if(cate ==0)if (jets.at(jetn2[minM]).m()>0 ) {
            drbb=jets.at(jetn1[minM]).delta_R(jets.at(jetn2[minM]));//
            detabb=abs(jets.at(jetn1[minM]).eta()-jets.at(jetn2[minM]).eta());
         }
         cout<<"plots vbf "<<basicLeptons.size()<<" "<<MTnunu<<endl;
         //double met = MET;
         double monovarlep[numvar3] = {
          leptons.at(0).pt(),leptons.at(0).eta(), //2
          leptons.at(1).pt(),leptons.at(1).eta(), //2
          (leptons.at(0)+leptons.at(1)).m(),(leptons.at(0)+leptons.at(1)).pt(),detalep , // 2
	  (H2+leptons.at(0)+leptons.at(1)).m(), detabb,drlep,MTnunu
         }; //25
         for (unsigned int j=0; j<numvar33; j++) basicLeptons[j]->Fill(monovarlep[j],weight); // weight
       } // close if correct btag
     btagselected->Fill(nbtag,weight); 
  }  // close if fill 
     // histos -- met and Mee Ptee 
  //} // close if the leptons pass the cuts
return found; // close if 2 tags
}//close WWbb
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool analyse_4b(
	vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag,
        vector<int> vbftag, int Xmass){
  // pair the jets that was not VBF tagged
  bool found = false; 
  PseudoJet Xres;
  // now we separate analysis
  // number of fattags
  int nfat=0, nbtag=0, nmistag=0;
  double massDiff; 
  int minM,minM2,cate=-2;
  // do a vector with btagged jets and fat btagged jets
  vector<int> bjets, fatbjets, misbjets, misfatbjets, listNonTag; // list non fat b jets
  vector<double> bweights; double final_weight;
  unsigned int jsize = jets.size(), fsize = fattag.size(), vbf1 = vbftag[0], vbf2 = vbftag[1];
  for(unsigned int i=0;i<jsize;i++) { 
    nbtag = nbtag + btag[i]; nmistag = nmistag + bmistag[i];
    if(btag[i]>0) { // think on mistag later
      int teste=0;
      for(unsigned int j=0;j<fsize;j++) {
       unsigned int ff = fattag[j]; if(ff==i) {
         fatbjets.push_back(i);teste=1;
         if(btag[i]==2) bweights.push_back(subjet2b); // calculate distances to decide
         if(btag[i]==1) bweights.push_back(fatjet2b); 
      }}  
      if(teste==0) {bjets.push_back(i); bweights.push_back(normalb);}
    } else {listNonTag.push_back(i); bweights.push_back(1);}
  } 
  // change the analysis to take only b--tagged 
  nfat = fattag.size();// nfat + fattag[i];
  PseudoJet H1,H2;
  float Hmin1 = higgs_mass*(1-toleranceH2);
  float Hmax1 = higgs_mass*(1+toleranceH1);
  float Hmin2 = higgs_mass*(1-toleranceH2);
  float Hmax2 = higgs_mass*(1+toleranceH2);
  double drbb =-1, drbb2 =-1;
  std::vector<double> a1,a2; 
  std::vector< int > jetn1, jetn2,jetn3, jetn4; // to keep the pairs resolved
  std::vector<double> a3; 
  std::vector< int > jetn11, jetn21; // to keep the pairs 1 tag
  if(jsize > 3 && nfat >1) { // if 2 tag
    //std::cout<<"2 tag! "<<std::endl;
      H1=jets.at(fattag[0]);
      H2=jets.at(fattag[1]);
    // quality requirements   
    massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
    double rapDiff = abs(H1.eta() - H2.eta());
    Xres = H1 +H2; 
    if(1>0
       && H1.m() > minMH && H2.m() > minMH 
       && H1.m() > Hmin1 && H1.m() < Hmax1 
       && H2.m() > Hmin2 && H2.m() < Hmax2
       && H1.pt() > H1_ptmin && H2.pt() > H2_ptmin  
       && Xres.pt() > HH_ptmin 
       && (H1+H2).m() >MHH
       && Xres.m() < Xmass*(1.+toleranceX)
       && Xres.m() > Xmass*(1.-toleranceX) 
       && rapDiff<DetaHH
      ){  found=true;
	cate=2; 
        final_weight = bweights[fattag[0]]*bweights[fattag[0]];
        //cout << " category "<<cate<<" mass "<<Xres.m()<<endl;
        //if(Xres.m()<250) cout << " category "<<2<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
      }   
  } // close 2 tag
   else if(nfat >0 && jsize > 4 &&  found==false) { // if 1 tag
   //std::cout<<"1 tag! "<<std::endl;
   if(nfat>1) {cout<<"ops! "<<endl;}
   //int nj; for(nj=0; nj< jets.size(); nj++) if(fattag[nj]==1) {
   H1=jets.at(fattag[0]);//}
   // pair H2 the jets by the minimum invariant mass difference with H1
   double invmassB = H1.m();
   unsigned int fat1 = fattag[0];
   for(unsigned int nj1=0; nj1< jsize; nj1++)  if(nj1 != vbf1 && nj1 != vbf2 && nj1!=fat1) 
     for(unsigned int nj2=nj1+1; nj2< jsize; nj2++) if(nj2 != vbf1 && nj2 != vbf2 && nj2!=fat1) { 
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
       1>0
       && H1.m() > minMH && H2.m() > minMH 
       && H1.m() > Hmin1 && H1.m() < Hmax1 
       && H2.m() > Hmin2 && H2.m() < Hmax2
       && H1.pt() > H1_ptmin && H2.pt() > H2_ptmin  
       && (H1+H2).m() >MHH
       && Xres.pt() > HH_ptmin 
       && Xres.m()< Xmass*(1+toleranceX)
       && Xres.m()> Xmass*(1-toleranceX) 
       && rapDiff<DetaHH
       //&& abs(jets.at(jetn1[minM]).eta()+jets.at(jetn2[minM]).eta())<DetaH
      ){ found=true; 
        //if(Xres.m()<250)
	cate=1; 
        final_weight = bweights[fattag[0]]*bweights[jetn11[minM]]*bweights[jetn21[minM]];
      } // close quality
  } // close if 1 tag
  else if(jsize > 5 && found==false) { // resolved
   // pair the jets by the minimum invariant mass difference
   //std::cout<<"resolved! "<<std::endl;
   for(unsigned int nj1=0; nj1< jsize; nj1++)   
     for(unsigned int nj2=0; nj2< jsize; nj2++) 
       for(unsigned int nj3=0; nj3< jsize; nj3++)     
	 for(unsigned int nj4=0; nj4< jsize; nj4++)   
           { 
           if( 1>0
           && nj1 != vbf1 && nj1 != vbf2
	   && nj2 != vbf1 && nj2 != vbf2 
           && nj3 != vbf1 && nj3 != vbf2
           && nj4 != vbf1 && nj4 != vbf2 
           && nj1 !=nj2 && nj1 !=nj3 && nj1 !=nj4 
           && nj2 !=nj3 && nj2 !=nj4 
           && nj1 !=nj4  
              ){
	      double invmassA =  (jets.at(nj1)+jets.at(nj2)).m();
	      double invmassB =  (jets.at(nj3)+jets.at(nj4)).m();
              if(invmassA< 50 || invmassB<50) a2.push_back(100000000); else
                 a2.push_back((invmassA-invmassB)*(invmassA-invmassB)); 
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
       1>0
       && H1.m() > minMH && H2.m() > minMH 
       && H1.m() > Hmin1 && H1.m() < Hmax1 
       && H2.m() > Hmin2 && H2.m() < Hmax2
       && H1.pt() > H1_ptmin && H2.pt() > H2_ptmin         
       && (H1+H2).m() >MHH
       && Xres.pt() > HH_ptmin 
       && Xres.m()< Xmass*(1.+toleranceX)
       && Xres.m()> Xmass*(1.-toleranceX) 
       && rapDiff<DetaHH
       //&&        H2.m() > 15
      ){ //std::cout<<"getting there"<<std::endl; 
	found=true;
	cate=0;  
        final_weight = bweights[jetn1[minM]]*bweights[jetn2[minM]]*bweights[jetn3[minM]]*bweights[jetn4[minM]];
      } // close quality
  } 
  //////////////////////////////////////
  //if(massDiff > 0.1 && nbtag>3) cout << "failed , btags = "<<H1.m()<<" "<<H2.m()<<" "<<cate<<endl;
  //cout << "Mass hyp = "<<Xmass<<endl;
  if(found){ 
       ///////////////////////////
       // fill the histos
       ///////////////////////////
       //cout<<"fat tag = " <<nfat <<" number of plots "<<basicHiggses.size()<<endl;
       const int numvar1=16;
       unsigned int numvar11=16;
       if( 
         (nbtag >cat && cate==0 ) //|| // && (H1.m()< 50 || H2.m()<50)
         //(nbtag >cat-1 && cate==1) ||
         //(nbtag >cat-2 && cate==2)
         ){
        //if(H1.m()>130) cout << " category "<<cate<<" mass "<<Xres.m()<<" nbtag "<<nbtag<<endl;
         Cat->Fill(cate,weight);
         // cout << " category "<<cate<<" mass "<<Xres.m()<<endl;
         double monovar[numvar1] = {
          H1.m(),H1.pt(),H1.eta(),H1.phi(), //4
          H2.m(),H2.pt(),H2.eta(),H2.phi(), // 4
          Xres.m(),Xres.pt(),Xres.eta(),Xres.phi(), // 4
	  nfat,nbtag,nmistag,abs(H1.eta() - H2.eta())
         }; //25
         for (unsigned int j=0; j<numvar11; j++) basicHiggses[j]->Fill(monovar[j],final_weight); // weight
         ///////////////////////////
         // fill the jet histos
         ///////////////////////////
         const int numvar2=10; 
         unsigned int numvar22=10; 
         PseudoJet vbfmass = jets.at(vbf1) +jets.at(vbf2); 
         double Deta = abs(jets.at(vbf1).eta() -jets.at(vbf2).eta());
         //cout<<"plots vbf "<<basicvbf.size()<<endl;
         double monovarvbf[numvar2] = {
          jets.at(vbf1).pt(),jets.at(vbf1).eta(), //2
          jets.at(vbf2).pt(),jets.at(vbf2).eta(), //2
          vbfmass.m(),vbfmass.pt(),Deta, // 2
	  nfat,nbtag,nmistag
         }; //25
         for (unsigned int j=0; j<numvar22; j++) basicvbf[j]->Fill(monovarvbf[j],final_weight); // weight
        /////////////////////////////////
         // fill leptons histos with 0
         ////////////////////////////////
         //basicLeptons
         const int numvar3=12;
         unsigned int numvar33=12;
          if(minM2 >0) { 
            drbb=jets.at(jetn1[minM2]).delta_R(jets.at(jetn2[minM2]));//
            drbb2=jets.at(jetn3[minM2]).delta_R(jets.at(jetn4[minM2]));
         } else if(minM >0) {drbb=jets.at(jetn11[minM]).delta_R(jets.at(jetn21[minM])); drbb2=0;} 
         double monovarlep[numvar3] = {0,0,0,0,0,0,0,0,0,drbb,drbb2,0}; //25
         for (unsigned int j=0; j<numvar33; j++) basicLeptons[j]->Fill(monovarlep[j],weight); // weight
       } // close if correct btag
     btagselected->Fill(nbtag,weight); 
  }  // close if fill 
  //else cout << " category "<<nfat<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
return found; // close if 2 tags
} // close 4b analysis
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool analyse_4b_prior(
	vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag,
        vector<int> vbftag, int Xmass){
  // pair the jets using btag
  bool found = false, foundvbf =false; 
  PseudoJet Xres;
  // now we separate analysis
  // number of fattags
  int nfat=0, nbtag=0, nmistag=0;
  double massDiff; 
  int minM,minM2,cate=-2;
  // do a vector with btagged jets and fat btagged jets
  vector<int> bjets; // list non fat b jets
  vector<int> fatbjets;
  vector<int> listNonTag; // list the non tagged
  unsigned int jsize = jets.size(), fsize = fattag.size();
  for(unsigned int i=0;i<jsize;i++) { 
    nbtag = nbtag + btag[i]; nmistag = nmistag + bmistag[i];
    if(1>0 
       && btag[i]>0 || (btag[i]==0 && bmistag[i] >0)
    ) { // think on weight later
      int teste=0;
      for(unsigned int j=0;j<fsize;j++) {
        unsigned int ff = fattag[j];
        if(ff==i) {fatbjets.push_back(i);teste=1;} 
      }
      if(teste==0) bjets.push_back(i);
    } else listNonTag.push_back(i);
  } unsigned int lsize = listNonTag.size(), fbsize = fatbjets.size(), bsize = bjets.size();
  //nbfat = fatbjets.size();// nfat + fattag[i];
  PseudoJet H1,H2;
  float Hmin1 = higgs_mass*(1-toleranceH1);
  float Hmax1 = higgs_mass*(1+toleranceH1);
  float Hmin2 = higgs_mass*(1-toleranceH2);
  float Hmax2 = higgs_mass*(1+toleranceH2);
  double drbb =-1, drbb2 =-1;
  std::vector<double> a1,a2; 
  std::vector< int > jetn1, jetn2,jetn3, jetn4; // to keep the pairs resolved
  std::vector<double> a3; 
  std::vector< int > jetn11, jetn21; // to keep the pairs 1 tag
  //////////////////////////////////////////////////////////////
  if(fbsize > 1) { // if 2 tag
    //std::cout<<"2 tag! "<<std::endl;
    H1=jets.at(fatbjets[0]);
    H2=jets.at(fatbjets[0]);
    // quality requirements   
    massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
    double rapDiff = abs(H1.eta() - H2.eta());
    Xres = H1 +H2;    
    if( //massDiff1 < tolerance 
       //&& massDiff2 < tolerance //rapDiff < deltaEtaHH &&
       1>0
       && H1.m() > minMH && H2.m() > minMH 
       && H1.m() > Hmin1 && H1.m() < Hmax1 
       && H2.m() > Hmin2 && H2.m() < Hmax2
       && H1.pt() > H1_ptmin && H2.pt() > H2_ptmin  
       && Xres.pt() > HH_ptmin 
       && (H1+H2).m() >MHH
       && Xres.m() < Xmass*(1.+toleranceX)
       && Xres.m() > Xmass*(1.-toleranceX) 
       && rapDiff<DetaHH
      ){ 
        found=true; cate=2;     
        // then fill the WBF from the non tagged 
        if(lsize>1) foundvbf = findVBFsimple(jets, vbftag, listNonTag); 
       }
        //cout << " category "<<cate<<" mass "<<Xres.m()<<endl;
        //if(Xres.m()<250) cout << " category "<<2<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
  } // close 2 tag
  //////////////////////////////////////////////////////////////////////////
  if(found==false && fbsize > 0 && bsize>1) { // if 1 tag
   //std::cout<<"1 tag! "<<std::endl; 
   // if(nfat>1) {cout<<"ops! "<<endl;}
   //int nj; for(nj=0; nj< jets.size(); nj++) if(fattag[nj]==1) {
   H1=jets.at(fatbjets[0]); 
   // pair H2 the jets by the minimum invariant mass difference with H1
   double invmassB = 125;// H1.m(); 
   if(bsize==2) {H2=jets.at(bjets[0])+jets.at(bjets[1]);}
   else if (bsize>2) {
    for(unsigned int nj1=0; nj1< bsize; nj1++)  if(bjets[nj1]!=fatbjets[0]) 
     for(unsigned int nj2=nj1+1; nj2< bsize; nj2++) if( bjets[nj2]!=fatbjets[0] && bjets[nj2]!=bjets[nj1]) { 
	   //std::cout<<nj1<<nj2<<" "<<nj3<<nj4<<std::endl; 
	   double invmassA =  (jets.at(bjets[nj1])+jets.at(bjets[nj2])).m();
	   a3.push_back((invmassA-invmassB)*(invmassA-invmassB)); 
	   jetn11.push_back(bjets[nj1]);jetn21.push_back(bjets[nj2]); 
           // we also what to keep the nj...           
     } // loop on jets
      //int minM;
      //Find the minumum value of the vector (iterator version)
      minM = TMath::LocMin(a3.size(), &a3[0]);
      //std::cout<<"hi, the jets pairs are !!!! "<<jetn1[minM]<<jetn2[minM]<<" "
	//<<fattag[0]<<std::endl;
      H2=jets.at(jetn11[minM])+jets.at(jetn21[minM]);
    // quality requirements   
    // massDiff = abs(2*(H1.m() - H2.m())/(H1.m() + H2.m()));
    } // close if >2 btag
    Xres = H1 +H2; 
    double rapDiff = abs(H1.eta() - H2.eta());
    if( //massDiff < tolerance && //rapDiff < deltaEtaHH &&
       1>0
       && H1.m() > minMH && H2.m() > minMH 
       && H1.m() > Hmin1 && H1.m() < Hmax1 
       && H2.m() > Hmin2 && H2.m() < Hmax2
       && H1.pt() > H1_ptmin && H2.pt() > H2_ptmin  
       //&& (H1+H2).m() >MHH
       && Xres.pt() > HH_ptmin 
       && Xres.m()< Xmass*(1+toleranceX)
       && Xres.m()> Xmass*(1-toleranceX) 
       && rapDiff<DetaHH
//       && abs(jets.at(jetn1[minM]).eta()+jets.at(jetn2[minM]).eta())<DetaH
      ){ found=true; 
        //if(Xres.m()<250)
	cate=1; 
//if(H2.m() > 110 && H2.m() < 130) cout<<"ops! "<<bjets.size()<<" "<<jets.size()<<" fatjet "<<fattag.size()<<" fatbjet "<<fatbjets.size() <<endl;
        if(listNonTag.size()>1) foundvbf = findVBFsimple(jets, vbftag, listNonTag); 
      } // close quality
      // then fill the WBF from the non tagged 
      //vector<int> listNonTag; // list the non tagged
      //for(int nj1=0; nj1< jets.size(); nj1++)
      //  if(nj1 != fattag[0] || nj1 != jetn11[minM] || nj1 != jetn21[minM]) listNonTag.push_back(nj1); 
      // if not fat tagged keep
  } // close if 1 tag
  //////////////////////////////////////////////////////////////////////
  if(jsize > 5 && found==false) { // resolved
   // pair the jets by the minimum invariant mass difference
   //std::cout<<"resolved! "<<std::endl; //cout<<"ops! "<<bjets.size()<<endl; 
  if (bsize>3) { 
    for(unsigned int nj1=0; nj1< bsize; nj1++)   
     for(unsigned int nj2=0; nj2< bsize; nj2++) 
       for(unsigned int nj3=0; nj3< bsize; nj3++)     
	 for(unsigned int nj4=0; nj4< bsize; nj4++) {
           if( 1>0
           && nj1 !=nj2 && nj1 !=nj3 && nj1 !=nj4 
           && nj2 !=nj3 && nj2 !=nj4 
           && nj1 !=nj4  
              ){
	      double invmassA =  (jets.at(bjets[nj1])+jets.at(bjets[nj2])).m();
	      double invmassB =  (jets.at(bjets[nj3])+jets.at(bjets[nj4])).m();
              a2.push_back((invmassA-invmassB)*(invmassA-invmassB)); 
	      jetn1.push_back(bjets[nj1]);jetn2.push_back(bjets[nj2]); // we also what to keep the nj...
	      jetn3.push_back(bjets[nj3]);jetn4.push_back(bjets[nj4]);
           } }// loop on jets
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
       1>0
       && H1.m() > minMH && H2.m() > minMH 
       && H1.m() > Hmin1 && H1.m() < Hmax1 
       && H2.m() > Hmin2 && H2.m() < Hmax2
       && H1.pt() > H1_ptmin && H2.pt() > H2_ptmin         
       && (H1+H2).m() >MHH
       && Xres.pt() > HH_ptmin 
       && Xres.m()< Xmass*(1.+toleranceX)
       && Xres.m()> Xmass*(1.-toleranceX) 
       && rapDiff<DetaHH
       //&&        H2.m() > 15
      ){ //std::cout<<"getting there"<<std::endl; 
	found=true;
	cate=0;  
        if(listNonTag.size()>1) foundvbf = findVBFsimple(jets, vbftag, listNonTag); 
      } // close quality
      // then fill the WBF from the non tagged 
      //vector<int> listNonTag; // list the non tagged
      //for(int nj1=0; nj1< jets.size(); nj1++)
      //  if(nj1 != jetn1[minM] || nj1 != jetn2[minM] || nj1 != jetn3[minM] || nj1 != jetn4[minM]) 
      //    listNonTag.push_back(nj1); // if not fat tagged keep
    } //if not 4 btag 
  } 
  //////////////////////////////////////
  //if(massDiff > 0.1 && nbtag>3) cout << "failed , btags = "<<H1.m()<<" "<<H2.m()<<" "<<cate<<endl;
  //cout << "Mass hyp = "<<Xmass<<endl;
  if(found && foundvbf){ 
       ///////////////////////////
       // fill the histos
       ///////////////////////////
       //cout<<"fat tag = " <<nfat <<" number of plots "<<basicHiggses.size()<<endl;
       const int numvar1=16;
       if( 
         (cate==0 ) ||
         (cate==1 ) ||
         (cate==2 )
         ){
        //if(H1.m()>130) cout << " category "<<cate<<" mass "<<Xres.m()<<" nbtag "<<nbtag<<endl;
         Cat->Fill(cate,weight);
         //cout << " category "<<cate<<"higgses mass "<<H1.m()<<" "<<H2.m()<<endl;
         double monovar[numvar1] = {
          H1.m(),H1.pt(),H1.eta(),H1.phi(), //4
          H2.m(),H2.pt(),H2.eta(),H2.phi(), // 4
          Xres.m(),Xres.pt(),Xres.eta(),Xres.phi(), // 4
	  nfat,nbtag,nmistag,abs(H1.eta() - H2.eta())
         }; //25
         unsigned int bh = basicHiggses.size(), bl = basicLeptons.size(), bv = 10; // basicvbf.size();
         for (unsigned int j=0; j<bh; j++) basicHiggses[j]->Fill(monovar[j],weight); // weight
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
         for (unsigned int j=0; j<bv; j++) basicvbf[j]->Fill(monovarvbf[j],weight); // weight
        /////////////////////////////////
         // fill leptons histos with 0
         ////////////////////////////////
         //basicLeptons
         const int numvar3=12;
          if(minM2 >0) { 
            drbb=jets.at(jetn1[minM2]).delta_R(jets.at(jetn2[minM2]));//
            drbb2=jets.at(jetn3[minM2]).delta_R(jets.at(jetn4[minM2]));
         } else if(minM >0) {drbb=jets.at(jetn11[minM]).delta_R(jets.at(jetn21[minM])); drbb2=0;} 
         double monovarlep[numvar3] = {0,0,0,0,0,0,0,0,0,drbb,drbb2,0}; //25
         for (unsigned int j=0; j<bl; j++) basicLeptons[j]->Fill(monovarlep[j],weight); // weight
       } // close if correct btag
     btagselected->Fill(nbtag,weight); 
  }  // close if fill 
  //else cout << " category "<<nfat<<" mass "<<Xres.m()<<" btags "<<nbtag <<endl;
return found; // close if 2 tags
} // close 4b analysis prior higgs reco
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
bool findVBFsimple(vector<PseudoJet> jets, vector<int> & vbftag, vector<int> & listNonTag){
    // find the hightest inv mass pair among jets 
    std::vector<double> a1; 
    std::vector< int > jetn1, jetn2; // to keep the pairs
    for(unsigned int nj1=0; nj1< listNonTag.size(); nj1++) 
	for(unsigned int nj2=nj1+1; nj2< listNonTag.size(); nj2++) { // we also what to keep the nj...
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
        && jets.at(vbftag[0]).delta_R(jets.at(vbftag[1])) > DeltaRVBF
        && (jets.at(vbftag[0])+jets.at(vbftag[1])).pt() > PTjj
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
} // close VBF selection
/////////////////////////////////////////////////////////////////
bool findVBF(vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag, vector<int> & vbftag){
  // highest invariant mass among the non-tagged
  vector<int> listNonTag; // list the non tagged
  unsigned int jsize = jets.size();
  for(unsigned int nj1=0; nj1< jsize; nj1++) { 
  //cout<<"fat tagged = " <<fattag[nj1]<<" b-quarks mistagged = 
  //" <<bmistag[nj1] <<" b-quark = " <<btag[nj1] <<endl;
	if(1>0
          // && fattag[nj1]==0 
          // && btag[nj1] ==0 // if not b tagged see
           //&& bmistag[nj1]==0 // if not b mis tagged see
           ) {
               int found=0;
               //for(int nj2=0; nj2< fattag.size(); nj2++)  
               //   {if(fattag[nj2]==nj1){found=1; break;}}
               if(found==0) listNonTag.push_back(nj1); // if not fat tagged keep
             }
  } 
  unsigned int nsize = listNonTag.size();  // cout<<"enter vbf! "<<nsize<<endl;
  if(nsize>1){   // find the hightest inv mass pair 
    int nfat=0, nbtag=0, nmistag=0;
    for(unsigned int i=0;i<nsize;i++) 
	{//nfat = nfat + fattag[listNonTag[i]]; 
         nbtag = nbtag + btag[listNonTag[i]]; nmistag = nmistag + bmistag[listNonTag[i]];}
    std::vector<double> a1; 
    std::vector< int > jetn1, jetn2; // to keep the pairs
    for(unsigned int nj1=0; nj1< nsize; nj1++) 
	for(unsigned int nj2=nj1+1; nj2< nsize; nj2++) { // we also what to keep the nj...
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
        && jets.at(vbftag[0]).delta_R(jets.at(vbftag[1])) > DeltaRVBF
        && (jets.at(vbftag[0])+jets.at(vbftag[1])).pt() > PTjj
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
    double const ptmin=0.0; // if parton jet pt min zero
    Selector jet_selector_parton = SelectorPtMin(ptmin);
    jets_akt = sorted_by_pt(jet_selector_parton(cs_akt.inclusive_jets()));
  }
  unsigned int njets = jets_akt.size();
  //cout<<njets<<endl;
  Njets_passing_kLooseID->Fill(njets,weight);
  isbtagged(jets_akt, btag, bmistag); // check wheather the b(c)jet is b--(mis)tagable
  ///////////////////// check tag
  JetDefinition CA10(cambridge_algorithm, Rsb);
  // Filter definition to improve mass resolution
  Filter filter(JetDefinition(cambridge_algorithm, Rfilt), SelectorNHardest(n_subjet));
  PseudoJet tagged_jet;
  for (unsigned int i = 0; i < njets; i++) { // to each akt jet
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
  //jets = jets_akt;
  //////////////////////////// see n-subjetiness
  /* for (int i = 0; i < jets_akt.size(); i++) { // to each akt jet
    // want the 31 and 21, by now the basic usage of 21
    double beta = 1.0; // beta
    UnnormalizedMeasure measureSpec1(beta);
    OnePass_WTA_KT_Axes     axisMode1;
    //Njettiness::AxesMode axisMode;
    //axisMode = Njettiness::onepass_kt_axes;
    NsubjettinessRatio nSubRatio21(2, 1, axisMode1, measureSpec1);
    double tau21 = nSubRatio21(jets_akt[i]);
  } // to each kt jet */
  return njets;
} // close cluster jets
////////////////////////////////////////////////////////////////////////////////////////////////
void isbtagged(vector<PseudoJet> jets, vector<int> & btag, vector<int> & bmistag){ 
  unsigned int jsize = jets.size();
  for (unsigned int i=0; i<jsize; i++) { // check wheter jet have inside a b's are taggable 
  int see=0,see2=0;
     vector<PseudoJet> constitu=jets.at(i).constituents();
     unsigned int csize = constitu.size();
     for (unsigned int j=0; j<csize; j++) {
       //cout<<"constituents flavour "<<constitu.at(j).user_index()<<endl;
       if((constitu.at(j).user_index() == 5 || constitu.at(j).user_index() == -5) // work !!
            && constitu.at(j).pt() > bjetpt
	    && constitu.at(j).eta() < etab
          ) {see++;}//btag.push_back(1);} else btag.push_back(0); 
       if( abs(constitu.at(j).user_index()) == 4  // work !! 
            && constitu.at(j).pt() > bjetpt
	    && constitu.at(j).eta() < etab
         ) {see2++;}// bmistag.push_back(1);} bmistag.push_back(0);
     } // close constituents
     //bmistag.push_back(see2);
     btag.push_back(see); //else btag.push_back(0); // count all tag/jet
//     if(see>0) btag.push_back(1); else btag.push_back(0); // count only one tag/jet
     if(see2>0) bmistag.push_back(1); else bmistag.push_back(0); // count only one tag/jet
     //cout<<"b-quarks mistagged = " <<bmistag[i] <<" b-quark = " <<btag[i] <<endl;
  } // close for each jet
} // close isbtagged
/////////////////////////////////////////////////////////////////////////
// save the histos
int save_hist(int nmass, bool resonant,bool bkg, bool fourb){
  const char* Mass;
  if(resonant && !bkg) {
    if(fourb) Mass = Form("spin0/Control_shower_%d.root",nmass);
    if(!fourb) Mass = Form("bulk_graviton_mad_WWbb/Control_shower_%d.root",nmass);
    //Mass = Form("histos/Madgraph0_0137/Control_shower_%d.root",nmass);
  }
  else if(!bkg && fourb) Mass = Form("nonresonant/Control_shower_%d.root",nmass);
  else if(!bkg && !fourb) Mass = Form("nonresWWbb/Control_shower_%d.root",nmass);
  else Mass = Form("4bsbkg/Control_shower_%d.root",nmass);
  TFile f1(Mass, "recreate");
  f1.cd();
  Njets_passing_kLooseID->Write();
  Cat->Write();
  btagselected->Write();
  gen_higgs->Write();
  unsigned int bl = basicLeptons.size(), bh = basicHiggses.size(), bv = 10;// basicvbf.size();
  for (unsigned int j=0; j<bl; j++){basicLeptons[j]->Write();} //
  for (unsigned int j=0; j<bh; j++){basicHiggses[j]->Write();} //
  for (unsigned int i=0; i<bv; i++) {basicvbf[i]->Write();}
  f1.Close();
  //
  Njets_passing_kLooseID->Reset();
  Cat->Reset();
  gen_higgs->Reset();
  btagselected->Reset();
 basicLeptons.clear();
 basicHiggses.clear();
// basicvbf.clear();
//  for (unsigned int i=0; i<bl; i++) basicLeptons[i]->Reset();
//  for (unsigned int i=0; i<bh; i++) basicHiggses[i]->Reset();
  for (unsigned int i=0; i<bv; i++) basicvbf[i]->Reset();
  return 0;
}
///////////////////////////////////////////////////////////////////////////
// declare the histos
int decla(int mass){
   const char* label="i";

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
	//H1hist->GetYaxis()->SetTitle("Events              .");
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
		80, 0, 800);
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
		80, 20, 2000);
	j1histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	j1histpt->GetXaxis()->SetTitle("vbf j1 P_T (GeV)");
	basicvbf[0]= j1histpt; 

	TH1D *j1histeta = new TH1D("j1histeta",  
		label, 
		32, -6, 6);
	j1histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	j1histeta->GetXaxis()->SetTitle("vbf #eta_{j1} (GeV)");
	basicvbf[1]=j1histeta; 

	//

	TH1D *j2histpt = new TH1D("j2histpt",  
		label, 
		80, 20, 400);
	j2histpt->GetYaxis()->SetTitle("Events/ 2 GeV");
	j2histpt->GetXaxis()->SetTitle("vbf j2 P_T (GeV)");
	basicvbf[2]=j2histpt; 

	TH1D *j2histeta = new TH1D("j2histeta",  
		label, 
		32, -6, 6);
	j2histeta->GetYaxis()->SetTitle("Events/ 2 GeV");
	j2histeta->GetXaxis()->SetTitle("vbf #eta_{j2} (GeV)");
	basicvbf[3]=j2histeta; 

	//

	TH1D *RRadMassvbf = new TH1D("RadMass_ct4vbf",  
		label, 
		80, 50, 5000);
	RRadMassvbf->GetXaxis()->SetTitle("vbf M_{j j} (GeV)");
	basicvbf[4]=RRadMassvbf;

	//

	TH1D *RRadPtvbf = new TH1D("RadPt_ct4vbf",  
		label, 
		90, 50, 850);
	RRadPtvbf->GetXaxis()->SetTitle("vbf pt_{j j} (GeV)");
	basicvbf[5]=RRadPtvbf;

	//

	TH1D *Detavbf = new TH1D("Deta",  
		label, 
		30, 0, 10);
	Detavbf->GetXaxis()->SetTitle("vbf #Delta #eta");
	basicvbf[6]=Detavbf;
	
	//

	TH1D *Nfattagvbf = new TH1D("fattag_ct4vbf",  
		label, 
		6, -1.5, 4.5);
	Nfattagvbf->GetYaxis()->SetTitle("");
	Nfattagvbf->GetXaxis()->SetTitle("Number of fat tags vbf");
	basicvbf[7]=Nfattagvbf;

	TH1D *Nbtagvbf = new TH1D("btag_ct4vbf",  
		label, 
		7, -1.5, 5.5);
	Nbtagvbf->GetYaxis()->SetTitle("");
	Nbtagvbf->GetXaxis()->SetTitle("Number b-jets vbf");
	basicvbf[8]=Nbtagvbf;

	TH1D *Nbmistagvbf = new TH1D("bmistag_ct4vbf",  
		label, 
		7, -1.5, 5.5);
	Nbmistagvbf->GetYaxis()->SetTitle("");
	Nbmistagvbf->GetXaxis()->SetTitle("Number mistagged b jets's vbf");
	basicvbf[9]=Nbmistagvbf;

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
		30, 0, 6);
	DRll->GetXaxis()->SetTitle("#Delta R ll");
	basicLeptons.push_back (DRll);

	TH1D *DRbb = new TH1D("DRbb",  
		label, 
		30, 0, 6);
	DRbb->GetXaxis()->SetTitle("#Delta R bb");
	basicLeptons.push_back (DRbb);

	TH1D *MetMass1 = new TH1D("MetMass1_ct4",  
		label, 
		50, 0, 300);
	MetMass1->GetXaxis()->SetTitle("MET (GeV)");
	basicLeptons.push_back (MetMass1);

return 0;
}


