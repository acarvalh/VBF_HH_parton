//////////////////
// to run:
// make HH_VBF
// ./HH_VBF
/////////////////
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <fstream>
using namespace fastjet;
using namespace std;
#include "Functions.h"
#include "choices.h"
int main() {
  srand( time(NULL) );
  hello();
  // input
  vector<string> filename;
  string file, path,data;
/*
  int points;
  if(resonant && !bkg) {
      //path="bulk_graviton/RSG_WBF_hh-Mhh"; points = masses;
      if(fourb) path="/afs/cern.ch/work/a/acarvalh/phenoHH/model_LHEfiles/spin0_hh_13tev/Mh2_bulk_"; points = masses;
      if(!fourb) path="bulk_graviton_mad_WWbb/MGraviton_"; points = masses;       
  }
  else if (!bkg && fourb) {path="nonresonant/pp_hh_vbf_"; points = parameters;}//nonresWWbb
  else if (!bkg && !fourb) {path="nonresWWbb/pp_hh_vbf_"; points = parameters;}//
  else {path="4bsbkg/"; points = components;}
  data = ".lhe.decayed.shower";
*/
//  path="bulk_graviton_mad/MGraviton_"; 
//  unsigned int points = masses;
//  data = ".lhe.decayed";
  path="spin0/Mh2_bulk_"; 
  unsigned int points = masses;
  data = ".lhe.decayed.shower";
  //  
  for(unsigned int i=0;i<points;i++){
    if(resonant && !bkg) {
      ostringstream o;
      o << "" << mass[i];//filenames
      file = path + o.str() + data;
    } else if(!bkg) file = path + filenames[i] + data;
    else file = path + bkgfilenames[i] + data;
    filename.push_back(file);
    cout<<" "<<points<<" "<<i<<" "<<filename.at(i)<<endl;
  } // close filename
  //////////////////////////////////////////////////////////////////////////////////////
  ifstream in1;
  // declare the root plots to be done
  decla(0);
  for(unsigned int i=1;i<2;i++){ // for each mass
  //  for(int i=0;i<2;i++){ // for each mass
  //  for(int i=8;i<9;i++){ // for each mass
  //  for(int i=11;i<15;i++){ // for each mass
    cout<<"\n\n reading file = "<<filename.at(i)<<endl;
    in1.open(filename.at(i).c_str());
    for(unsigned int ievent=0;ievent<nevent;ievent++){ // for each event  // 
    //for(int ievent=0;ievent<20;ievent++){ // for each event  // 
        //cout<<"----------------------------------------------------"<<endl;
        string c;
        in1>>c;
        double Px, Py , Pz, E;
        int pID;
        unsigned int nparticles;
        vector<PseudoJet> particles; 
        vector<PseudoJet> neutrinos;
        vector<PseudoJet> leptons;
        vector<PseudoJet> photons;   
        vector<PseudoJet> higgses;                    
        int nb = 0;
        int nvbf = 0;
        in1>>nparticles; unsigned int counter=0,counterh=0,counterl=0,countern=0;
   //cout<<nparticles<<" wtf! "<<endl; 
        for(unsigned int ipart=0;ipart<nparticles;ipart++){ // loop on particles
          in1 >> pID >> Px >> Py >> Pz >> E ;//>> idup;
          //cout<< pID <<endl;   
          if (abs(pID) < 6 || pID==21){  // if a quark/gluon -- neglect hadrons
		particles.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); 
                particles.at(counter).set_user_index(pID); 
		if(abs(pID) == 5) nb++; else nvbf++; // count b's and no-b's
                //cout<<"particle flavour "<< particles.at(counter).user_index()<<endl;
		counter++;
	  } else if (pID==25) {
		higgses.push_back(fastjet::PseudoJet(Px,Py,Pz,E));
		counterh++;
          } else if (abs(pID)==11 || abs(pID)==13) {
		leptons.push_back(fastjet::PseudoJet(Px,Py,Pz,E));
		counterl++;
          } else if (abs(pID)==12 || abs(pID)==14) {
		neutrinos.push_back(fastjet::PseudoJet(Px,Py,Pz,E));
		countern++;
          }
        } // finish to read all particles
        //if(nb!=4) {cout<<"not correct b quarks number "<<nb<<endl; exit(-10); }
	if(counterh==2) genhiggs(counterh, higgses);
        // construct the jets
        //cout<<"number of b particles = "<< nb<<" number of leptons = "<< counterl<<" number of neutrinos = "<< countern<<" number of jets = "<< nvbf<<" event "<< ievent<<endl;
        //cout<<"number of non-b particles = "<< nvbf<<endl;
	vector<PseudoJet> jets; 
        // check tags
	vector<int> btag; vector<int> bmistag; // number of btags/bmistags of each jet entry 
	vector<int> fattag;
        int njets = recojets(particles, jets,btag,bmistag,fattag);
        particles.clear();
	//cout<<"njets "<<njets<<endl; 
        //istagged(jets, fattag); // check wheather the jet is fat -- by now only invariant mass
	// VBF combinatorics and cuts
	vector<int> vbftag; // keeps the entries of jets that are the VBF, have two entries
	bool VBF =  findVBF(jets, fattag, btag, bmistag, vbftag);
        //cout<<"number of b particles = "<< nb<<" number of leptons = "<< counterl<<" number of neutrinos = "<< countern<<" number of jets = "<< nvbf<<" event "<< ievent<<endl;
	if (VBF){ // pass VBF
	  bool reco=false; 
	  if(fourb) {reco = analyse_4b(jets,fattag,btag,bmistag,vbftag,mass[i]); } //cout<<"here 3"<<endl;  
          else if (counterl>1 && countern>1) reco = analyse_2b2w(jets,fattag,btag,bmistag,vbftag,leptons,neutrinos); else cout<<"number of b particles = "<< nb<<" number of leptons = "<< counterl<<" number of neutrinos = "<< countern<<" number of jets = "<< nvbf<<" event "<< ievent<<endl;   
        } // close pass VBF
        //
    jets.clear();
    fattag.clear();
    btag.clear();
    } // close for event
    if(resonant && !bkg) save_hist(mass[i],resonant,bkg,fourb); 
    else if(!bkg) save_hist(parameter[i],resonant,bkg,fourb);
    else save_hist(comp[i],resonant,bkg,fourb);
    in1.close();
  } // close for each mass
}
