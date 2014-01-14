//////////////////
// to run:
// make HH_VBF.cc
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
  path="/home/xanda/Documents/ggAnalysis/Parton/VBF_HH_parton/bulk_graviton/RSG_WBF_hh-Mhh";
  data = ".lhe.decayed";
  for(unsigned i=0;i<masses;i++){
    ostringstream o;
    o << "" << mass[i];
    file = path + o.str() + data;
    filename.push_back(file);
    cout<<filename.at(i)<<endl;
  } // close filename
  //////////////////////////////////////////////////////////////////////////////////////
  ifstream in1;
  for(int i=0;i<masses;i++){ // for each mass
    cout<<"\n\n reading file = "<<filename.at(i)<<endl;
    in1.open(filename.at(i).c_str());
    // declare the root plots to be done
    decla(mass[i]);
    for(unsigned ievent=0;ievent<nevent;ievent++){ // for each event  // 
//    for(unsigned ievent=0;ievent<20;ievent++){ // for each event  // 
        string c;
        in1>>c;
        double Px, Py , Pz, E;
        int pID;
        unsigned nparticles;
        vector<PseudoJet> vbfjets;
        vector<PseudoJet> bjets; 
        vector<PseudoJet> neutrinos;          
        vector<PseudoJet> leptons;
        vector<PseudoJet> photons;                    
        int nb = 0;
        int nvbf = 0;
        in1>>nparticles;  
        for(unsigned ipart=0;ipart<nparticles;ipart++){ // loop on particles
          in1 >> pID >> Px >> Py >> Pz >> E ;//>> idup;
          //cout<< pID <<endl;
          // photon counter      
          // if a b quark
          if (pID == 5 || pID == -5){bjets.push_back(fastjet::PseudoJet(Px,Py,Pz,E));nb++;} 
	  else {vbfjets.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); vbfjets.at(nvbf).set_user_index(pID); nvbf++;}           
        } // finish to read all particles
        if(nb!=4) cout<<"b-quarks = " <<nb <<" light-quarks = " <<nvbf <<endl;
        // test analysis
	vector<PseudoJet> H1;
	vector<PseudoJet> H2;
        for(int j=0;j<2;j++) H1.push_back(bjets.at(j));
        for(int j=2;j<4;j++) H2.push_back(bjets.at(j));
        // bvf jets
        PseudoJet j1 = vbfjets.at(0); 
	PseudoJet j2 = vbfjets.at(1);
	fill_hist(H1, H2, j1, j2);
    } // close for event
    save_hist(mass[i]);
  } // close for each mass

}
