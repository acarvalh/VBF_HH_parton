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
  int points;
  if(resonant) {path="bulk_graviton/Madgraph0/MGraviton_"; points = masses;}
  else {path="nonresonant/pp_hh_vbf_"; points = parameters;}
  data = ".lhe.decayed";
  
  for(unsigned i=0;i<points;i++){
    if(resonant) {
      ostringstream o;
      o << "" << mass[i];//filenames
      file = path + o.str() + data;
    } else file = path + filenames[i] + data;
    filename.push_back(file);
    cout<<" "<<points<<" "<<i<<" "<<filename.at(i)<<endl;
  } // close filename
  //////////////////////////////////////////////////////////////////////////////////////
  ifstream in1;
  // declare the root plots to be done
  decla(0);
  for(int i=0;i<points;i++){ // for each mass
    cout<<"\n\n reading file = "<<filename.at(i)<<endl;
    in1.open(filename.at(i).c_str());
    for(unsigned ievent=0;ievent<nevent;ievent++){ // for each event  // 
//    for(unsigned ievent=0;ievent<2;ievent++){ // for each event  // 
        //cout<<"----------------------------------------------------"<<endl;
        string c;
        in1>>c;
        double Px, Py , Pz, E;
        int pID;
        unsigned nparticles;
        vector<PseudoJet> particles; 
        vector<PseudoJet> neutrinos;          
        vector<PseudoJet> leptons;
        vector<PseudoJet> photons;                    
        int nb = 0;
        int nvbf = 0;
        in1>>nparticles;  int counter=0; 
        for(unsigned ipart=0;ipart<nparticles;ipart++){ // loop on particles
          in1 >> pID >> Px >> Py >> Pz >> E ;//>> idup;
          //cout<< pID <<endl;   
          if (abs(pID) < 6 ){  // if a quark/gluon -- neglect hadrons
		particles.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); 
                particles.at(counter).set_user_index(pID); 
		if(abs(pID) == 5) nb++; else nvbf++; // count b's and no-b's
                //cout<<"particle flavour "<< particles.at(counter).user_index()<<endl;
		counter++;
	  }           
        } // finish to read all particles
        // construct the jets
        //cout<<"number of b particles = "<< nb<<endl;
        //cout<<"number of non-b particles = "<< nvbf<<endl;
	vector<PseudoJet> jets; 
        // check tags
	vector<int> btag; vector<int> bmistag; // number of btags/bmistags of each jet entry 
	vector<int> fattag;
        int njets = recojets(particles, jets,btag,bmistag);
	//cout<<"njets "<<njets<<endl; 
        istagged(jets, fattag); // check wheather the jet is fat -- by now only invariant mass
	// VBF combinatorics and cuts
	vector<int> vbftag; // keeps the entries of jets that are the VBF, have two entries
	bool VBF = findVBF(jets, fattag, btag, bmistag, vbftag);
	if (VBF){ // pass VBF
	  //vector<int> H1;
	  //vector<int> H2;
	  bool reco=false;
	  if(fourb) reco = analyse_4b(jets,fattag,btag,bmistag,vbftag);
        } // close pass VBF
    } // close for event
    if(resonant) save_hist(mass[i],resonant); else save_hist(parameter[i],resonant);
    in1.close();
  } // close for each mass
}
