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
  path="bulk_graviton/kkgraviton_cg";
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
        cout<<"----------------------------------------------------"<<endl;
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
		particles.push_back(fastjet::PseudoJet(Px,Py,Pz,E)); particles.at(counter).set_user_index(pID); 
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
	cout<<"njets "<<njets<<endl; 
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
    save_hist(mass[i]);
  } // close for each mass

}
