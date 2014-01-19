#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <TH1F.h>

// Test function
void hello();
///////////////////////////////////////////////////////////
// histos
int decla(int);
int fill_hist(const vector<fastjet::PseudoJet> ,const vector<fastjet::PseudoJet> ,
		const fastjet::PseudoJet ,const fastjet::PseudoJet);
int save_hist(int,bool);
//////////////////////////////////////////////////////////
// tags
int recojets(vector<PseudoJet> particles,vector<PseudoJet> & jets, vector<int> & btag, vector<int> & bmistag, vector<int> & fattag);
//void istagged(vector<PseudoJet> jets, vector<int> & fattag);
bool findVBF(vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag,vector<int> & vbftag);
/////////////////////////////////////////////////////////////
// analysis
bool analyse_4b(
	vector<PseudoJet> jets, vector<int> fattag, vector<int> btag, vector<int> bmistag,vector<int> vbftag);
bool jets_semi_hadronic( ); // deal with the 2 jets of a bbXX decay, where X is not q/g
bool analyse_2b2w( );
////////////////////////////////////////////////////////////
