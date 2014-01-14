#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include <TH1F.h>

// Test function
void hello();
///////////////////////////////////////////////////////////
int decla(int);
int fill_hist(const vector<fastjet::PseudoJet> ,const vector<fastjet::PseudoJet> ,
		const fastjet::PseudoJet ,const fastjet::PseudoJet);
int save_hist(int);
//////////////////////////////////////////////////////////
