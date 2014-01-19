//#ifndef root_plots_h
//#define root_plots_h
/*
#include <TROOT.h>
#include <TChain.h>

#include <TH2F.h>
#include <TFile.h>
#include <vector>
#include <iostream>
#include <fstream>
*/
#include <TH1F.h>

//int decla();
//int fill_hist(const vector<fastjet::PseudoJet> ,const vector<fastjet::PseudoJet> ,int);
//int save_hist();

	std::vector<TH1F *> basicvbf;
	std::vector<TH1F *> basicHiggses;
	TH1F *Njets_passing_kLooseID;
        TH1F *Cat;
