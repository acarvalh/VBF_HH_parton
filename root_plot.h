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
#include <TH1D.h>
#include <vector>
using namespace std;
//int decla();
//int fill_hist(const vector<fastjet::PseudoJet> ,const vector<fastjet::PseudoJet> ,int);
//int save_hist();

	TH1D * basicvbf[10];
	vector<TH1D *> basicHiggses;
	vector<TH1D *> basicLeptons;
//	vector<TH1D *> basicvbf;
//	vector<TH1D *> basicHiggses;
//	vector<TH1D *> basicLeptons;
	TH1D *Njets_passing_kLooseID;
        TH1D *Cat;
        TH1D *gen_higgs;
        TH1D *btagselected;
