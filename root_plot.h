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

	std::vector<TH1F *> basic;
        TH1F* NJets;

	TH1F *Njets_passing_kLooseID;
	TH1F* DiPhotonMass;
	TH1F* DiJetMass;
	TH1F* TotalMass;

	TH1F* DeltaPhilead;
	TH1F* DeltaRlead;

	TH1F* PtPh1;
	TH1F* PtPh2;

	TH1F* EtaPh1;
	TH1F* EtaPh2;

	TH1F* PhiPh1;
	TH1F* PhiPh2;

	TH1F* PtJ1;
	TH1F* PtJ2;

	TH1F* EtaJ1;
	TH1F* EtaJ2;

	TH1F* PhiJ1;
	TH1F* PhiJ2;
	TH1F* DeltaRph;
	TH1F* leadsel;
	TH1F* subleadsel;
