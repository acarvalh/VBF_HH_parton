/////////////////////////////////
// cuts
double const  weight =1./20000;
double Mjj = 0;//400; 
double DeltayVBF = 0;//3;
bool shower=false;
// To be applied only to hadron level events
double const jet_ptmin=0.0; // for jet reconstruction
double const rapmax=50.0; // for jet reconstruction
double const jet_ptminvbf=0.0; // We cut on all jets below 50 GeV
double const bjetpt = 10;
double const etab = 2.5;
double const etaj=45;
double const higgs_mass = 125.0;
int const cat =3; // minimum number of btag
double const RR =0.5;
///////////////////////////////////
// for substructure
// mass drop
double const Rsb = 1.1; // CA recluster
double const mu = 0.67;
double const ycut = 0.09;
double const Mfat =120;
//
double const Rfilt = 0.1;
int const n_subjet =3;
///////////////////////////////////
// on the 4b's analysis
double const tolerance=100;
double const HThiggses = 0;
double const MHH = 0; // minimum
double const DetaHH = 100;//1.3;
double const DetaH = 150;//1.5;
//////////////////////////////////
// on the wwbb analysis
double const ptlepton = 20.0;
// hadronic dijet mass resolution
double const mass_resolution = 0.05;
