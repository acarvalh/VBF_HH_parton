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
double const jet1_ptminvbf=0.0; // We cut on all jets below 50 GeV
double const jet2_ptminvbf=0.0; // We cut on all jets below 50 Ge
double const bjetpt = 10;
double const etab = 2.5;
double const etaj=5;
double const higgs_mass = 125.0;
int const cat =3; // minimum number of btag
double const RR =0.001;
/////////////////////////////////////////
double const H1_ptmin=0.0; // We cut on all jets below 50 GeV
double const H2_ptmin=0.0; // We cut on all jets below 50 GeV
double const HH_ptmin=0.0; // We cut on all jets below 50 GeV
double const MHH = 0; // minimum
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
double const tolerance=10;
double const toleranceX=10;
double const HThiggses = 0;
double const DetaHH = 100;//1.3;
double const DetaH = 150;//1.5;
//////////////////////////////////
// on the wwbb analysis
double const ptlepton = 0.0;
double const lepiso = 0.0;
double const MeeMax = 5000.0;
double const MnunuMax = 5000.0;
