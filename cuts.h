/////////////////////////////////
// cuts
double const  weight =1./50000;
double HTVBF = 400;//400; 
double DeltayVBF = 3;//3;
bool shower=false;
// To be applied only to hadron level events
double const jet_ptmin=25.0; // We cut on all jets below 50 GeV
double const jet_ptminvbf=25.0; // We cut on all jets below 50 GeV
double const rapmax=5;
double const higgs_mass = 125.0;
double const etab = 2.5;
///////////////////////////////////
// for substructure
// mass drop
double const Rsb = 1.1; // CA recluster
double const mu = 0.67;
double const ycut = 0.09;
double const Mfat =110;
//
double const Rfilt = 0.3;
int const n_subjet =3;
///////////////////////////////////
// on the 4b's analysis
double const tolerance=10;
double const HThiggses = 0;
double const MHH = 200;
double const DetaHH = 10;//1.3;
double const DetaH = 15;//1.5;

