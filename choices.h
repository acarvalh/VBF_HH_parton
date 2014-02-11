// inputs
// Analysis
bool fourb=true;
bool resonant=true;
bool bkg = false;
//
int const nevent=20000;
// for bkg
int const components =11;
int comp[components] = {0,1,2,3,4,5,6,7,8,9,10}; 
string bkgfilenames[components]={
"zbbbbjj_unweighted_events", // 10k
"zbbzbbjj_unweighted_events", // 10k
"dd__4bdd", // 10k
"gg__4bddbar", // 10k
"gg__4buubar", // 10k
"ud__4bud", // 10k
"BBBBjj_run3", // 10k
"uu__4buu", // 2.3k events
"4b2j_alpgen", //500 k
"t_wbt_wbjj", // 100kevents
"wwbbjj_ptb20" // 10k
};
// for nonresonant
int const parameters=9;
int parameter[parameters]={
101010,1000,1001,
51010,151010,
100010,102010,
101000,101020
};
string filenames[parameters]={
"BSM_13tev_VBFcuts_CV_p1p0_C2V_p1p0_C3_p1p0",
"SM_13tev_VBFcuts",
"SM_13tev_nocuts",
//
"BSM_13tev_VBFcuts_CV_p0p5_C2V_p1p0_C3_p1p0",
"BSM_13tev_VBFcuts_CV_p1p5_C2V_p1p0_C3_p1p0",
//
"BSM_13tev_VBFcuts_CV_p1p0_C2V_p0p0_C3_p1p0",
"BSM_13tev_VBFcuts_CV_p1p0_C2V_p2p0_C3_p1p0",
//
"BSM_13tev_VBFcuts_CV_p1p0_C2V_p1p0_C3_p0p0",
"BSM_13tev_VBFcuts_CV_p1p0_C2V_p1p0_C3_p2p0"
};
// for resonant
int const masses=26;
int mass[masses]={260,300 ,350 ,400 ,450 ,500 ,
		      550 ,600 ,650 ,700 ,750,
		      800 ,850 ,900 ,950 ,1000,
                      1050,1100,1150,1200,1250,
		      1300,1350,1400,1450,1500
                 };

//unsigned const masses=2;
//int mass[masses]={260,500};
/////////////////////////////////////

