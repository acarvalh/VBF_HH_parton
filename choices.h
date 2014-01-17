// inputs
unsigned const nevent=20000;
// Analysis
bool fourb=true;
bool resonant=false;
// for nonresonant
unsigned const parameters=9;
int parameter[parameters]={
1001,1000,101010,
51010,151010,
100010,102010,
101000,101020
};
string filenames[parameters]={
"SM_13tev_nocuts",
"SM_13tev_VBFcuts",
"BSM_13tev_VBFcuts_CV_p1p0_C2V_p1p0_C3_p1p0",
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
unsigned const masses=11;
int mass[masses]={260,300,350,400,450,
                  500,550,600,650,700,750};


/////////////////////////////////////

