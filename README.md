VBF_HH_parton
=============

Les Houches

only needs root / fastjet

uses the output of the main11.cc commited here

.................................................................

source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.00/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh

mkdir histos

make HH_VBF (fix your root and fastjet paths on the Makefile)

./HH_VBF

root -l merge_hist.C

.................................................................

you can find a template file here
/afs/cern.ch/work/a/acarvalh/public/RSG_WBF_hh-Mhh260.lhe.decayed

if you do:

..................................................................

mkdir bulk_graviton

cp RSG_WBF_hh-Mhh260.lhe.decayed bulk_graviton

.................................................................

this recipe is out of envelope 




