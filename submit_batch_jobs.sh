rm -rf results/output*;
cd condor/;
root -q -b -l condor.C;
sh submit.sh >> /dev/null;
