# FibresCalo
Geant4 simulation of a fibres "SPACAL" calorimeter


## Instructions for first time setup

### Instructions for installation on a personal PC

To donwnload and compile the package (assuming you have a local working release of Geant4 and CLHEP installed), execute the following commands:
```
git clone https://github.com/abenagli/FibresCalo.git
cd FibresCalo/build
cmake ../
make -j
cd -
```

### Instructions for installation on lxplus

To donwnload and compile the package on lxplus, these commands need to be run at each login:
```
. /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6/setup.sh
cd /afs/cern.ch/sw/lcg/external/geant4/10.3.p01/x86_64-slc6-gcc49-opt/bin/
source ./geant4.sh
cd -
```
The installation and compilation of the package then follows similarly:
```
git clone https://github.com/abenagli/FibresCalo.git
cd FibresCalo/build
../COMPILE_lxplus.sh 
make -j
cd -
```


## Instructions to run the code

To run the code in interactive mode (i.e. with visualization), simply execute
```
./build/FibresCalo template.cfg
```
This mode reads the visualization and beam parameters from file `vis.mac`

To produce multiple events and save the output in a root file, execute
```
./build/FibresCalo template.cfg outFileName
```
In this case, the beam parameters are read from file `gps.mac`

The file `template.cfg` contains a number of parameters to configure the calorimeter layout.
Most of them should be self-explanatory.
