//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: PrimaryGeneratorAction.cc,v 1.6 2006-06-29 17:54:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

//#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "TH1F.h"
#include "TRandom3.h"

#include "CreateTree.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(const G4ThreeVector& posCentre , G4long myseed, ConfigFile config)
// : G4VUserPrimaryGeneratorAction(),
// fParticleGun(0)
{

  gRandom->SetSeed(myseed);
  primaries = config.read("primaries",1);
  int material = config.read("opticalMaterial",0); // 0 = GAGG_Ce_Mg - 1 = YAG_Ce

  if(material == 0)
  {
    MakeEnergyHisto_GaGG_Ce_Mg();
    MakeTimeHisto_GaGG_Ce_Mg();
  }
  else
  {
    MakeEnergyHisto_YAG_Ce();
    MakeTimeHisto_YAG_Ce();
  }


  G4GeneralParticleSource* gps = new G4GeneralParticleSource();

  //gps->GetCurrentSource()->GetEneDist()->SetMonoEnergy(0.511*MeV);
  //gps->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(0.,-1.,0.));
  // gps->GetCurrentSource()->GetPosDist()->SetCentreCoords(posCentre);

  gun = gps;
  // G4int n_particle = 1;
  // fParticleGun = new G4ParticleGun(n_particle);
  // //
  // G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  // G4ParticleDefinition* particle = particleTable->FindParticle("opticalphoton");
  // //
  // fParticleGun->SetParticleDefinition(particle);
  // fParticleGun->SetParticleTime(0.0*ns);
  // //
  // //
  // //
  // fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  // fParticleGun->SetParticleEnergy(3*eV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete gun;
  // delete fParticleGun;
}

void PrimaryGeneratorAction::MakeEnergyHisto_GaGG_Ce_Mg()
{
  const G4int N_GaGG_Ce_Mg = 132;
  G4double FAST_Energy_GaGG_Ce_Mg[N_GaGG_Ce_Mg]    =
  {2.6327, 2.6188, 2.60504, 2.59143, 2.57796, 2.56463, 2.55144, 2.53838,
2.52546, 2.51266, 2.5, 2.48746, 2.47505, 2.46276, 2.45059, 2.43854,
2.42661, 2.4148, 2.4031, 2.39151, 2.38004, 2.36867, 2.35741, 2.34626,
2.33522, 2.32427, 2.31343, 2.30269, 2.29205, 2.28151, 2.27106, 2.26071,
2.25045, 2.24029, 2.23022, 2.22023, 2.21034, 2.20053, 2.19081, 2.18118,
2.17163, 2.16216, 2.15278, 2.14347, 2.13425, 2.12511, 2.11604, 2.10705,
2.09814, 2.0893, 2.08054, 2.07185, 2.06323, 2.05468, 2.0462, 2.0378,
2.02946, 2.02119, 2.01299, 2.00485, 1.99678, 1.98877, 1.98083, 1.97295,
1.96513, 1.95738, 1.94969, 1.94205, 1.93448, 1.92696, 1.9195, 1.9121,
1.90476, 1.89748, 1.89024, 1.88307, 1.87595, 1.86888, 1.86186, 1.8549,
1.84799, 1.84113, 1.83432, 1.82756, 1.82085, 1.81419, 1.80758, 1.80102,
1.7945, 1.78803, 1.78161, 1.77523, 1.7689, 1.76262, 1.75637, 1.75018,
1.74402, 1.73791, 1.73184, 1.72582, 1.71983, 1.71389, 1.70799, 1.70213,
1.69631, 1.69052, 1.68478, 1.67908, 1.67341, 1.66779, 1.6622, 1.65665,
1.65113, 1.64565, 1.64021, 1.63481, 1.62943, 1.6241, 1.6188, 1.61353,
1.6083, 1.6031, 1.59794, 1.59281, 1.58771, 1.58264, 1.57761, 1.57261,
1.56764, 1.5627, 1.55779, 1.55291};


  G4double FAST_COMPONENT_GaGG_Ce_Mg[N_GaGG_Ce_Mg] =
  {0.0348132, 0.041346, 0.0462213, 0.055994, 0.0670797, 0.0809222, 0.0962054, 0.119476, 0.147153,
0.182414, 0.217745, 0.256009, 0.299913, 0.34846, 0.391284, 0.436086, 0.487627, 0.538087, 0.591007,
0.64409, 0.68982, 0.740803, 0.785007, 0.825553, 0.855418, 0.87689, 0.893165, 0.909735, 0.904133,
0.904835, 0.904515, 0.896044, 0.885449, 0.872733, 0.857353, 0.837832, 0.829267, 0.810443, 0.78659,
0.767397, 0.745877, 0.730203, 0.705547, 0.682982, 0.665159, 0.645723, 0.621589, 0.595471, 0.576994,
0.556964, 0.533997, 0.507302, 0.484761, 0.463904, 0.439846, 0.418253, 0.404643, 0.384661, 0.360237,
0.339777, 0.324817, 0.304963, 0.287939, 0.271009, 0.252333, 0.235159, 0.226375, 0.212185, 0.196952,
0.184789, 0.1707, 0.159992, 0.147724, 0.140217, 0.130929, 0.121987, 0.113251, 0.102937, 0.0964418,
0.0913118, 0.0863073, 0.0811948, 0.0732201, 0.0693917, 0.0642312, 0.0575242, 0.0533638, 0.0488039,
0.0474854, 0.0466519, 0.0428456, 0.0431906, 0.0402339, 0.0363488, 0.0340234, 0.0328062, 0.0313949,
0.0292326, 0.0259991, 0.0247748, 0.0240735, 0.0210389, 0.021114, 0.0211829, 0.0192813, 0.0187306,
0.0185142, 0.0171277, 0.0167855, 0.0169199, 0.016618, 0.0150835, 0.0137522, 0.0137759, 0.0131568,
0.0134337, 0.0135361, 0.0132448, 0.0122409, 0.0119099, 0.011421, 0.0115286, 0.0103885, 0.0105544,
0.0118189, 0.0109937, 0.0102842, 0.0103254, 0.0101991, 0.00931661, 0.00833321, 0.00805133};


  enHisto = new TH1F("enHisto","enHisto",10*N_GaGG_Ce_Mg,FAST_Energy_GaGG_Ce_Mg[0],FAST_Energy_GaGG_Ce_Mg[N_GaGG_Ce_Mg-1]);
  for(int en = 0; en < N_GaGG_Ce_Mg ; en++)
  {
    enHisto->Fill(FAST_Energy_GaGG_Ce_Mg[en],FAST_COMPONENT_GaGG_Ce_Mg[en]);
  }
  CreateTree::Instance()->enHisto = enHisto;

}

void PrimaryGeneratorAction::MakeEnergyHisto_YAG_Ce()
{
  //YAG_Ce
  const G4int N_YAG_Ce = 137;
  G4double FAST_Energy_YAG_Ce[N_YAG_Ce]    =
  {2.72827, 2.71335, 2.69859, 2.68398, 2.66954, 2.65525, 2.64111, 2.62712,
  2.61328, 2.59958,
    2.58603, 2.57261, 2.55934, 2.5462, 2.5332, 2.52033, 2.50758, 2.49497,
  2.48248, 2.47012,
    2.45788, 2.44576, 2.43376, 2.42188, 2.41011, 2.39845, 2.38691, 2.37548,
  2.36416, 2.35294,
    2.34183, 2.33083, 2.31993, 2.30912, 2.29842, 2.28782, 2.27732, 2.26691,
  2.2566, 2.24638,
    2.23625, 2.22621, 2.21626, 2.20641, 2.19663, 2.18695, 2.17735, 2.16783,
  2.1584, 2.14905,
    2.13978, 2.13058, 2.12147, 2.11244, 2.10348, 2.09459, 2.08579, 2.07705,
  2.06839, 2.0598,
    2.05128, 2.04283, 2.03445, 2.02614, 2.0179, 2.00972, 2.00161, 1.99357,
  1.98559, 1.97767,
    1.96982, 1.96203, 1.95429, 1.94662, 1.93901, 1.93146, 1.92397, 1.91654,
  1.90916, 1.90184,
    1.89458, 1.88737, 1.88021, 1.87311, 1.86606, 1.85907, 1.85213, 1.84524,
  1.8384, 1.83161,
    1.82487, 1.81818, 1.81154, 1.80495, 1.7984, 1.79191, 1.78546, 1.77905,
  1.77269, 1.76638,
    1.76011, 1.75389, 1.74771, 1.74157, 1.73548, 1.72943, 1.72342, 1.71745,
  1.71153, 1.70564,
    1.69979, 1.69399, 1.68822, 1.6825, 1.67681, 1.67116, 1.66555, 1.65997,
  1.65444, 1.64894,
    1.64347, 1.63804, 1.63265, 1.6273, 1.62198, 1.61669, 1.61144, 1.60622,
  1.60103, 1.59588,
    1.59076, 1.58568, 1.58062, 1.5756, 1.57061, 1.56566, 1.56073};

  G4double FAST_COMPONENT_YAG_Ce[N_YAG_Ce] =
  {0.017312, 0.0096453, 0.00790803, 0.00810917, 0.0088267, 0.0100443, 0.0119759, 0.013618, 0.0158952,
  0.0205216, 0.0276587,
    0.0355978, 0.047533, 0.063484, 0.0846373, 0.110553, 0.142739, 0.179262, 0.220629, 0.265465,
  0.316428, 0.373613, 0.435236,
    0.503049, 0.573143, 0.645777, 0.716553, 0.781489, 0.830571, 0.868259, 0.892522, 0.905547, 0.90921,
  0.904013, 0.894359,
    0.879165, 0.863961, 0.848755, 0.832279, 0.814811, 0.796789, 0.781315, 0.763563, 0.741355,
  0.724224, 0.705468, 0.687985,
    0.668306, 0.65034, 0.62812, 0.608418, 0.588743, 0.569351, 0.545681, 0.521421, 0.500633, 0.47682,
  0.455015, 0.432658,
    0.408202, 0.385464, 0.363256, 0.341688, 0.321938, 0.302229, 0.283179, 0.264294, 0.244777,
  0.227309, 0.208308, 0.193561,
    0.178659, 0.163398, 0.14938, 0.137934, 0.126634, 0.116962, 0.106936, 0.0974455, 0.0876783,
  0.0803322, 0.073087, 0.0664133,
    0.0601629, 0.0542353, 0.0497304, 0.0440716, 0.0397871, 0.035986, 0.0335943, 0.0305447, 0.027748,
  0.0251957, 0.0232941,
    0.0210163, 0.0181465, 0.0167443, 0.014946, 0.0138611, 0.0129008, 0.0123556, 0.0118624, 0.011176,
  0.0102985, 0.00944158,
    0.00946978, 0.00862405, 0.00771896, 0.00694318, 0.00643815, 0.00590943, 0.00537023, 0.00446022,
  0.00436113, 0.0043437,
    0.00408437, 0.00385623, 0.0032281, 0.00318978, 0.00319752, 0.00311753, 0.00256044, 0.00262745,
  0.00253195, 0.00259158,
    0.00233708, 0.00213198, 0.00234155, 0.0024483, 0.00209324, 0.0018311, 0.00175891, 0.00183815,
  0.00192849, 0.00163606,
    0.00136759, 0.00170809};

  enHisto = new TH1F("enHisto","enHisto",N_YAG_Ce,FAST_Energy_YAG_Ce[0],FAST_Energy_YAG_Ce[N_YAG_Ce-1]);
  for(int en = 0; en < N_YAG_Ce ; en++)
  {
    enHisto->Fill(FAST_Energy_YAG_Ce[en],FAST_COMPONENT_YAG_Ce[en]);
  }
  CreateTree::Instance()->enHisto = enHisto;

}

void PrimaryGeneratorAction::MakeTimeHisto_GaGG_Ce_Mg()
{
  G4double tr = 0.05;  //ns
  G4double td1 = 50.0; //ns
  G4double td2 = 200.0; //ns

  float pulseLength = 4000.0; //ns
  int steps = 4000000;
  tHisto = new TH1F("tHisto","tHisto",steps,0,pulseLength);// ns
  for(int t = 1; t < steps ; t++)
  {
    double tStamp = tHisto->GetBinCenter(t);
    double pvalue = (0.4*exp(-tStamp/td2)+0.6*exp(-tStamp/td1))*(1.0 - exp(-tStamp/tr));
    tHisto->SetBinContent(t,pvalue);
  }
  CreateTree::Instance()->tHisto = tHisto;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



void PrimaryGeneratorAction::MakeTimeHisto_YAG_Ce()
{

  G4double tr = 0.5;  //ns
  G4double td = 90.0; //ns
  float pulseLength = 1800.0; //ns
  int steps = 1800000;
  tHisto = new TH1F("tHisto","tHisto",steps,0,pulseLength);// ns
  for(int t = 0; t < steps ; t++)
  {
    double tStamp = tHisto->GetBinCenter(t);
    double pvalue = exp(-tStamp/td)*(1.0 - exp(-tStamp/tr));
    tHisto->SetBinContent(t,pvalue);
  }
  CreateTree::Instance()->tHisto = tHisto;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // std::cout << "New Primary" << std::endl;

  for(int p = 0 ; p < primaries ; p++)
  {
    gun->SetParticleTime((tHisto->GetRandom())*ns);
    // gun->GetCurrentSource()->GetEneDist()->SetMonoEnergy((2.0+(1.0/10000)*(p+1))*eV);
    gun->GetCurrentSource()->GetEneDist()->SetMonoEnergy((enHisto->GetRandom())*eV);
    gun->GeneratePrimaryVertex(anEvent);
  }

  // float sourcex = 0.;
  // float sourcey = 0.;
  // float sourcez = -30 *cm;
  // G4ThreeVector directionVector = G4ThreeVector(0.,0.,-1.);
  // //
  // //
  // fParticleGun->SetParticlePosition(G4ThreeVector(sourcex,sourcey,sourcez));
  // fParticleGun->SetParticleMomentumDirection(directionVector);
  // fParticleGun->GeneratePrimaryVertex(anEvent);
}
