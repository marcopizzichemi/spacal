#ifndef CreateTree_H
#define CreateTree_H 1

#include <iostream>
#include <vector>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"



class CreateTree
{
private:

  TTree*  ftree ;
  TString fname ;

public:

  CreateTree (TString name,
              const std::vector<float>& attL) ;
  ~CreateTree () ;

  TTree*             GetTree  () const { return ftree ; } ;
  TString            GetName  () const { return fname ; } ;
  void               AddEnergyDeposit (int index, float deposit, std::map<int,float>& depositAtt) ;
  // void               AddEnergyDepositABS (int Aindex, float depositA, std::map<int,float>) ;
 // void               AddEnergyDeposit_1st_Section (int indexion, float deposit) ;
 // void               AddEnergyDeposit_2nd_Section (int index2, float deposit) ;
  void               AddScintillationPhoton (int index) ;
  void               AddCerenkovPhoton (int index) ;
  int                Fill     () ;
  bool               Write    (TFile *) ;
  void               Clear    () ;
  static CreateTree* Instance () { return fInstance ; } ;
  static CreateTree* fInstance ;

  int   Event ;

  std::vector<float> * inputMomentum ; // Px Py Pz E
  std::vector<float> * inputInitialPosition ; // x, y, z

  float depositedEnergyTotal ;
  float depositedEnergy_1st_Section ;
  float depositedEnergy_2nd_Section ;
  float depositedEnergyFibres ;
  float depositedEnergyFibres_1st_Section ;
  float depositedEnergyFibres_2nd_Section ;
  float depositedEnergyFibresCross;
  float depositedEnergyFibresCenter;
  float depositedEnergyFibresCorners;
  float depositedEnergy_2nd_Sect_FibresCross;
  float depositedEnergy_2nd_Sect_FibresCenter;
float depositedEnergy_2nd_Sect_FibresCorners;
  float depositedEnergyAbsorber ;
  float depositedEnergyCell1 ;
  float depositedEnergyCell2 ;
  float depositedEnergyCell3 ;
  float depositedEnergyCell4 ;
  float depositedEnergyCell5 ;
  float depositedEnergyCell6 ;
  float depositedEnergyCell7 ;
  float depositedEnergyCell8 ;
  float depositedEnergyCell9 ;
  float depositedEnergyCell10 ;
  float depositedEnergyCell11 ;
  float depositedEnergyCell12 ;
  float depositedEnergyCell13 ;
  float depositedEnergyCell14 ;
  float depositedEnergyCell15 ;
  float depositedEnergyCell16 ;
  float depositedEnergyCell17 ;
  float depositedEnergyCell18 ;
  float depositedEnergyAbsorber_1st_Section ;
  float depositedEnergyAbsorber_2nd_Section ;
  float depositedEnergyFibres_post ;
  float depositedEnergyAbsorber_post ;
  float depositedEnergyWorld ;
  std::vector<float> * depositedEnergyFibresAtt ;

  float totalTrackLengthFibres ;
  float totalTrackLengthOverThFibres ;

  int tot_phot_cer;
  int tot_phot_cer_post;
  int tot_gap_phot_cer;
  int tot_det_phot_cer;
  std::vector<float> * tot_gap_photFast_cer;
  std::vector<float> * tot_det_photFast_cer;

  // energy deposited in each fibre of a tower
  std::vector<float> * depositedEnergies ;
  //std::vector<float> * depositedEnergiesABS ;
  // std::vector<float> * depositedEnergies_1st_Section;
  // std::vector<float> * depositedEnergies_2nd_Section;
  std::vector<std::vector<float> > * depositedEnergiesAtt ;
  // index of the fibre where the deposit happens
  std::vector<int> * depositFibres ;
  //std::vector<int> * depositAbsorber ;

  // std::vector<int> * depositFibres_1st_Section;
  // std::vector<int> * depositFibres_2nd_Section;

  // scintillation photons produced in each fibre of a tower
  std::vector<int> * scintillationPhotons ;
  // index of the fibre where the deposit happens
  std::vector<int> * scintillationFibres ;

  // cerenkov photons produced in each fibre of a tower
  std::vector<int> * cerenkovPhotons ;
  // index of the fibre where the deposit happens
  std::vector<int> * cerenkovFibres ;

  float Radial_stepLength;
  float Longitudinal_stepLength;
  float Radial_ion_energy_absorber[5000];
  float Longitudinal_ion_energy_absorber[5000];

  float PrimaryParticleX[1000];
  float PrimaryParticleY[1000];
  float PrimaryParticleZ[1000];
  float PrimaryParticleE[1000];

  // histograms
  //TH1F* h_phot_cer_lambda;
  //TH1F* h_phot_cer_E;
  //TH1F* h_phot_cer_time;

  //TH1F* h_phot_cer_gap_lambda;
  //TH1F* h_phot_cer_gap_E;
  //TH1F* h_phot_cer_gap_time;

  //std::vector<TH1F*> h_photFast_cer_gap_lambda;
  //std::vector<TH1F*> h_photFast_cer_gap_E;
  //std::vector<TH1F*> h_photFast_cer_gap_time;

  // to be filled at the beginning of the event generation only
  TNtuple * fibresPosition ;
  TNtuple * fibresPosition_1st_Section;
  TNtuple * fibresPosition_2nd_Section;
  TTree * attenuationLengths ;

  TTree * photons;
  int   eventNumber      ;
  int   front_back      ;
  int   pmt_number      ;
  int   module_number       ;

  int   phPerPMT[2][9];

  float vertX         ;
  float vertY         ;
  float vertZ         ;
  float PositionX     ;
  float PositionY     ;
  float PositionZ     ;
  float PreMomentumX  ;
  float PreMomentumY  ;
  float PreMomentumZ  ;
  float PostMomentumX ;
  float PostMomentumY ;
  float PostMomentumZ ;
  float globalTime    ;
  float PhotonEnergy  ;

  TTree * photonsAbsPoint;
  // int abs_eventNumber;
  float abs_x         ;
  float abs_y          ;
  float abs_z           ;
  float abs_globalTime   ;
  float abs_PhotonEnergy  ;

  std::vector<float> * attLengths;

  TH1F *enHisto;
  TH1F  *tHisto;
};

#endif
