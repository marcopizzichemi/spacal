#include "CreateTree.hh"
#include <algorithm>

using namespace std ;

CreateTree* CreateTree::fInstance = NULL ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::CreateTree (TString name,
                        const std::vector<float>& attL)
{
  if ( fInstance )
  {
    return ;
  }

  this->fInstance = this ;
  this->fname     = name ;
  this->ftree     = new TTree (name,name) ;


  this->attLengths = new std::vector<float>;
  for(unsigned int it = 0; it < attL.size(); ++it)
    this->attLengths->push_back( attL.at(it) );

  this->GetTree ()->Branch ("Event", &this->Event, "Event/I") ;


  inputMomentum = new vector<float> (4, 0.) ;
  inputInitialPosition = new vector<float> (3, 0.) ;
  depositedEnergyFibresAtt = new vector<float> ();
  depositedEnergies = new vector<float> () ;
  //  depositedEnergiesABS = new vector<float> () ;
  // depositedEnergies_1st_Section = new vector<float> () ;
  // depositedEnergies_2nd_Section = new vector<float> () ;
  depositedEnergiesAtt = new vector<vector<float> > () ;
  tot_gap_photFast_cer = new vector<float> ();
  tot_det_photFast_cer = new vector<float> ();
  // depositFibres_1st_Section = new vector<int> () ;
  // depositFibres_2nd_Section = new vector<int> () ;
  depositFibres = new vector<int> () ;
  // depositAbsorber = new vector<int> () ;
  cerenkovPhotons = new vector<int> () ;
  cerenkovFibres = new vector<int> () ;
  scintillationPhotons = new vector<int> () ;
  scintillationFibres = new vector<int> () ;

  this->GetTree ()->Branch ("inputMomentum","vector<float>",&inputMomentum) ;
  this->GetTree ()->Branch ("inputInitialPosition","vector<float>",&inputInitialPosition) ;

  this->GetTree ()->Branch ("depositedEnergyTotal",        &this->depositedEnergyTotal,                "depositedEnergyTotal/F") ;

  this->GetTree ()->Branch ("depositedEnergy_1st_Section",        &this->depositedEnergy_1st_Section,                "depositedEnergy_1st_Section/F") ;

  this->GetTree ()->Branch ("depositedEnergy_2nd_Section",        &this->depositedEnergy_2nd_Section,                "depositedEnergy_2nd_Section/F") ;



  this->GetTree ()->Branch ("depositedEnergyFibres",       &this->depositedEnergyFibres,              "depositedEnergyFibres/F") ;

  //  ************    *****************   **********
  this->GetTree ()->Branch ("depositedEnergyFibres_1st_Section",       &this->depositedEnergyFibres_1st_Section,              "depositedEnergyFibres_1st_Section/F") ;

  this->GetTree ()->Branch ("depositedEnergyFibres_2nd_Section",       &this->depositedEnergyFibres_2nd_Section,              "depositedEnergyFibres_2nd_Section/F") ;


  // ************* ************   ***************

  //////////////////////////////////////////////
 this->GetTree ()->Branch ("depositedEnergyFibresCross",       &this->depositedEnergyFibresCross,              "depositedEnergyFibresCross/F") ;
  this->GetTree ()->Branch ("depositedEnergyFibresCenter",       &this->depositedEnergyFibresCenter,              "depositedEnergyFibresCenter/F") ;
 this->GetTree ()->Branch ("depositedEnergyFibresCorners",       &this->depositedEnergyFibresCorners,              "depositedEnergyFibresCorners/F") ;

  // **********   **********   *******************
 this->GetTree ()->Branch ("depositedEnergy_2nd_Sect_FibresCross",       &this->depositedEnergy_2nd_Sect_FibresCross,              "depositedEnergy_1st_Sect_FibresCross/F") ;
  this->GetTree ()->Branch ("depositedEnergy_2nd _Sect_FibresCenter",       &this->depositedEnergy_2nd_Sect_FibresCenter,              "depositedEnergy_2nd_Sect_FibresCenter/F") ;
 this->GetTree ()->Branch ("depositedEnergy_2nd _Sect_FibresCorners",       &this->depositedEnergy_2nd_Sect_FibresCorners,              "depositedEnergy_2nd_Sect_FibresCorners/F") ;


  // ***********   ************  *****************

  ////////////////////////////////////////////////////////////////
  this->GetTree ()->Branch ("depositedEnergyAbsorber",     &this->depositedEnergyAbsorber,          "depositedEnergyAbsorber/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell1",     &this->depositedEnergyCell1,          "depositedEnergyCell1/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell2",     &this->depositedEnergyCell2,          "depositedEnergyCell2/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell3",     &this->depositedEnergyCell3,          "depositedEnergyCell3/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell4",     &this->depositedEnergyCell4,          "depositedEnergyCell4/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell5",     &this->depositedEnergyCell5,          "depositedEnergyCell5/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell6",     &this->depositedEnergyCell6,          "depositedEnergyCell6/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell7",     &this->depositedEnergyCell7,          "depositedEnergyCell7/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell8",     &this->depositedEnergyCell8,          "depositedEnergyCell8/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell9",     &this->depositedEnergyCell9,          "depositedEnergyCell9/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell10",     &this->depositedEnergyCell10,          "depositedEnergyCell10/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell11",     &this->depositedEnergyCell11,          "depositedEnergyCell11/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell12",     &this->depositedEnergyCell12,          "depositedEnergyCell12/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell13",     &this->depositedEnergyCell13,          "depositedEnergyCell13/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell14",     &this->depositedEnergyCell14,          "depositedEnergyCell14/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell15",     &this->depositedEnergyCell15,          "depositedEnergyCell15/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell16",     &this->depositedEnergyCell16,          "depositedEnergyCell16/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell17",     &this->depositedEnergyCell17,          "depositedEnergyCell17/F") ;
  this->GetTree ()->Branch ("depositedEnergyCell18",     &this->depositedEnergyCell18,          "depositedEnergyCell18/F") ;

  this->GetTree ()->Branch ("depositedEnergyAbsorber_1st_Section",     &this->depositedEnergyAbsorber_1st_Section,          "depositedEnergyAbsorber_1st_section/F") ;
  this->GetTree ()->Branch ("depositedEnergyAbsorber_2nd_Section",     &this->depositedEnergyAbsorber_2nd_Section,          "depositedEnergyAbsorber_2nd_Section/F") ;
  this->GetTree ()->Branch ("depositedEnergyFibres_post",  &this->depositedEnergyFibres_post,    "depositedEnergyFibres_post/F") ;
  this->GetTree ()->Branch ("depositedEnergyAbsorber_post",&this->depositedEnergyAbsorber_post,"depositedEnergyAbsorber_post/F") ;
  this->GetTree ()->Branch ("depositedEnergyWorld",        &this->depositedEnergyWorld,       "         depositedEnergyWorld/F") ;
  this->GetTree ()->Branch ("depositedEnergyFibresAtt", "vector<float>",&depositedEnergyFibresAtt);

  this->GetTree ()->Branch ("totalTrackLengthFibres",       &this->totalTrackLengthFibres,             "totalTrackLengthFibres/F") ;
  this->GetTree ()->Branch ("totalTrackLengthOverThFibres", &this->totalTrackLengthOverThFibres, "totalTrackLengthOverThFibres/F") ;

  this->GetTree ()->Branch ("tot_phot_cer",         &this->tot_phot_cer,                 "tot_phot_cer/I") ;
  this->GetTree ()->Branch ("tot_phot_cer_post",    &this->tot_phot_cer_post,       "tot_phot_cer_post/I") ;
  this->GetTree ()->Branch ("tot_gap_phot_cer",     &this->tot_gap_phot_cer,         "tot_gap_phot_cer/I") ;
  this->GetTree ()->Branch ("tot_det_phot_cer",     &this->tot_det_phot_cer,         "tot_det_phot_cer/I") ;
  this->GetTree ()->Branch ("tot_gap_photFast_cer", "vector<float>",&this->tot_gap_photFast_cer) ;
  this->GetTree ()->Branch ("tot_det_photFast_cer", "vector<float>",&this->tot_det_photFast_cer) ;

  // this->GetTree ()->Branch ("depositedEnergies_1st_Section","vector<float>",&depositedEnergies_1st_Section) ;
  // this->GetTree ()->Branch ("depositedEnergies_2nd_Section","vector<float>",&depositedEnergies_2nd_Section) ;
  this->GetTree ()->Branch ("depositedEnergies","vector<float>",&depositedEnergies) ;
  // this->GetTree ()->Branch ("depositedEnergiesABS","vector<float>",&depositedEnergiesABS) ;
  this->GetTree ()->Branch ("depositedEnergiesAtt","vector<vector<float> >",&depositedEnergiesAtt) ;
  this->GetTree ()->Branch ("depositFibres","vector<int>",&depositFibres) ;
  //  this->GetTree ()->Branch ("depositAbsorber","vector<int>",&depositAbsorber) ;
 // this->GetTree ()->Branch ("depositFibres_1st_Section","vector<int>",&depositFibres_1st_Section) ;
 // this->GetTree ()->Branch ("depositFibres_2nd_Section","vector<int>",&depositFibres_2nd_Section) ;

  this->GetTree ()->Branch ("cerenkovPhotons","vector<int>",&cerenkovPhotons) ;
  this->GetTree ()->Branch ("cerenkovFibres","vector<int>",&cerenkovFibres) ;

  this->GetTree ()->Branch ("scintillationPhotons","vector<int>",&scintillationPhotons) ;
  this->GetTree ()->Branch ("scintillationFibres","vector<int>",&scintillationFibres) ;

  this->GetTree ()->Branch ("Radial_stepLength",               &Radial_stepLength,                                     "Radial_stepLength/F");
  this->GetTree ()->Branch ("Longitudinal_stepLength",         &Longitudinal_stepLength,                         "Longitudinal_stepLength/F");
  this->GetTree ()->Branch ("Radial_ion_energy_absorber",       Radial_ion_energy_absorber,             "Radial_ion_energy_absorber[5000]/F");
  this->GetTree ()->Branch ("Longitudinal_ion_energy_absorber", Longitudinal_ion_energy_absorber, "Longitudinal_ion_energy_absorber[5000]/F");

  this->GetTree()->Branch("PrimaryParticleX",PrimaryParticleX,"PrimaryParticleX[1000]/F");
  this->GetTree()->Branch("PrimaryParticleY",PrimaryParticleY,"PrimaryParticleY[1000]/F");
  this->GetTree()->Branch("PrimaryParticleZ",PrimaryParticleZ,"PrimaryParticleZ[1000]/F");
  this->GetTree()->Branch("PrimaryParticleE",PrimaryParticleE,"PrimaryParticleE[1000]/F");


  //h_phot_cer_lambda = new TH1F("h_phot_cer_lambda","",1000,250.,1250.);
  //h_phot_cer_E = new TH1F("h_phot_cer_E","",1000,0.,5.);
  //h_phot_cer_time = new TH1F("h_phot_cer_time","",10000,0.,10000.);

  //h_phot_cer_gap_lambda = new TH1F("h_phot_cer_gap_lambda","",1000,250.,1250.);
  //h_phot_cer_gap_E = new TH1F("h_phot_cer_gap_E","",1000,0.,5.);
  //h_phot_cer_gap_time = new TH1F("h_phot_cer_gap_time","",10000,0.,10000.);

  //for(unsigned int it = 0; it < this->attLengths.size(); ++it)
  //{
  //  int attLength = int( attLengths.at(it) );
  //  h_photFast_cer_gap_lambda[attLength] = new TH1F(Form("h_photFast_cer_gap_lambda_attLength%04d",attLength),"",1000,250.,1250.);
  //  h_photFast_cer_gap_E[attLength] = new TH1F(Form("h_photFast_cer_gap_E_attLength%04d",attLength),"",1000,0.,5.);
  //  h_photFast_cer_gap_time[attLength] = new TH1F(Form("h_photFast_cer_gap_time_attLength%04d",attLength),"",10000,0.,10000.);
  //}

  fibresPosition = new TNtuple ("fibresPosition", "fibresPosition", "N:x:y:z") ;
  fibresPosition_1st_Section = new TNtuple ("fibresPosition_1st_Section", "fibresPosition_1st_Section", "N:x:y") ;
  fibresPosition_2nd_Section = new TNtuple ("fibresPosition_2nd_Section", "fibresPosition_2nd_Section", "N:x:y") ;

  photons = new TTree ("photons", "photons");
  // photons->Branch("PrimaryParticleX",PrimaryParticleX,"PrimaryParticleX[1000]/F");

  photons->Branch("event"         , &this->eventNumber   , "event/I");
  photons->Branch("front_back"    , &this->front_back    , "front_back/I");
  photons->Branch("pmt_number"    , &this->pmt_number    , "pmt_number/I");
  photons->Branch("module_number"    , &this->module_number    , "module_number/I");
  photons->Branch("phPerPMT"      , &this->phPerPMT      , "phPerPMT/I");
  photons->Branch("vertX"         , &this->vertX         , "vertX/F");
  photons->Branch("vertY"         , &this->vertY         , "vertY/F");
  photons->Branch("vertZ"         , &this->vertZ         , "vertZ/F");
  photons->Branch("PositionX"     , &this->PositionX     , "PositionX/F");
  photons->Branch("PositionY"     , &this->PositionY     , "PositionY/F");
  photons->Branch("PositionZ"     , &this->PositionZ     , "PositionZ/F");
  photons->Branch("PreMomentumX"  , &this->PreMomentumX  , "PreMomentumX/F");
  photons->Branch("PreMomentumY"  , &this->PreMomentumY  , "PreMomentumY/F");
  photons->Branch("PreMomentumZ"  , &this->PreMomentumZ  , "PreMomentumZ/F");
  photons->Branch("PostMomentumX" , &this->PostMomentumX , "PostMomentumX/F");
  photons->Branch("PostMomentumY" , &this->PostMomentumY , "PostMomentumY/F");
  photons->Branch("PostMomentumZ" , &this->PostMomentumZ , "PostMomentumZ/F");
  photons->Branch("globalTime"    , &this->globalTime    , "globalTime/F");
  photons->Branch("PhotonEnergy"  , &this->PhotonEnergy  , "PhotonEnergy/F");

  photonsAbsPoint = new TTree ("photonsAbsPoint", "photonsAbsPoint");
  photonsAbsPoint->Branch("event"         , &this->eventNumber   , "event/I");
  photonsAbsPoint->Branch("x"             , &this->abs_x             , "x/F");
  photonsAbsPoint->Branch("y"             , &this->abs_y             , "y/F");
  photonsAbsPoint->Branch("z"             , &this->abs_z             , "z/F");
  photonsAbsPoint->Branch("globalTime"    , &this->abs_globalTime    , "globalTime/F");
  photonsAbsPoint->Branch("PhotonEnergy"  , &this->abs_PhotonEnergy  , "PhotonEnergy/F");

  attenuationLengths = new TTree("attenuationLengths", "attenuationLenghts");
  attenuationLengths -> Branch("attLengths","vector<float>",&attLengths);


  this->Clear () ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::~CreateTree ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void
CreateTree::AddEnergyDeposit (int index, float deposit, std::map<int,float>& depositAtt)
{
  // find if it exists already
  vector<int>::const_iterator where = find (depositFibres->begin (),
                                            depositFibres->end (), index) ;

  if (depositFibres->end () == where)
  {
    depositFibres->push_back (index) ;
    depositedEnergies->push_back (deposit) ;
    int i = 0;
    for(std::map<int,float>::const_iterator it = depositAtt.begin(); it != depositAtt.end(); ++it)
    {
      (depositedEnergiesAtt->at(i)).push_back( (depositAtt[it->first]) ) ;
      ++i;
    }
  }
  else
  {
    depositedEnergies->at (where - depositFibres->begin ()) += deposit ;
    int i = 0;
    for(std::map<int,float>::const_iterator it = depositAtt.begin(); it != depositAtt.end(); ++it)
    {
      (depositedEnergiesAtt->at(i)).at (where - depositFibres->begin ()) += depositAtt[it->first] ;
      ++i;
    }
  }

  return ;
}
// ********************************************************
// void
// CreateTree::AddEnergyDepositABS (int Aindex, float depositA, std::map<int,float>)
// {
//   // find if it exists already
//   vector<int>::const_iteratorA where = find (depositAbsorber->begin (),
//                                             depositAbsorber->end (), Aindex) ;

//   if (depositAbsorber->end () == where)
//   {
//     depositAbsorber->push_back (Aindex) ;
//     depositedEnergiesABS->push_back (depositA) ;
//     int i = 0;

//   }
//   else
//   {
//     depositedEnergiesABS->at (where - depositAbsorber->begin ()) += depositA ;
//     int i = 0;

//   }

//   return ;
//}
// ********************************************************

  // ***************  *******************  *****************

// void
// CreateTree::AddEnergyDeposit_1st_Section (int indexion, float deposit1)
// {
//   // find if it exists already
//   vector<int>::const_iterator where = find (depositFibres_1st_Section->begin (),
//                                             depositFibres_1st_Section->end (), indexion) ;

//   if (depositFibres_1st_Section->end () == where)
//   {
//     depositFibres_1st_Section->push_back (indexion) ;
//     depositedEnergies_1st_Section->push_back (deposit1) ;

//   }
//   return;
// }

// void
// CreateTree::AddEnergyDeposit_2nd_Section (int index2, float deposit2)
// {
//  vector<int>::const_iterator where = find (depositFibres_2nd_Section->begin (),
//                                             depositFibres_2nd_Section->end (), index2) ;

//   if (depositFibres_2nd_Section->end () == where)
//   {
//     depositFibres_2nd_Section->push_back (index2) ;
//     depositedEnergies_2nd_Section->push_back (deposit2) ;

//   }
//   return;
// }

  // *************** ******************** ****************


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void
CreateTree::AddScintillationPhoton (int index)
{
  // find if it exists already
  vector<int>::const_iterator where = find (scintillationFibres->begin (),
                                            scintillationFibres->end (), index) ;
  if (scintillationFibres->end () == where)
    {
      scintillationFibres->push_back (index) ;
      scintillationPhotons->push_back (1) ;
    }
  else
    {
      scintillationPhotons->at (where - scintillationFibres->begin ()) += 1 ;
    }
  return ;
}



// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void
CreateTree::AddCerenkovPhoton (int index)
{
  // find if it exists already
  vector<int>::const_iterator where = find (cerenkovFibres->begin (),
                                            cerenkovFibres->end (), index) ;
  if (cerenkovFibres->end () == where)
    {
      cerenkovFibres->push_back (index) ;
      cerenkovPhotons->push_back (1) ;
    }
  else
    {
      cerenkovPhotons->at (where - cerenkovFibres->begin ()) += 1 ;
    }
  return ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int CreateTree::Fill ()
{
  return this->GetTree ()->Fill () ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


bool CreateTree::Write (TFile * outfile)
{
  outfile->cd () ;
  ftree->Write () ;
  fibresPosition->Write () ;
  fibresPosition_1st_Section->Write () ;
  fibresPosition_2nd_Section->Write () ;
  attenuationLengths->Write() ;
  photons->Write();
  photonsAbsPoint->Write();
  enHisto->Write();
  tHisto->Write();
  //h_phot_cer_lambda->Write();
  //h_phot_cer_E->Write();
  //h_phot_cer_time->Write();
  //h_phot_cer_gap_lambda->Write();
  //h_phot_cer_gap_E->Write();
  //h_phot_cer_gap_time->Write();
  //for(unsigned int it = 0; it < attLengths.size(); ++it)
  //{
  //  int attLength = int( attLengths.at(it) );
  //  h_photFast_cer_gap_lambda[attLength]->Write();
  //  h_photFast_cer_gap_E[attLength]->Write();
  //  h_photFast_cer_gap_time[attLength]->Write();
  //}
  return true ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void CreateTree::Clear ()
{
  Event	= 0 ;

  depositedEnergyTotal = 0. ;
  depositedEnergyFibres = 0. ;
  depositedEnergyFibresCross = 0. ;
  depositedEnergyFibresCenter = 0. ;
  depositedEnergyFibresCorners = 0. ;
  depositedEnergyAbsorber = 0. ;
 depositedEnergyCell1 = 0. ;
 depositedEnergyCell2 = 0. ;
 depositedEnergyCell3 = 0. ;
 depositedEnergyCell4 = 0. ;
 depositedEnergyCell5 = 0. ;
 depositedEnergyCell6 = 0. ;
 depositedEnergyCell7 = 0. ;
 depositedEnergyCell8 = 0. ;
 depositedEnergyCell9 = 0. ;
 depositedEnergyCell10 = 0. ;
 depositedEnergyCell11 = 0. ;
 depositedEnergyCell12 = 0. ;
 depositedEnergyCell13 = 0. ;
 depositedEnergyCell14 = 0. ;
 depositedEnergyCell15 = 0. ;
 depositedEnergyCell16 = 0. ;
 depositedEnergyCell17 = 0. ;
 depositedEnergyCell18 = 0. ;
  depositedEnergyFibres_post = 0. ;
  depositedEnergyAbsorber_post = 0. ;
  depositedEnergyWorld = 0. ;
  depositedEnergyFibresAtt->clear() ;
  depositedEnergyFibresAtt->resize(attLengths->size(),0);
   depositedEnergy_1st_Section = 0. ;
   depositedEnergy_2nd_Section = 0. ;
  depositedEnergyFibres_1st_Section = 0. ;
  depositedEnergyFibres_2nd_Section = 0. ;
  depositedEnergy_2nd_Sect_FibresCross = 0.;
  depositedEnergy_2nd_Sect_FibresCenter = 0.;
depositedEnergy_2nd_Sect_FibresCorners = 0.;

  depositedEnergyAbsorber_1st_Section = 0. ;
  depositedEnergyAbsorber_2nd_Section = 0.;

  totalTrackLengthFibres = 0.;
  totalTrackLengthOverThFibres = 0.;

  tot_phot_cer = 0;
  tot_phot_cer_post = 0;
  tot_det_phot_cer = 0;
  tot_gap_phot_cer = 0;
  tot_gap_photFast_cer->clear();
  tot_gap_photFast_cer->resize(attLengths->size(),0);
  tot_det_photFast_cer->clear();
  tot_det_photFast_cer->resize(attLengths->size(),0);

  for (int i = 0 ; i < 4 ; ++i)
  {
    inputMomentum->at (i) = 0. ;
  }
  for (int i = 0 ; i < 3 ; ++i)
  {
    inputInitialPosition->at (i) = 0. ;
  }

  depositedEnergies->clear () ;
  depositedEnergiesAtt->clear () ;
  depositedEnergiesAtt->resize(attLengths->size());
  depositFibres->clear () ;
 // depositFibres_1st_Section->clear () ;
 // depositFibres_2nd_Section->clear () ;
  scintillationPhotons->clear () ;
  scintillationFibres->clear () ;
  cerenkovPhotons->clear () ;
  cerenkovFibres->clear () ;

  Radial_stepLength = 0.;
  Longitudinal_stepLength = 0.;
  for(int i = 0; i < 5000; ++i)
  {
    Radial_ion_energy_absorber[i] = 0.;
    Longitudinal_ion_energy_absorber[i] = 0.;
  }

  for(int i = 0; i < 1000; ++i)
  {
    PrimaryParticleX[i] = 0.;
    PrimaryParticleY[i] = 0.;
    PrimaryParticleZ[i] = 0.;
    PrimaryParticleE[i] = 0.;
  }

  front_back    = 0 ;
  pmt_number    = 0 ;
  module_number    = 0 ;
  vertX         = 0 ;
  vertY         = 0 ;
  vertZ         = 0 ;
  PositionX     = 0 ;
  PositionY     = 0 ;
  PositionZ     = 0 ;
  PreMomentumX  = 0 ;
  PreMomentumY  = 0 ;
  PreMomentumZ  = 0 ;
  PostMomentumX = 0 ;
  PostMomentumY = 0 ;
  PostMomentumZ = 0 ;
  globalTime    = 0 ;
  PhotonEnergy  = 0 ;

  // abs_eventNumber  = 0;
  abs_x            = 0;
  abs_y            = 0;
  abs_z            = 0;
  abs_globalTime   = 0;
  abs_PhotonEnergy = 0;


}
