void fast(TString fileName, int primaries = 10000)
{

  // front == 0
  // back == 1

  TFile *_file0 = TFile::Open(fileName);
  TTree *photons = (TTree*) _file0->Get("photons");
  // photons->Draw("PositionY:PositionX","front_back == 0 && event == 5","COLZ");


  int   event   ;
  int   front_back    ;
  int   pmt_number    ;
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

  TBranch *b_event   ;
  TBranch *b_front_back    ;
  TBranch *b_pmt_number    ;
  TBranch *b_phPerPMT      ;
  TBranch *b_vertX         ;
  TBranch *b_vertY         ;
  TBranch *b_vertZ         ;
  TBranch *b_PositionX     ;
  TBranch *b_PositionY     ;
  TBranch *b_PositionZ     ;
  TBranch *b_PreMomentumX  ;
  TBranch *b_PreMomentumY  ;
  TBranch *b_PreMomentumZ  ;
  TBranch *b_PostMomentumX ;
  TBranch *b_PostMomentumY ;
  TBranch *b_PostMomentumZ ;
  TBranch *b_globalTime    ;
  TBranch *b_PhotonEnergy  ;



  photons->SetBranchAddress("event"         , &event         , &b_event);
  photons->SetBranchAddress("front_back"    , &front_back    , &b_front_back);
  photons->SetBranchAddress("pmt_number"    , &pmt_number    , &b_pmt_number);
  photons->SetBranchAddress("vertX"         , &vertX         , &b_vertX);
  photons->SetBranchAddress("vertY"         , &vertY         , &b_vertY);
  photons->SetBranchAddress("vertZ"         , &vertZ         , &b_vertZ);
  photons->SetBranchAddress("PositionX"     , &PositionX     , &b_PositionX);
  photons->SetBranchAddress("PositionY"     , &PositionY     , &b_PositionY);
  photons->SetBranchAddress("PositionZ"     , &PositionZ     , &b_PositionZ);
  photons->SetBranchAddress("PreMomentumX"  , &PreMomentumX  , &b_PreMomentumX);
  photons->SetBranchAddress("PreMomentumY"  , &PreMomentumY  , &b_PreMomentumY);
  photons->SetBranchAddress("PreMomentumZ"  , &PreMomentumZ  , &b_PreMomentumZ);
  photons->SetBranchAddress("PostMomentumX" , &PostMomentumX , &b_PostMomentumX);
  photons->SetBranchAddress("PostMomentumY" , &PostMomentumY , &b_PostMomentumY);
  photons->SetBranchAddress("PostMomentumZ" , &PostMomentumZ , &b_PostMomentumZ);
  photons->SetBranchAddress("globalTime"    , &globalTime    , &b_globalTime);
  photons->SetBranchAddress("PhotonEnergy"  , &PhotonEnergy  , &b_PhotonEnergy);

  // std::vector<data_t> data;
  // for(int i = 0; i < points ; i++)
  // {
  //   data_t temp_data;
  //   temp_data.z = 0;
  //
  // }

  int phPerPMT[2][9];
  TH1F *phPMT = new TH1F("phPMT","phPMT",100,-100,0);
  TH1F *phOUT = new TH1F("phOUT","phOUT",100,-100,0);

  for(int i = 0; i < 2; i++)
  {
    for(int j = 0; j < 9; j++)
    {
       phPerPMT[i][j] = 0;
    }
  }

  std::vector<float> zPos;

  for(int i = 0 ; i < photons->GetEntries(); i++)
  {
    photons->GetEvent(i);
    if(front_back == 0)
    {
      if(pmt_number == 4)
      {
        phPMT->Fill(vertZ);
        bool zIsThere = false;
        for(unsigned int iz = 0; iz < zPos.size() ; iz++ )
        {
          if(zPos[iz] == vertZ)
          {
            zIsThere = true;
          }
        }
        if(!zIsThere)
        {
          zPos.push_back(vertZ);
        }
      }
      else
      {
        phOUT->Fill(vertZ);
      }
    }
    phPerPMT[front_back][pmt_number]++;
  }

  std::vector<float> fX;
  std::vector<float> fY;
  std::vector<float> efX;
  std::vector<float> efY;

  std::vector<float> eff_X;
  std::vector<float> eff_Y;
  std::vector<float> e_eff_X;
  std::vector<float> e_eff_Y;




  for(int i = 1; i < phPMT->GetNbinsX()+1; i++)
  {
    if(phPMT->GetBinContent(i) != 0)
    {

      fX.push_back(phOUT->GetBinCenter(i) + 100.05);
      eff_X.push_back(phOUT->GetBinCenter(i) + 100.05);

      float fraction = (phOUT->GetBinContent(i))/(phPMT->GetBinContent(i) + phOUT->GetBinContent(i));
      float eFraction = fraction*sqrt( pow( (phOUT->GetBinError(i)/phOUT->GetBinContent(i)) ,2) + pow( (phPMT->GetBinError(i)/phPMT->GetBinContent(i)) ,2)   );

      float ext_eff = ((float) phPMT->GetBinContent(i))/((float) primaries);
      float e_ext_eff = ((float) phPMT->GetBinError(i))/((float) primaries);

      fY.push_back(fraction);
      efX.push_back(0);
      efY.push_back(eFraction);

      eff_Y.push_back(ext_eff);
      e_eff_X.push_back(0);
      e_eff_Y.push_back(e_ext_eff);
    }
  }

  // float firstT[21];
  // // float *firstT = new float(fX.size());
  // for (unsigned int i = 0 ; i < zPos.size() ; i++)
  // {
  //   std::cout << zPos[i] << " ";
  //   firstT[i] = INFINITY;
  // }
  // std::cout << std::endl;

  // for (unsigned int i = 0 ; i < fX.size() ; i++)
  // {
  //   std::vector<float> temp_v;
  //   tStamp.push_back(temp_v);
  // }
  std::vector<float> tStamp;

  for(int i = 0 ; i < photons->GetEntries(); i++)
  {
    photons->GetEvent(i);
    for(unsigned int iz = 0; iz < zPos.size(); iz++)
    {
      if(vertZ == zPos[iz]) // get the origin z
      {
        float dt = fabs(zPos[iz] - zPos[0])/(3.0e5);
        tStamp.push_back(globalTime + dt);
        // if(globalTime < firstT[iz])
        // {
        //   firstT[iz] = globalTime;
        // }
        // tStamp[iz].push_back(globalTime);
      }
    }
  }

  std::sort(tStamp.begin(),tStamp.end());
  //
  // for (unsigned int i = 0 ; i < zPos.size() ; i++)
  // {
  //   std::cout << firstT[i] << " ";
  // }
  std::cout << tStamp[0] <<std::endl;

  // for (unsigned int i = 0 ; i < tStamp.size() ; i++)
  // {
  //   std::cout << tStamp[i][0] << " ";
  // }


  TCanvas *c_fr = new TCanvas("c_fr","c_fr",1200,800);
  c_fr->SetGrid();
  TGraphErrors *g_fr = new TGraphErrors(eff_X.size(),&eff_X[0],&eff_Y[0],&e_eff_X[0],&e_eff_Y[0]);
  g_fr->SetName("g_fr");
  g_fr->SetTitle("Light collection efficiency");

  g_fr->GetXaxis()->SetTitle("Position in Fibre [mm]");
  g_fr->GetYaxis()->SetTitle("Efficiency");
  g_fr->GetYaxis()->SetRangeUser(0,0.1);
  g_fr->GetXaxis()->SetRangeUser(-40,120);
  g_fr->SetMarkerColor(4);
  g_fr->SetMarkerStyle(21);
  g_fr->Draw("AP");

  TCanvas *c_g = new TCanvas("c_g","c_g",1200,800);
  c_g->SetGrid();
  TGraphErrors *g_Light_Loss = new TGraphErrors(fX.size(),&fX[0],&fY[0],&efX[0],&efY[0]);
  g_Light_Loss->SetName("g_Light_Loss");
  g_Light_Loss->SetTitle("Optical crosstalk [PMT is on the left]");

  g_Light_Loss->GetXaxis()->SetTitle("Position in Fibre [mm]");
  g_Light_Loss->GetYaxis()->SetTitle("Fraction of optical crosstalk");
  g_Light_Loss->GetYaxis()->SetRangeUser(0,1);
  g_Light_Loss->GetXaxis()->SetRangeUser(-40,120);
  g_Light_Loss->SetMarkerColor(4);
  g_Light_Loss->SetMarkerStyle(21);
  // g_Light_Loss->Draw("ALP");
  g_Light_Loss->Draw("AP");

  new TCanvas;
  phPMT->Draw();
  new TCanvas;
  phOUT->Draw();

  gStyle->SetOptStat(0);
  TCanvas *c_h2front = new TCanvas("c_h2front","c_h2front",800,800);
  TH2F *h2front = new TH2F("h2front","",50,-30,30,50,30,30);
  h2front->SetTitle("Photons on PMTs");
  h2front->GetXaxis()->SetTitle("x [mm]");
  h2front->GetYaxis()->SetTitle("y [mm]");
  photons->Draw("PositionY:PositionX >> h2front","front_back == 0 && event == 0","COLZ");


  TCanvas *c_h2middle = new TCanvas("c_h2middle","c_h2middle",800,800);
  TH2F *h2middle = new TH2F("h2middle","",50,-30,30,50,-30,30);
  h2middle->SetTitle("Photons on PMTs");
  h2middle->GetXaxis()->SetTitle("x [mm]");
  h2middle->GetYaxis()->SetTitle("y [mm]");
  photons->Draw("PositionY:PositionX >> h2middle","front_back == 0 && event == 10","COLZ");

  TCanvas *c_h2back = new TCanvas("c_h2back","c_h2back",800,800);
  TH2F *h2back = new TH2F("h2back","",50,-30,30,50,-30,30);
  h2back->SetTitle("Photons on PMTs");
  h2back->GetXaxis()->SetTitle("x [mm]");
  h2back->GetYaxis()->SetTitle("y [mm]");
  photons->Draw("PositionY:PositionX >> h2back","front_back == 0 && event == 20","COLZ");


  // new TCanvas;














}
