void dp_vs_th()
{  
  Int_t runlist[5] = {3892, 3893, 3894, 4074, 4075}; 
  Double_t shmin=        0.2;
  Double_t ymin=         -0.03,        ymax=         0.03;
  Double_t thmin=        -0.03,        thmax=        0.03;
  Double_t phmin=        -0.02,        phmax=        0.02;
  Double_t dpmin=        -0.0,         dpmax=        0.03;
  Double_t xbmin=        2.96,         xbmax=        3.06;
  Double_t ymin_SIMC=    ymin*100.,    ymax_SIMC=    ymax*100.;
  Double_t thmin_SIMC=   thmin,        thmax_SIMC=   thmax;
  Double_t phmin_SIMC=   phmin,        phmax_SIMC=   phmax;
  Double_t dpmin_SIMC=   dpmin*100.,   dpmax_SIMC=   dpmax*100.;
  Double_t xbmin_SIMC=   xbmin,        xbmax_SIMC=   xbmax;

  Int_t    nbins=        50;
  Int_t    nevts_SIMC=   50000;
  Double_t normfac=      0.242967e10;
  Double_t scale1=       0.3,         scale2=   0.3;

  TChain *T = new TChain("T");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3892.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3893.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3894.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4073.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4074*.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4075*.root");
      
  //Read in SIMC data.
  TChain *SIMC = new TChain("h666");
  SIMC->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_y20_dp3_x30.root");

  T->SetBranchStatus("*",                0);
  T->SetBranchStatus("EKL.x_bj",         1);
  T->SetBranchStatus("L.tr.tg_y",        1);
  T->SetBranchStatus("L.tr.n",           1);
  T->SetBranchStatus("L.cer.asum_c",     1);
  T->SetBranchStatus("L.prl1.e",         1);
  T->SetBranchStatus("L.prl2.e",         1);
  T->SetBranchStatus("L.tr.p",           1);
  T->SetBranchStatus("DBB.evtypebits",   1);
  T->SetBranchStatus("L.cer.trx",        1);
  T->SetBranchStatus("L.cer.try",        1);
  T->SetBranchStatus("L.tr.tg_th",       1);
  T->SetBranchStatus("L.tr.tg_ph",       1);
  T->SetBranchStatus("L.tr.tg_dp",       1);
  T->SetBranchStatus("L.tr.p",           1);
      
  Double_t x_bj=0.,L_tr_tg_y[21],L_cer=0.,L_prl1_e=0.,L_prl2_e=0.,L_tr_p[21],L_tr_n=0.,evtypebits;
      
  T->SetBranchAddress("EKL.x_bj",         &x_bj);
  T->SetBranchAddress("L.tr.tg_y",        L_tr_tg_y);
  T->SetBranchAddress("L.tr.n",           &L_tr_n);
  T->SetBranchAddress("L.cer.asum_c",     &L_cer);
  T->SetBranchAddress("L.prl1.e",         &L_prl1_e);
  T->SetBranchAddress("L.prl2.e",         &L_prl2_e);
  T->SetBranchAddress("L.tr.p",           L_tr_p);
  T->SetBranchAddress("DBB.evtypebits",   &evtypebits);
      
  Int_t nevts = T->GetEntries();
  cout<<"nevts = "<<nevts<<endl;
  //TH1D *h1 = new TH1D("h1", "Xbj", 1000, 0., 4.);
  //TH1D *h2 = new TH1D("h2", "Xbj", 1000, 0., 4.);

  SIMC->SetBranchStatus("*",         0);
  SIMC->SetBranchStatus("xbj",       1);
  SIMC->SetBranchStatus("ssytar",    1);
  SIMC->SetBranchStatus("ssyptar",   1);
  SIMC->SetBranchStatus("ssxptar",   1);
  SIMC->SetBranchStatus("ssdelta",   1);
  SIMC->SetBranchStatus("Weight",    1);

  nevts_SIMC = SIMC->GetEntries();
  cout<<"nevts_SIMC = "<<nevts_SIMC<<endl;

  //Define some cuts.
  TCut ct_ep=        Form("(L.prl1.e+L.prl2.e)/(1000*L.tr.p)>%f",     shmin);
  TCut ct_1tr=       Form("L.tr.n==1"                                      );
  TCut ct_th=        Form("L.tr.tg_th>%f&&L.tr.tg_th<%f",       thmin,thmax);
  TCut ct_ph=        Form("L.tr.tg_ph>%f&&L.tr.tg_ph<%f",       phmin,phmax);
  TCut ct_dp=        Form("L.tr.tg_dp>%f&&L.tr.tg_dp<%f",       dpmin,dpmax);
  TCut ct_T3=        Form("(DBB.evtypebits&1<<3)==1<<3"                    );
  TCut ct_T7=        Form("(DBB.evtypebits&1<<7)==1<<7"                    );
  TCut ct_xb=        Form("EKL.x_bj>%f&&EKL.x_bj<%f",           xbmin,xbmax);
  TCut ct_y=         Form("L.tr.tg_y>%f&&L.tr.tg_y<%f",           ymin,ymax);

  TCut ct_y_SIMC=    Form("ssytar>%f&&ssytar<%f",       ymin_SIMC,ymax_SIMC);
  TCut ct_th_SIMC=   Form("ssxptar>%f&&ssxptar<%f",   thmin_SIMC,thmax_SIMC);
  TCut ct_ph_SIMC=   Form("ssyptar>%f&&ssyptar<%f",   phmin_SIMC,phmax_SIMC);
  TCut ct_dp_SIMC=   Form("ssdelta>%f&&ssdelta<%f",   dpmin_SIMC,dpmax_SIMC);
  TCut ct_xb_SIMC=   Form("xbj>%f&&xbj<%f",           xbmin_SIMC,xbmax_SIMC);

  //Draw experimental elastic band from dp vs. theta (scattering angle).
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  //T->Draw("L.cer.trx:L.cer.try>>h1(200,-1.,1.,200,-1.,1.)",ct_T3&&ct_1tr&&ct_th&&ct_ph&&ct_dp&&ct_ep,"colz");
  T->Draw(Form("L.tr.tg_dp:L.tr.tg_ph>>h1(%d,%f,%f,%d,%f,%f)",nbins,phmin,phmax,nbins,dpmin,dpmax),ct_y&&ct_1tr&&ct_T3&&ct_ep&&ct_th&&ct_ph&&ct_dp,"colz");

  //Draw SIMC elastic band dp vs. theta.
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();
  h666->Draw(Form("ssdelta:ssyptar>>h2(%d,%f,%f,%d,%f,%f)",nbins,phmin,phmax,nbins,dpmin*100.,dpmax*100.),Form("Weight*%f/%d",normfac,nevts_SIMC),"colz");

  //Compare experimental dp to SIMC dp.
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();

  TH1D *hscaled_dp = new TH1D("hscaled_dp","Scaled SIMC dp and experimental dp" ,nbins,dpmin,dpmax);
  //h666->Draw(Form("ssdelta/100.>>h4(%d,%f,%f)",nbins,dpmin,dpmax),Form("Weight*%f/%d",normfac,nevts_SIMC),"colz");
  T->Draw(Form("L.tr.tg_dp>>h3(%d,%f,%f)",nbins,dpmin,dpmax),ct_y&&ct_1tr&&ct_T3&&ct_ep&&ct_th&&ct_ph&&ct_dp&&ct_xb);
  h666->Draw(Form("ssdelta/100.>>h4(%d,%f,%f)",nbins,dpmin,dpmax),Form("Weight*%f/%d",normfac,nevts_SIMC)*ct_y_SIMC*ct_ph_SIMC*ct_th_SIMC*ct_dp_SIMC,"same");
  hscaled_dp->Add(h4,scale1);
  hscaled_dp->Draw("same");
  hscaled_dp->SetLineColor(2);
  h4->SetLineColor(0);

  //Compare experimental scattering angle to SIMC scattering angle.
  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();

  TH1D *hscaled_ph = new TH1D("hscaled_ph","Scaled SIMC ph and experimental ph" ,nbins,phmin,phmax);
  //h666->Draw(Form("ssdelta/100.>>h4(%d,%f,%f)",nbins,dpmin,dpmax),Form("Weight*%f/%d",normfac,nevts_SIMC),"colz");
  T->Draw(Form("L.tr.tg_ph>>h5(%d,%f,%f)",nbins,phmin,phmax),ct_y&&ct_1tr&&ct_T3&&ct_ep&&ct_th&&ct_ph&&ct_dp&&ct_xb);
  h666->Draw(Form("ssyptar>>h6(%d,%f,%f)",nbins,phmin,phmax),Form("Weight*%f/%d",normfac,nevts_SIMC)*ct_y_SIMC*ct_ph_SIMC*ct_th_SIMC*ct_dp_SIMC,"same");
  hscaled_ph->Add(h6,scale2);
  hscaled_ph->Draw("same");
  hscaled_ph->SetLineColor(2);
  h6->SetLineColor(0);
}
