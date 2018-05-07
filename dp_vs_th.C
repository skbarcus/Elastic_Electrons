void dp_vs_th()
{  
  Int_t runlist[5] = {3892, 3893, 3894, 4074, 4075}; 
  Double_t shmin=0.65;
  Double_t thmin=-0.07, thmax=0.07;
  Double_t phmin=-0.03, phmax=0.03;
  Double_t dpmin=-0.0, dpmax=0.03;
  Double_t xbmin=2.96, xbmax=3.06;
  Int_t nbins =50;
  Double_t normfac=0.242967e10;
  Int_t nevts_SIMC=50000;
  
      
  for(Int_t i=0;i<5;i++)
    {
      Int_t run=runlist[i];
      //============  Reading the Rootfile =======================//
      
      const TString rootfilePath = "/home/skbarcus/Tritium/Analysis/He3/Rootfiles/";
      std::ostringstream str;
      str << rootfilePath<<"e08014_"<<run;
      TString basename = str.str().c_str();
      TString rootfile = basename + ".root";
      TChain* T;
      TChain*ev;
      
      if(run < 20000)
	{
	  T = new TChain("T");
	  ev = new TChain("RIGHT");
	}
      if(run > 20000)
	{
	  T = new TChain("T");
	  ev = new TChain("EVRIGHT");
	}
      
      Long_t split=0;
      
      char* file = 0;
      
      //====adding splits rootfiles =======================//
      
      Long_t u=0;
      while ( !gSystem->AccessPathName(rootfile.Data()) ) 
	{
	  T->Add(rootfile.Data());
	  ev->Add(rootfile.Data());
	  cout << "ROOT file " << rootfile << " added to TChain." << endl;
	  u++;
	  rootfile = basename + "_" + u + ".root";
	}
    }
      if(!T->GetEntries())
	{
	  cerr<< "Not root file was found" << endl;
	  return;
	}
      //==finish adding splits rootfiles=====================//
      
      T->SetBranchStatus("*",0);
      T->SetBranchStatus("EKL.x_bj",1);
      T->SetBranchStatus("L.tr.tg_y",1);
      T->SetBranchStatus("L.tr.n",1);
      T->SetBranchStatus("L.cer.asum_c",1);
      T->SetBranchStatus("L.prl1.e",1);
      T->SetBranchStatus("L.prl2.e",1);
      //T->SetBranchStatus("L.tr.p",1);
      T->SetBranchStatus("DBB.evtypebits",1);
      T->SetBranchStatus("L.cer.trx",1);
      T->SetBranchStatus("L.cer.try",1);
      T->SetBranchStatus("L.tr.tg_th",1);
      T->SetBranchStatus("L.tr.tg_ph",1);
      T->SetBranchStatus("L.tr.tg_dp",1);
      T->SetBranchStatus("L.tr.p",1);
      
      Double_t x_bj=0.,L_tr_tg_y[21],L_cer=0.,L_prl1_e=0.,L_prl2_e=0.,L_tr_p[21],L_tr_n=0.,evtypebits;
      
      T->SetBranchAddress("EKL.x_bj",&x_bj);
      T->SetBranchAddress("L.tr.tg_y",L_tr_tg_y);
      T->SetBranchAddress("L.tr.n",&L_tr_n);
      T->SetBranchAddress("L.cer.asum_c",&L_cer);
      T->SetBranchAddress("L.prl1.e",&L_prl1_e);
      T->SetBranchAddress("L.prl2.e",&L_prl2_e);
      //T->SetBranchAddress("L.tr.p",L_tr_p);
      T->SetBranchAddress("DBB.evtypebits",&evtypebits);
      
      Int_t nevts = T->GetEntries();
      cout<<"nevts = "<<nevts<<endl;
      //TH1D *h1 = new TH1D("h1","Xbj" , 1000, 0., 4.);
      //TH1D *h2 = new TH1D("h2","Xbj" , 1000, 0., 4.);

      //Define some cuts.
      TCut ct_ep = Form("(L.prl1.e+L.prl2.e)/(1000*L.tr.p)>%f",shmin);
      TCut ct_1tr = Form("L.tr.n==1");
      TCut ct_th = Form("L.tr.tg_th>%f&&L.tr.tg_th<%f",thmin,thmax);
      TCut ct_ph = Form("L.tr.tg_ph>%f&&L.tr.tg_ph<%f",phmin,phmax);
      TCut ct_dp = Form("L.tr.tg_dp>%f&&L.tr.tg_dp<%f",dpmin,dpmax);
      TCut ct_T3 = Form("(DBB.evtypebits&1<<3)==1<<3");
      TCut ct_T7 = Form("(DBB.evtypebits&1<<7)==1<<7");
      TCut ct_xb = Form("EKL.x_bj>%f&&EKL.x_bj<%f",xbmin,xbmax);

      //Draw experimental elastic band from dp vs. theta (scattering angle).
      TCanvas* c1=new TCanvas("c1");
      c1->SetGrid();
      //T->Draw("L.cer.trx:L.cer.try>>h1(200,-1.,1.,200,-1.,1.)",ct_T3&&ct_1tr&&ct_th&&ct_ph&&ct_dp&&ct_ep,"colz");
      T->Draw(Form("L.tr.tg_dp:L.tr.tg_ph>>h1(%d,%f,%f,%d,%f,%f)",nbins,phmin,phmax,nbins,dpmin,dpmax),ct_1tr&&ct_T3&&ct_ep&&ct_th&&ct_ph&&ct_dp&&ct_xb,"colz");

      //Read in SIMC data.
      TChain *SIMC = new TChain("h666");
      SIMC->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic.root");
      SIMC->SetBranchStatus("*",0);
      SIMC->SetBranchStatus("xbj",1);
      SIMC->SetBranchStatus("ssytar",1);
      SIMC->SetBranchStatus("ssyptar",1);
      SIMC->SetBranchStatus("ssxptar",1);
      SIMC->SetBranchStatus("ssdelta",1);
      SIMC->SetBranchStatus("Weight",1);

      nevts_SIMC = SIMC->GetEntries();
      cout<<"nevts_SIMC = "<<nevts_SIMC<<endl;

      //Draw SIMC elastic band dp vs. theta.
      TCanvas* c2=new TCanvas("c2");
      c2->SetGrid();
      h666->Draw(Form("ssdelta:ssyptar>>h2(%d,%f,%f,%d,%f,%f)",nbins,phmin,phmax,nbins,dpmin*100.,dpmax*100.),Form("Weight*%f/%d",normfac,nevts_SIMC),"colz");

      //Compare experimental dp to SIMC dp.
      TCanvas* c3=new TCanvas("c3");
      c3->SetGrid();

      TH1D *hscaled_dp = new TH1D("hscaled_dp","Scaled SIMC dp and experimental dp" ,nbins,dpmin,dpmax);
      //h666->Draw(Form("ssdelta/100.>>h4(%d,%f,%f)",nbins,dpmin,dpmax),Form("Weight*%f/%d",normfac,nevts_SIMC),"colz");
      T->Draw(Form("L.tr.tg_dp>>h3(%d,%f,%f)",nbins,dpmin,dpmax),ct_1tr&&ct_T3&&ct_ep&&ct_th&&ct_ph&&ct_dp&&ct_xb,"colz");
      h666->Draw(Form("ssdelta/100.>>h4(%d,%f,%f)",nbins,dpmin,dpmax),Form("Weight*%f/%d",normfac,nevts_SIMC),"same colz");
      hscaled_dp->Add(h4,0.75);
      hscaled_dp->Draw("same");
      hscaled_dp->SetLineColor(2);
      h4->SetLineColor(0);

      //Compare experimental scattering angle to SIMC scattering angle.
      TCanvas* c4=new TCanvas("c4");
      c4->SetGrid();

      TH1D *hscaled_ph = new TH1D("hscaled_ph","Scaled SIMC ph and experimental ph" ,nbins,phmin,phmax);
      //h666->Draw(Form("ssdelta/100.>>h4(%d,%f,%f)",nbins,dpmin,dpmax),Form("Weight*%f/%d",normfac,nevts_SIMC),"colz");
      T->Draw(Form("L.tr.tg_ph>>h5(%d,%f,%f)",nbins,phmin,phmax),ct_1tr&&ct_T3&&ct_ep&&ct_th&&ct_ph&&ct_dp&&ct_xb,"colz");
      h666->Draw(Form("ssyptar>>h6(%d,%f,%f)",nbins,phmin,phmax),Form("Weight*%f/%d",normfac,nevts_SIMC),"same colz");
      hscaled_ph->Add(h6,0.50);
      hscaled_ph->Draw("same");
      hscaled_ph->SetLineColor(2);
      h6->SetLineColor(0);
}
