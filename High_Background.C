#include "Riostream.h"
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>
#include <math.h>

//#define theta1 21.;
Int_t use_histo_method = 1;
Int_t show_histos = 1;
Int_t runlist[6] = {3892, 3893, 3894, 4073, 4074, 4075};
//Int_t runlist[1] = {4075};

Double_t Normfac = 0.242967e10;
Int_t nevts_SIMC = 50000;

Double_t charge = 21.2708;     //Scale the Al background to the charge of the production runs.
Double_t thickness = 0.1979;   //Scale the Al background down to the thickness of 3He cell.
Double_t rc_dummy = 1./0.548506, rc_walls = 1./0.681787;  //Scale Al background by the RC for dummy and 3He cells.
Double_t xb_nbins = 200.;
Double_t xbmin = 0., xbmax = 4.;
Double_t xb_binwidth = 0.;
//Double_t xmin = 2.92, xmax = 3.15;
Double_t xmin = 2.95, xmax = 3.10;
Double_t xmin_no_elastics = 2.8, xmax_no_elastics = 3.15;
//Double_t fitmin = 2.8, fitmax = 3.2;
Double_t fitmin = 2.3, fitmax = 4.;
//Double_t ymin = -0.028, ymax = 0.028;
Double_t ymin = -0.03, ymax = 0.03;      //I think this should be changed to a L.tr.vz cut instead.
Double_t thmin = -0.04, thmax = 0.055;
Double_t phmin = -0.0, phmax = 0.03;
Double_t dpmin = -0.03, dpmax = 0.03;
Double_t ymin_SIMC = ymin*100., ymax_SIMC = ymax*100.;
Double_t thmin_SIMC = TMath::ATan(thmin), thmax_SIMC = TMath::ATan(thmax);
Double_t phmin_SIMC = TMath::ATan(phmin), phmax_SIMC = TMath::ATan(phmax);
Double_t dpmin_SIMC = dpmin*100., dpmax_SIMC = dpmax*100.;
Double_t density_cut = -0.01469; //Y target position where the density profile changes.

void High_Background() 
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  TH1D *h1 = new TH1D("h1","Xbj" , xb_nbins, -phmax, phmax);
  TH1D *h2 = new TH1D("h2","Xbj" , xb_nbins, -phmax, phmax);
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  c1->SetLogy();

  for(Int_t m=0;m<6;m++)
    {
      Int_t run = runlist[m];

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
      
      if(!T->GetEntries())
	{
	  cerr<< "No root file was found" << endl;
	  return;
	}
      //==finish adding splits rootfiles=====================//

      //Read in the files containing beam trip cuts.
      FILE *fp1 = fopen(Form("/home/skbarcus/Tritium/Analysis/He3/Luminosity/%d_u1.txt",run),"r");
      FILE *fp2 = fopen(Form("/home/skbarcus/Tritium/Analysis/He3/Luminosity/%d_d1.txt",run),"r");

      Int_t skip1 = 0;                          //Gives number of lines to skip at top of data file. 
      Int_t nlines1 = 0;                        //Counts number of lines in the data file. 
      Int_t ncols1;                             //Set how many columns of data we have in the data file.
      char* string1[1000];                          //Variable to read lines of the data file.
      Float_t rising_temp_u1=0.,falling_temp_u1=0.,rising_temp_evt_u1=0.,falling_temp_evt_u1=0.;
      Float_t rising_u1[10]={},falling_u1[10]={},rising_evt_u1[10]={},falling_evt_u1[10]={};

      Int_t skip2 = 0;                          //Gives number of lines to skip at top of data file. 
      Int_t nlines2 = 0;                        //Counts number of lines in the data file. 
      Int_t ncols2;                             //Set how many columns of data we have in the data file.
      char* string2[1000];                          //Variable to read lines of the data file.
      Float_t rising_temp_d1=0.,falling_temp_d1=0.,rising_temp_evt_d1=0.,falling_temp_evt_d1=0.;
      Float_t rising_d1[10]={},falling_d1[10]={},rising_evt_d1[10]={},falling_evt_d1[10]={};

      //Read in data.
      while (1) 
	{
	  //Skips the first skip lines of the file. 
	  if (nlines1 < skip1)
	    {
	      fgets(string1,1000,fp1);
	      nlines1++;
	    }
	//Reads the two columns of data into x and y.
	  else
	    {
	      //Read in the number of columns of data in your data file. 
	      ncols1 = fscanf(fp1,"%f %f %f %f",&rising_temp_evt_u1,&rising_temp_u1,&falling_temp_evt_u1,&falling_temp_u1);
	      if (ncols1 < 0) break;  
	      rising_evt_u1[nlines1-skip1] = rising_temp_evt_u1;
	      rising_u1[nlines1-skip1] = rising_temp_u1;
	      falling_evt_u1[nlines1-skip1] = falling_temp_evt_u1;
	      falling_u1[nlines1-skip1] = falling_temp_u1; 
	      nlines1++;
	    }
	}
      
      //Read in data.
      while (1) 
	{
	  //Skips the first skip lines of the file. 
	  if (nlines2 < skip2)
	    {
	      fgets(string2,1000,fp2);
	      nlines2++;
	    }
	//Reads the two columns of data into x and y.
	  else
	    {
	      //Read in the number of columns of data in your data file. 
	      ncols2 = fscanf(fp2,"%f %f %f %f",&rising_temp_evt_d1,&rising_temp_d1,&falling_temp_evt_d1,&falling_temp_d1);
	      if (ncols2 < 0) break;   
	      rising_evt_d1[nlines2-skip2] = rising_temp_evt_d1;
	      rising_d1[nlines2-skip2] = rising_temp_d1;
	      falling_evt_d1[nlines2-skip2] = falling_temp_evt_d1;
	      falling_d1[nlines2-skip2] = falling_temp_d1; 
	      nlines2++;
	    }
	}

      //Print rising and falling edges if beam trip cut for u1.
      for(Int_t i=0;i<nlines1;i++)
	{
	  cout<<"rising_evt_u1["<<i<<"] = "<<rising_evt_u1[i]<<"   rising_u1["<<i<<"] = "<<rising_u1[i]<<"   falling_evt_u1["<<i<<"] = "<<falling_evt_u1[i]<<"   falling_u1["<<i<<"] = "<<falling_u1[i]<<endl;
	}
      cout<<"********************************************************"<<endl;

      //Print rising and falling edges if beam trip cut for d1.
      for(Int_t i=0;i<nlines2;i++)
	{
	  cout<<"rising_evt_d1["<<i<<"] = "<<rising_evt_d1[i]<<"   rising_d1["<<i<<"] = "<<rising_d1[i]<<"   falling_evt_d1["<<i<<"] = "<<falling_evt_d1[i]<<"   falling_d1["<<i<<"] = "<<falling_d1[i]<<endl;
	}
      
      T->SetBranchStatus("*",0);
      T->SetBranchStatus("EKL.x_bj",1);
      T->SetBranchStatus("ExtEKL.x_bj",1);
      T->SetBranchStatus("L.tr.tg_y",1);
      T->SetBranchStatus("L.tr.n",1);
      T->SetBranchStatus("L.cer.asum_c",1);
      T->SetBranchStatus("L.prl1.e",1);
      T->SetBranchStatus("L.prl2.e",1);
      T->SetBranchStatus("L.tr.tg_th",1);
      T->SetBranchStatus("L.tr.tg_ph",1);
      T->SetBranchStatus("L.tr.tg_dp",1);
      //T->SetBranchStatus("L.tr.p",1);
      T->SetBranchStatus("DBB.evtypebits",1);
      T->SetBranchStatus("right_bcm_u1c",1);
      T->SetBranchStatus("right_bcm_d1c",1);
      T->SetBranchStatus("right_clkcount",1);
      
      Double_t x_bj=0.,Ext_x_bj=0.,L_tr_tg_y[21],L_tg_th[21],L_tg_ph[21],L_tg_dp[21],L_cer=0.,L_prl1_e=0.,L_prl2_e=0.,L_tr_p=0.,L_tr_n=0.,evtypebits;
      Double_t u1c=0.,d1c=0.,rclk=0.;
      
      T->SetBranchAddress("EKL.x_bj",&x_bj);
      T->SetBranchAddress("ExtEKL.x_bj",&Ext_x_bj);
      T->SetBranchAddress("L.tr.tg_y",&L_tr_tg_y);
      T->SetBranchAddress("L.tr.n",&L_tr_n);
      T->SetBranchAddress("L.cer.asum_c",&L_cer);
      T->SetBranchAddress("L.prl1.e",&L_prl1_e);
      T->SetBranchAddress("L.prl2.e",&L_prl2_e);
      T->SetBranchAddress("L.tr.tg_th",L_tg_th);
      T->SetBranchAddress("L.tr.tg_ph",L_tg_ph);
      T->SetBranchAddress("L.tr.tg_dp",L_tg_dp);
      //T->SetBranchAddress("L.tr.p",L_tr_p);
      T->SetBranchAddress("DBB.evtypebits",&evtypebits);
      T->SetBranchAddress("right_bcm_u1c",&u1c);
      T->SetBranchAddress("right_bcm_d1c",&d1c);
      T->SetBranchAddress("right_clkcount",&rclk);
      
      Int_t nevts = T->GetEntries();
      Int_t nevts_cuts = 0;   //Count number of events surviving all of the cuts.
      //Double_t xmin = 0., xmax = 0.;
      
      cout<<"nevts = "<<nevts<<endl;

      //Create and fill histogram with Xbj values with various cuts applied to the data.
      //TH1D *h1 = new TH1D("h1","Xbj" , xb_nbins, xbmin, xbmax);
      Int_t evt_136=0,evt_8=0;
      Int_t maxtracks=0;
      for(Int_t i=0;i<nevts;i++) 
	{
	  T->GetEntry(i);
	  //Usual cuts without Xbj cut.
	  if(L_tr_tg_y[0]>ymin && L_tr_tg_y[0]<ymax && L_tr_n==1 && (evtypebits&1<<3)==1<<3 && L_prl1_e>(-L_prl2_e+2000) && L_tg_ph[0]>phmin && L_tg_ph[0]<phmax && L_tg_th[0]>thmin && L_tg_th[0]<thmax && L_tg_dp[0]>dpmin && L_tg_dp[0]<dpmax)
	    {
	      //Beam trip cuts. Cut out BCM scaler values in the range of the beam trips.
	      if((u1c>=rising_u1[0]&&u1c<=falling_u1[0])||(u1c>=rising_u1[1]&&u1c<=falling_u1[1])||(u1c>=rising_u1[2]&&u1c<=falling_u1[2])||(u1c>=rising_u1[3]&&u1c<=falling_u1[3])||(u1c>=rising_u1[4]&&u1c<=falling_u1[4])||(u1c>=rising_u1[5]&&u1c<=falling_u1[5])||(u1c>=rising_u1[6]&&u1c<=falling_u1[6])||(u1c>=rising_u1[7]&&u1c<=falling_u1[7])||(u1c>=rising_u1[8]&&u1c<=falling_u1[8])||(u1c>=rising_u1[9]&&u1c<=falling_u1[9]))
		{ 
		  h1->Fill(L_tg_ph[0]);
		  nevts_cuts++;
		}
	    }
	  //Usual cuts with added Xbj cut to isolate tail after elastic peak.
	  if(L_tr_tg_y[0]>ymin && L_tr_tg_y[0]<ymax && L_tr_n==1 && (evtypebits&1<<3)==1<<3 && L_prl1_e>(-L_prl2_e+2000) && L_tg_ph[0]>phmin && L_tg_ph[0]<phmax && L_tg_th[0]>thmin && L_tg_th[0]<thmax && L_tg_dp[0]>dpmin && L_tg_dp[0]<dpmax&&Ext_x_bj>3.1)
	    {
	      //Beam trip cuts. Cut out BCM scaler values in the range of the beam trips.
	      if((u1c>=rising_u1[0]&&u1c<=falling_u1[0])||(u1c>=rising_u1[1]&&u1c<=falling_u1[1])||(u1c>=rising_u1[2]&&u1c<=falling_u1[2])||(u1c>=rising_u1[3]&&u1c<=falling_u1[3])||(u1c>=rising_u1[4]&&u1c<=falling_u1[4])||(u1c>=rising_u1[5]&&u1c<=falling_u1[5])||(u1c>=rising_u1[6]&&u1c<=falling_u1[6])||(u1c>=rising_u1[7]&&u1c<=falling_u1[7])||(u1c>=rising_u1[8]&&u1c<=falling_u1[8])||(u1c>=rising_u1[9]&&u1c<=falling_u1[9]))
		{ 
		  h2->Fill(L_tg_ph[0]);
		  cout<<"phi = "<<L_tg_ph[0]<<"   dp = "<<L_tg_dp[0]<<"   th = "<<L_tg_th[0]<<endl;
		}
	    }

	}
      cout<<"Number of events surviving cuts = "<<nevts_cuts<<endl;
    }//End loop over runlist. 

  gStyle->SetOptStat(0);
  h1->Draw();
  h2->Draw("same");
  h2->SetLineColor(2);
  //cout<<"# evt_136 = "<<evt_136<<endl;
  //cout<<"# evt_8 = "<<evt_8<<endl;

  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
