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
Int_t use_split_density = 0;
Int_t match_data = 1;
Int_t gaus_or_total = 2;        //0-> match height of the data and SIMC elastic peaks using the max value of the Gaussian in the total fits. 1-> match height of the data and SIMC elasstic peaks using the max value (in the elastic peak region) of the total fit. 2-> match the areas of the Gaussian parts of the total data and total SIMC fits. 
Int_t subtract_SIMC = 0;             //0-> The exponential fit of the experimental background is used to fill the h_no_elastics histo. It uses the region of fitmin to fitmax. 1-> then the exponential background fit of the experimental data has the generated SIMC elastic events subtracted from it before the h_no_elastics histo is filled. 
Int_t use_lorentz = 0;
Double_t scale_SIMC = 1.;       //Scale factor for SIMC histogram.
Int_t runlist[6] = {3892, 3893, 3894, 4073, 4074, 4075};
//Int_t runlist[1] = {4075};
  
//Double_t Normfac1 = 0.278306e10, Normfac2 = 0.188633e10;
//Double_t Normfac1 = 0.272783e10, Normfac2 = 0.183444e10;
//Double_t Normfac1 = 0.282782e10, Normfac2 = 0.191004e10;  //XS * 0.59552
//Double_t Normfac1 = 0.856831e9;      //XS * 0.59552 0.01045 g/cm^3
//Double_t Normfac1 = 0.112683e10;      //XS * 0.59552 0.01375 g/cm^3
//Double_t Normfac1 = 0.856831e9;      //XS * 1.41330 0.01375 g/cm^3
//Double_t Normfac1 = 0.112683e10;      //XS * 1.08312 0.01375 g/cm^3
Double_t Normfac1 = 0.106612e10;      //XS * 1.08312 0.013 g/cm^3 or XS * 1.1881 or XS * 0.99206.
Int_t nevts_SIMC = 100000;
  
Double_t charge = 21.2708;     //Scale the Al background to the charge of the production runs.
Double_t thickness = 0.1979;   //Scale the Al background down to the thickness of 3He cell.
Double_t rc_dummy = 1./0.548506, rc_walls = 1./0.681787;  //Scale Al background by the RC for dummy and 3He cells.
Double_t xb_nbins = 200.;
Double_t xbmin = 0., xbmax = 4.;
Double_t xb_binwidth = fabs(xbmax-xbmin)/xb_nbins;;
//Double_t xmin = 2.92, xmax = 3.15;
Double_t xmin = 2.95, xmax = 3.10;
Double_t xmin_no_elastics = 2.8, xmax_no_elastics = 3.15;
//Double_t fitmin = 2.8, fitmax = 3.2;
Double_t fitmin = 2.3, fitmax = 3.25;    //2.3, 3.15
//Double_t ymin = -0.028, ymax = 0.028;
Double_t ymin = -0.03, ymax = 0.03;      //I think this should be changed to a L.tr.vz cut instead.
Double_t thmin = -0.035, thmax = 0.035;  //Dien's cuts.
Double_t phmin = -0.025, phmax = 0.025;
Double_t dpmin = -0.035, dpmax = 0.035;
//Double_t thmin = -0.042, thmax = 0.049;
//Double_t phmin = -0.03, phmax = 0.03;
//Double_t dpmin = -0.02, dpmax = 0.03;
Double_t ymin_SIMC = ymin*100., ymax_SIMC = ymax*100.;
Double_t thmin_SIMC = TMath::ATan(thmin), thmax_SIMC = TMath::ATan(thmax);
Double_t phmin_SIMC = TMath::ATan(phmin), phmax_SIMC = TMath::ATan(phmax);
Double_t dpmin_SIMC = dpmin*100., dpmax_SIMC = dpmax*100.;
Double_t density_cut = -0.01469; //Y target position where the density profile changes.
Double_t pr_yint = 2000.;
Double_t GC_eff = 203./196.;     //GC efficiency correction to the elastic electron peak height. 
Double_t DT_correction = 1./0.9527;

void Xbj_New_Fit() 
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);
  
  TH1D *h1 = new TH1D("h1","Xbj" , xb_nbins, xbmin, xbmax);
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  c1->SetLogy();

  for(Int_t m=0;m<6;m++)
    {
      run = runlist[m];

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
      
      /*
	TChain *T = new TChain("T");
	T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3892.root");
	T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3893.root");
	T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3894.root");
	T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4073.root");
	T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4074*.root");
	T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4075*.root");
      */
      
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
      
      //Print u1c entries.
      /*
	for(Int_t i=0;i<100;i++)
	{
	T->GetEntry(i);
	cout<<"u1c = "<<u1c<<endl;
	}
      */
      
      //T->Draw("EKL.x_bj>>h2(1000,0.5,3.5)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1");
      
      //T->Draw("EKL.x_bj>>h2(1000,0.5,3.5)","L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1");
      
      //Create and fill histogram with Xbj values with various cuts applied to the data.
      //TH1D *h1 = new TH1D("h1","Xbj" , xb_nbins, xbmin, xbmax);
      Int_t evt_136=0,evt_8=0;
      Int_t maxtracks=0;
      for(Int_t i=0;i<nevts;i++) 
	{
	  T->GetEntry(i);
	  //T->Show(i);
	  /*if( (evtypebits&1<<3)==1<<3 )
	    {
	    cout<<evtypebits<<endl;
	    }*/
	  /*if(L_tr_n>maxtracks)
	    {
	    maxtracks=L_tr_n;
	    cout<<"Max Track Number = "<<maxtracks<<" at event "<<i<<"."<<endl;
	    }*/
	  /*
	    if(evtypebits==136)
	    {
	    evt_136++;
	    }
	    if(evtypebits==8)
	    {
	    evt_8++;
	    }*/
	  //if(L_tr_tg_y[0]>-0.028 && L_tr_tg_y[0]<0.028 && L_tr_n==1 && (evtypebits&1<<3)==1<<3)
	  //if(L_tr_tg_y[0]>-0.028 && L_tr_tg_y[0]<0.028 && L_tr_n==1 && (evtypebits&1<<3)==1<<3 && L_cer>150) //Add GC cut.
	  //if(L_tr_tg_y[0]>ymin && L_tr_tg_y[0]<ymax && L_tr_n==1 && (evtypebits&1<<3)==1<<3 && L_prl1_e>(-L_prl2_e+2000) && L_tg_th[0]>thmin && L_tg_th[0]<thmax) 
	  //if(L_tr_tg_y[0]>ymin && L_tr_tg_y[0]<ymax && L_tr_n==1 && (evtypebits&1<<3)==1<<3 && L_prl1_e>(-L_prl2_e+2000) && L_tg_ph[0]>phmin && L_tg_ph[0]<phmax)
	  //Various cuts.
	  if(L_tr_tg_y[0]>ymin && L_tr_tg_y[0]<ymax && L_tr_n==1 && (evtypebits&1<<3)==1<<3 && L_prl1_e>(-L_prl2_e+pr_yint) && L_tg_ph[0]>phmin && L_tg_ph[0]<phmax && L_tg_th[0]>thmin && L_tg_th[0]<thmax && L_tg_dp[0]>dpmin && L_tg_dp[0]<dpmax)
	    {
	      //Beam trip cuts. Cut out BCM scaler values in the range of the beam trips.
	      if((u1c>=rising_u1[0]&&u1c<=falling_u1[0])||(u1c>=rising_u1[1]&&u1c<=falling_u1[1])||(u1c>=rising_u1[2]&&u1c<=falling_u1[2])||(u1c>=rising_u1[3]&&u1c<=falling_u1[3])||(u1c>=rising_u1[4]&&u1c<=falling_u1[4])||(u1c>=rising_u1[5]&&u1c<=falling_u1[5])||(u1c>=rising_u1[6]&&u1c<=falling_u1[6])||(u1c>=rising_u1[7]&&u1c<=falling_u1[7])||(u1c>=rising_u1[8]&&u1c<=falling_u1[8])||(u1c>=rising_u1[9]&&u1c<=falling_u1[9]))
		{ 
		  h1->Fill(Ext_x_bj);
		  nevts_cuts++;
		}
	    }
	}
      cout<<"Number of events surviving cuts = "<<nevts_cuts<<endl;
    }//End loop over runlist. 

  gStyle->SetOptStat(0);
  h1->Draw();
  //cout<<"# evt_136 = "<<evt_136<<endl;
  //cout<<"# evt_8 = "<<evt_8<<endl;

  //Create an exponential fit function to fit the region with the elastic peak, but excluding the peak itself.
  Double_t fit_exp(Double_t *x,Double_t *par) 
  {
    Bool_t reject;
    //reject = kTRUE;
    reject = kFALSE;
    if (reject && x[0] > xmin && x[0] < xmax) 
      {
	TF1::RejectPoint();
	return 0;
      }
    Double_t fitval = TMath::Exp(par[0]+par[1]*x[0]);
    return fitval;
  }

  //Create an exponential fit function that skips a region.
  Double_t fit_exp_skip(Double_t *x,Double_t *par) 
  {
    Bool_t reject;
    reject = kTRUE;
    //reject = kFALSE;
    if (reject && x[0] > 3.3 && x[0] < 3.4) 
      {
	TF1::RejectPoint();
	return 0;
      }
    /*
      if (reject && x[0] > xmin && x[0] < xmax) 
      {
      TF1::RejectPoint();
      return 0;
      }
    */
    Double_t fitval = TMath::Exp(par[0]+par[1]*x[0]);
    return fitval;
  }

  //Create an exponential fit function to fit the background region, but exclude exponential peak and signifcant radiative tail of elastics. This fit will be binned to a hist and then added to the SIMC result to match experimental result.
  Double_t fit_exp_no_elastics(Double_t *x,Double_t *par) 
  {
    //Double_t xmin_no_elastics = 2.8, xmax_no_elastics = 3.15;
    Bool_t reject;
    reject = kTRUE;
    //reject = kFALSE;
    if (reject && x[0] > xmin_no_elastics && x[0] < xmax) 
      {
	TF1::RejectPoint();
	return 0;
      }
    Double_t fitval = TMath::Exp(par[0]+par[1]*x[0]);
    return fitval;
  }

  //Create Gaussian fit for the elastic peak.
  Double_t fit_gaus(Double_t *x,Double_t *par) 
  {
    Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((x[0]-par[1])/par[2]),2));
    return fitval;
  }
  //g1 = new TF1("g1","gaus",xmin,xmax);

  Double_t fit_lorentz(Double_t *x,Double_t *par)
  {
    Double_t fitval = (par[0]*par[2]) / ( pow(x[0]-par[1],2) + pow(0.5*par[2],2) );
    return fitval;
  }

  //Create a function that is the sum of the exponential background fit and the Gaussian elastic peak fit.
  Double_t fit_total(Double_t *x,Double_t *par) 
  {
    return fit_exp(x,par) + fit_gaus(x,&par[2]);
  }

  //Create a function that is the sum of the exponential background fit and the Gaussian elastic peak fit but skips over a region of the data. 
  Double_t fit_total_skip(Double_t *x,Double_t *par) 
  {
    Bool_t reject;
    reject = kTRUE;
    if (reject && x[0] > 2.8 && x[0] < 2.975) 
      {
	return fit_exp(x,par) + fit_gaus(x,&par[2]);
      }
    Double_t fitval = TMath::Exp(par[0]+par[1]*x[0]);
    return fitval;
  }

  //Create a function for fitting lines.
  Double_t fit_line(Double_t *x,Double_t *par) 
  {
    Bool_t reject;
    reject = kTRUE;
    //reject = kFALSE;
    if (reject && x[0] > 3.3 && x[0] < 3.4) 
      {
	TF1::RejectPoint();
	return 0;
      }
    return par[0]+par[1]*x[0];
  }

  //Create a function for fitting quadratic.
  Double_t fit_quad(Double_t *x,Double_t *par) 
  {
    Bool_t reject;
    reject = kTRUE;
    //reject = kFALSE;
    if (reject && x[0] > 3.3 && x[0] < 3.4) 
      {
	TF1::RejectPoint();
	return 0;
      }
    return par[0]+par[1]*x[0]+par[2]*pow(x[0],2);
  }

  //Create a function that is the sum of two exponentials (one before the elastic peak and one after it) and a Gaussian for the elastic peak.
  Double_t fit_total_new(Double_t *x,Double_t *par) 
  {
    return pow(TMath::Exp(par[0]+par[1]*x[0]),2) + fit_exp(x,&par[2]) + fit_gaus(x,&par[4]);
    //return TMath::Exp(par[0]+par[1]*x[0]) + TMath::Exp(par[2]+par[3]*x[0]) + fit_gaus(x,&par[4]);
    //return fit_exp(x,par) + fit_exp(x,&par[2]) + fit_gaus(x,&par[4]);
  }

 Double_t fit_total_lorentz(Double_t *x,Double_t *par) 
  {
    return fit_exp(x,par) + fit_lorentz(x,&par[2]);
  }

  //Create a function that draws the full fit including over the elastic peak region.
  Double_t fit_exp_full(Double_t *x,Double_t *par) 
  {
    Double_t fitval = TMath::Exp(par[0]+par[1]*x[0]);
    return fitval;
  }
  /*
  //Fit the histogram excluding the elastic peak.
  TF1 *func_exp = new TF1("func_exp",fit_exp,fitmin,fitmax,2);
  func_exp->SetLineColor(1);
  func_exp->SetParameter(0,12.);
  func_exp->SetParameter(1,-3.);
  h1->Fit("func_exp","R 0 M");
  cout<<"***** Exponential Fit: Chi^2 = "<<func_exp->GetChisquare()<<"   nDOF = "<<func_exp->GetNDF()<<"   Fit Probablility = "<<func_exp->GetProb()<<" *****"<<endl;
  
  //Fit the Gaussian to the elastic peak.
  TF1 *func_gaus = new TF1("func_gaus",fit_gaus,xmin,xmax,3);
  func_gaus->SetLineColor(3);
  func_gaus->SetParameter(0,1.);
  func_gaus->SetParameter(1,1.);
  func_gaus->SetParameter(2,1.);
  h1->Fit("func_gaus","R 0 M");
  func_gaus->Draw("same");
  //h1->Fit("g1","Rsame");
  //h1->GetFunction("g1")->SetLineColor(3);
  cout<<"***** Gaussian Fit: Chi^2 = "<<func_gaus->GetChisquare()<<"   nDOF = "<<func_gaus->GetNDF()<<"   Fit Probablility = "<<func_gaus->GetProb()<<" *****"<<endl;

  //Fit the combined background exponential and Gaussian elastic peak.
  TF1 *func_total = new TF1("func_total",fit_total,fitmin,fitmax,5);
  func_total->SetLineColor(2);
  func_total->SetNpx(1000);
  func_total->SetParameter(0,func_exp->GetParameter(0));
  func_total->SetParameter(1,func_exp->GetParameter(1));
  func_total->SetParameter(2,func_gaus->GetParameter(0));
  func_total->SetParameter(3,func_gaus->GetParameter(1));
  func_total->SetParameter(4,func_gaus->GetParameter(2));
  gStyle->SetOptFit(1111);
  h1->Fit("func_total","R M 0");
  func_total->SetLineWidth(6);
  func_total->Draw("same");
  cout<<"***** Combined Fit: Chi^2 = "<<func_total->GetChisquare()<<"   nDOF = "<<func_total->GetNDF()<<"   Fit Probablility = "<<func_total->GetProb()<<" *****"<<endl;
  */
  /*
  cout<<"*****Fit of exponential before elastic Peak*****"<<endl;
  //Fit the two exponentials and the Gaussian function to their respective regions.
  TF1 *func_exp_before = new TF1("func_exp_before",fit_exp,fitmin,xmin,2);
  func_exp_before->SetLineColor(3);
  func_exp_before->SetParameter(0,12.);
  func_exp_before->SetParameter(1,-3.);
  h1->Fit("func_exp_before","R M 0");
  func_exp_before->Draw("same");
  
  cout<<"*****Fit of exponential after elastic Peak*****"<<endl;
  TF1 *func_exp_after = new TF1("func_exp_after",fit_exp,xmax,fitmax,2);
  func_exp_after->SetLineColor(3);
  func_exp_after->SetParameter(0,12.);
  func_exp_after->SetParameter(1,-3.);
  h1->Fit("func_exp_after","R M 0");
  func_exp_after->Draw("same");

  
  //Fit the combined background exponential and Gaussian elastic peak.
  TF1 *func_total_new = new TF1("func_total_new",fit_total_new,fitmin,fitmax,7);
  func_total_new->SetLineColor(6);
  func_total_new->SetNpx(1000);
  func_total_new->SetParameter(0,func_exp_before->GetParameter(0));
  func_total_new->SetParameter(1,func_exp_before->GetParameter(1));
  func_total_new->SetParameter(2,func_exp_after->GetParameter(0));
  func_total_new->SetParameter(3,func_exp_after->GetParameter(1));
  func_total_new->SetParameter(4,func_gaus->GetParameter(0));
  func_total_new->SetParameter(5,func_gaus->GetParameter(1));
  func_total_new->SetParameter(6,func_gaus->GetParameter(2));
  //func_total_new->SetParLimits(4,func_gaus->GetParameter(0)-1,func_exp_after->GetParameter(0)+1);
  //func_total_new->SetParLimits(2,func_exp_after->GetParameter(0)-0.01,func_exp_after->GetParameter(0)+0.01);
  //func_total_new->SetParLimits(3,func_exp_after->GetParameter(1)-0.01,func_exp_after->GetParameter(1)+0.01);
  gStyle->SetOptFit(1111);
  h1->Fit("func_total_new","R M 0");
  func_total_new->Draw("same");
  */

  /*
    func_exp_before->SetParameter(0,func_total_new->GetParameter(0));
    func_exp_before->SetParameter(1,func_total_new->GetParameter(1));
    func_exp_before->SetLineColor(5);
    func_exp_before->Draw("same");
    func_exp_after->SetParameter(0,func_total_new->GetParameter(2));
    func_exp_after->SetParameter(1,func_total_new->GetParameter(3));
    func_exp_after->SetLineColor(5);
    func_exp_after->Draw("same");
  */
  /*
  //Plot the full exponential fit including the elstic peak if it was skipped over.
  TF1 *func_exp_full = new TF1("fit_exp_full",fit_exp_full,fitmin,fitmax,2);
  if(use_histo_method==1)
    {
      func_exp_full->SetParameter(0,func_exp->GetParameter(0));
      func_exp_full->SetParameter(1,func_exp->GetParameter(1));
      //func_exp_full->Draw("same");
      //cout<<"par[0] = "<<func_exp->GetParameter(0)<<"   par[1] = "<<func_exp->GetParameter(1)<<endl;
    }

  //Draw the exponential from the summed exponential background fit and Gaussian elastic peak fit.
  func_exp_full->SetParameter(0,func_total->GetParameter(0));
  func_exp_full->SetParameter(1,func_total->GetParameter(1));
  //func_exp_full->Draw("same");

  //Draw the Gaussian from the summed exponential background fit and Gaussian elastic peak fit.
  func_gaus->SetParameter(0,func_total->GetParameter(2));
  func_gaus->SetParameter(1,func_total->GetParameter(3));
  func_gaus->SetParameter(2,func_total->GetParameter(4));
  //func_gaus->Draw("same");
*/
  /*
  if(use_histo_method==1)
    {
      //Create another histogram with the same bin sizes in the region of the elastic peak with the same bin widths as the Xbj histo. This new histo can then be integrated and compare with the Xbj histogram.
      //Calculate bin width in the region of the elastic peak.
      TAxis *axis = h1->GetXaxis();
      Int_t bmin = axis->FindBin(xmin); 
      Int_t bmax = axis->FindBin(xmax);
      Double_t nbins = bmax-bmin;
      Double_t bin_width = (xmax-xmin)/nbins;
      cout<<"*********************************************"<<endl;
      cout<<"bmin = "<<bmin<<"   bmax = "<<bmax<<"   nbins = "<<nbins<<"   bin_width = "<<bin_width<<endl;
      TH1D *h2 = new TH1D("h2","Fit Xbj" , nbins, xmin, xmax);

      //Fill the fit histogram.
      for(Int_t i=0;i<nbins;i++)
	{
	  //h2->Fill(func_exp_full->Eval(h1->GetXaxis()->GetBinCenter(h1->GetXaxis()->FindBin(bmin+i))));
	  h2->SetBinContent(i, func_exp_full->Eval( h1->GetXaxis()->GetBinCenter(bmin+i) ) );
	  //cout<<"bin "<<bmin+i<<" has x value "<<h1->GetXaxis()->GetBinCenter(bmin+i)<<" and func_exp_full->Eval("<<h1->GetXaxis()->GetBinCenter(bmin+i)<<") = "<<func_exp_full->Eval( h1->GetXaxis()->GetBinCenter(bmin+i) )<<endl;
	}
      h2->Draw("same");
      h2->SetLineColor(kRed);
      
      //Find the integrals of the two histograms and subtract the fit histo integral from the Xbj histo to get the approximate number of events in the elastic peak.
      Double_t int_hist = h1->Integral(bmin,bmax);
      Double_t int_fit = h2->Integral(0,nbins-1);
      
      cout<<"The integral of the elastic peak region of the histogram is "<<int_hist<<". The integral of the fit function in the region of the elastic peak is "<<int_fit<<". Thus there are approximately "<<int_hist-int_fit<<" events in the elastic peak."<<endl;
      cout<<"*********************************************"<<endl;
    }
  
  //Print the integral of the Gaussian fitting the elastic peak to find the number of elastic electrons. Scale to number of bins.
  xb_binwidth = xb_nbins/(xbmax-xbmin);
  cout<<"xb_binwidth = "<<xb_binwidth<<endl;
  cout<<"There are "<<func_gaus->Integral(xmin,xmax)*xb_binwidth<<" elastic electrons in this fit."<<endl;
*/




  //Now subtract Al background from plot of Xbj.
  TChain *B = new TChain("T");
  B->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4076.root");
  B->SetBranchStatus("*",0);
  B->SetBranchStatus("EKL.x_bj",1);
  B->SetBranchStatus("ExtEKL.x_bj",1);
  B->SetBranchStatus("L.tr.tg_y",1);
  B->SetBranchStatus("L.tr.n",1);
  B->SetBranchStatus("L.cer.asum_c",1);
  B->SetBranchStatus("L.prl1.e",1);
  B->SetBranchStatus("L.prl2.e",1);
  B->SetBranchStatus("L.tr.tg_th",1);
  B->SetBranchStatus("L.tr.tg_ph",1);
  B->SetBranchStatus("L.tr.tg_dp",1);
  //B->SetBranchStatus("L.tr.p",1);
  B->SetBranchStatus("DBB.evtypebits",1);
  B->SetBranchStatus("right_bcm_u1c",1);
  B->SetBranchStatus("right_bcm_d1c",1);
  B->SetBranchStatus("right_clkcount",1);

  Double_t x_bj_Al=0.,Ext_x_bj_Al=0.,L_tr_tg_y_Al[21],L_tg_th_Al[21],L_tg_ph_Al[21],L_tg_dp_Al[21],L_cer_Al=0.,L_prl1_e_Al=0.,L_prl2_e_Al=0.,L_tr_p_Al=0.,L_tr_n_Al=0.,evtypebits_Al;
  Double_t u1c_Al=0.,d1c_Al=0.,rclk_Al=0.;
      
  B->SetBranchAddress("EKL.x_bj",&x_bj_Al);
  B->SetBranchAddress("ExtEKL.x_bj",&Ext_x_bj_Al);
  B->SetBranchAddress("L.tr.tg_y",&L_tr_tg_y_Al);
  B->SetBranchAddress("L.tr.n",&L_tr_n_Al);
  B->SetBranchAddress("L.cer.asum_c",&L_cer_Al);
  B->SetBranchAddress("L.prl1.e",&L_prl1_e_Al);
  B->SetBranchAddress("L.prl2.e",&L_prl2_e_Al);
  B->SetBranchAddress("L.tr.tg_th",L_tg_th_Al);
  B->SetBranchAddress("L.tr.tg_ph",L_tg_ph_Al);
  B->SetBranchAddress("L.tr.tg_dp",L_tg_dp_Al);
  //B->SetBranchAddress("L.tr.p",L_tr_p_Al);
  B->SetBranchAddress("DBB.evtypebits",&evtypebits_Al);
  B->SetBranchAddress("right_bcm_u1c",&u1c_Al);
  B->SetBranchAddress("right_bcm_d1c",&d1c_Al);
  B->SetBranchAddress("right_clkcount",&rclk_Al);
      
  Int_t nevts_Al = B->GetEntries();
  Int_t nevts_cuts_Al = 0;   //Count number of events surviving all of the cuts.
  cout<<"nevts_Al = "<<nevts_Al<<endl;

  //Create histogram of Al background. This histogram will be subtracted from the production Xbj histogram. Scale it to the charge of the combined production runs and also by the ratio of the radiative corrections for the Al dummy cell and the He3 cell Al walls.
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();
  c2->SetLogy();
  TH1D *hAl = new TH1D("hAl","Xbj Al" , xb_nbins, xbmin, xbmax);

  for(Int_t i=0;i<nevts_Al;i++) 
    {
      B->GetEntry(i);
      if(L_tr_tg_y_Al[0]>ymin && L_tr_tg_y_Al[0]<ymax && L_tr_n_Al==1 && (evtypebits_Al&1<<3)==1<<3 && L_prl1_e_Al>(-L_prl2_e_Al+pr_yint) && L_tg_ph_Al[0]>phmin && L_tg_ph_Al[0]<phmax && L_tg_th_Al[0]>thmin && L_tg_th_Al[0]<thmax && L_tg_dp_Al[0]>dpmin && L_tg_dp_Al[0]<dpmax)
	{
	  hAl->Fill(Ext_x_bj_Al);
	  nevts_cuts_Al++;
	}
    }
  cout<<"nevts_cuts_Al = "<<nevts_cuts_Al<<endl;
  hAl->Draw();

  //Now draw the production runs with the Al walls background subtracted out.
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();
  c3->SetLogy();

  TH1D *htot = new TH1D("htot","Xbj Minus Al Background" , xb_nbins, xbmin, xbmax);
  //h1->Add(hAl,-1.);
  htot->Add(h1,hAl,1.,-1.*charge*thickness*rc_dummy/rc_walls);
  htot->Draw();


  //Also read in SIMC elastic results to compare to the experimental elastic peak. Split into two ROOT files with different 3He densities to represent the boiling effects. Spliced together with a Y target cut.
  if(use_split_density==1)
    {
      //Add first SIMC ROOT file with initial density.
      TChain *SIMC1 = new TChain("h666");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_y50_dp6_x70_rho0.0345.root");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_final_cuts_rho0.0345_xs0.576893_5_29_18.root");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_final_cuts_rho0.0345_xs0.57299897225_5_29_18.root");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.0345_xs0.57299897225_5_30_18.root");
      SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.0345_xs0.59552_6_11_18.root");
      SIMC1->SetBranchStatus("*",0);
      SIMC1->SetBranchStatus("xbj",1);
      SIMC1->SetBranchStatus("ssytar",1);
      SIMC1->SetBranchStatus("ssytari",1);
      SIMC1->SetBranchStatus("ssyptar",1);
      SIMC1->SetBranchStatus("ssxptar",1);
      SIMC1->SetBranchStatus("ssdelta",1);
      SIMC1->SetBranchStatus("Weight",1);
      
      //Add second SIMC ROOT file with initial density.
      TChain *SIMC2 = new TChain("h666");
      //SIMC2->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_y50_dp6_x70_rho0.0233.root");
      //SIMC2->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_final_cuts_rho0.0233_xs0.576893_5_29_18.root");
      //SIMC2->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_final_cuts_rho0.0233_xs0.57299897225_5_29_18.root");
      //SIMC2->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.0233_xs0.57299897225_5_30_18.root");
      SIMC2->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.0233_xs0.59552_6_11_18.root");
      SIMC2->SetBranchStatus("*",0);
      SIMC2->SetBranchStatus("xbj",1);
      SIMC2->SetBranchStatus("ssytar",1);
      SIMC2->SetBranchStatus("ssytari",1);
      SIMC2->SetBranchStatus("ssyptar",1);
      SIMC2->SetBranchStatus("ssxptar",1);
      SIMC2->SetBranchStatus("ssdelta",1);
      SIMC2->SetBranchStatus("Weight",1);
      
      //Create histo for the first density ROOT file.
      TH1D *hSIMC1 = new TH1D("hSIMC1","SIMC First Density Xbj" , xb_nbins, xbmin, xbmax);
      hSIMC1->SetLineColor(3);
      if(show_histos==1)
	{
	  SIMC1->Draw("xbj>>hSIMC1",Form("Weight*%f/%d*(ssytar>%f&&ssytar<%f&&ssxptar>%f&&ssxptar<%f&&ssyptar>%f&&ssyptar<%f&&ssdelta>%f&&ssdelta<%f&&ssytari<%f)",Normfac1,nevts_SIMC,ymin_SIMC,ymax_SIMC,thmin_SIMC,thmax_SIMC,phmin_SIMC,phmax_SIMC,dpmin_SIMC,dpmax_SIMC,density_cut),"same");
	}
      else //Won't draw the histogram on the canvas but will still fill it.
	{
	  SIMC1->Draw("xbj>>hSIMC1",Form("Weight*%f/%d*(ssytar>%f&&ssytar<%f&&ssxptar>%f&&ssxptar<%f&&ssyptar>%f&&ssyptar<%f&&ssdelta>%f&&ssdelta<%f&&ssytari<%f)",Normfac1,nevts_SIMC,ymin_SIMC,ymax_SIMC,thmin_SIMC,thmax_SIMC,phmin_SIMC,phmax_SIMC,dpmin_SIMC,dpmax_SIMC,density_cut),"");
	  hSIMC1->Add(hSIMC1,1.);
	}
      
      //Create histo for the second density ROOT file.
      TH1D *hSIMC2 = new TH1D("hSIMC2","SIMC Second Density Xbj" , xb_nbins, xbmin, xbmax);
      hSIMC2->SetLineColor(8);
      if(show_histos==1)
	{
	  SIMC2->Draw("xbj>>hSIMC2",Form("Weight*%f/%d*(ssytar>%f&&ssytar<%f&&ssxptar>%f&&ssxptar<%f&&ssyptar>%f&&ssyptar<%f&&ssdelta>%f&&ssdelta<%f&&ssytari>%f)",Normfac2,nevts_SIMC,ymin_SIMC,ymax_SIMC,thmin_SIMC,thmax_SIMC,phmin_SIMC,phmax_SIMC,dpmin_SIMC,dpmax_SIMC,density_cut),"same");
	}
      else //Won't draw the histogram on the canvas but will still fill it.
    {
      SIMC2->Draw("xbj>>hSIMC2",Form("Weight*%f/%d*(ssytar>%f&&ssytar<%f&&ssxptar>%f&&ssxptar<%f&&ssyptar>%f&&ssyptar<%f&&ssdelta>%f&&ssdelta<%f&&ssytari>%f)",Normfac2,nevts_SIMC,ymin_SIMC,ymax_SIMC,thmin_SIMC,thmax_SIMC,phmin_SIMC,phmax_SIMC,dpmin_SIMC,dpmax_SIMC,density_cut),"");
      hSIMC2->Add(hSIMC2,1.);
    }
      
      cout<<"nevts hSIMC1 = "<<hSIMC1->GetEntries()<<"   nevts hSIMC2 = "<<hSIMC2->GetEntries()<<endl;
      
      //Now sum the two different density ROOT files together.
      TH1D *hSIMC = new TH1D("hSIMC","SIMC Combined Density Xbj" , xb_nbins, xbmin, xbmax);
      hSIMC->SetLineColor(kBlack);
      hSIMC->Add(hSIMC1,hSIMC2,1.,1.);
      hSIMC->Draw("same");   
    }
  else if       //If using only a single average density across the cell.
    {
      //Add SIMC ROOT file with average density.
      TChain *SIMC1 = new TChain("h666");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.01045_xs0.59552_6_19_18.root");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.01375_xs0.59552_6_19_18.root");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.01045_xs1.4133_6_20_18.root");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.01375_xs1.08312_6_20_18.root");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.013_xs1.08312_7_2_18.root");
      //SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.013_xs1.1881_7_2_18.root");
      SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_expanded_cuts_rho0.013_xs0.99206_7_5_18.root");
      SIMC1->SetBranchStatus("*",0);
      SIMC1->SetBranchStatus("xbj",1);
      SIMC1->SetBranchStatus("ssytar",1);
      SIMC1->SetBranchStatus("ssytari",1);
      SIMC1->SetBranchStatus("ssyptar",1);
      SIMC1->SetBranchStatus("ssxptar",1);
      SIMC1->SetBranchStatus("ssdelta",1);
      SIMC1->SetBranchStatus("Weight",1);

      //Create histo for the first density ROOT file.
      TH1D *hSIMC = new TH1D("hSIMC","SIMC Average Density Xbj" , xb_nbins, xbmin, xbmax);
      hSIMC->SetLineColor(3);
      if(show_histos==1)
	{
	  SIMC1->Draw("(xbj-0.02009)>>hSIMC",Form("Weight*%f/%d*(ssytar>%f&&ssytar<%f&&ssxptar>%f&&ssxptar<%f&&ssyptar>%f&&ssyptar<%f&&ssdelta>%f&&ssdelta<%f)",Normfac1,nevts_SIMC,ymin_SIMC,ymax_SIMC,thmin_SIMC,thmax_SIMC,phmin_SIMC,phmax_SIMC,dpmin_SIMC,dpmax_SIMC),"same");
	}
      else //Won't draw the histogram on the canvas but will still fill it.
	{
	  SIMC1->Draw("xbj>>hSIMC",Form("Weight*%f/%d*(ssytar>%f&&ssytar<%f&&ssxptar>%f&&ssxptar<%f&&ssyptar>%f&&ssyptar<%f&&ssdelta>%f&&ssdelta<%f)",Normfac1,nevts_SIMC,ymin_SIMC,ymax_SIMC,thmin_SIMC,thmax_SIMC,phmin_SIMC,phmax_SIMC,dpmin_SIMC,dpmax_SIMC),"");
	  hSIMC->Add(hSIMC,1.);
	}
    }
  
  //Now we want to fit the background subtracted histogram and find the number of elastic electrons.

  //Fit the histogram excluding the elastic peak.
  TF1 *func_exp_Al = new TF1("func_exp_Al",fit_exp,2.65,2.9,2);
  func_exp_Al->SetLineColor(1);
  func_exp_Al->SetParameter(0,12.);
  func_exp_Al->SetParameter(1,-3.);
  htot->Fit("func_exp_Al","R 0 M");
  cout<<"***** Exponential Fit: Chi^2 = "<<func_exp_Al->GetChisquare()<<"   nDOF = "<<func_exp_Al->GetNDF()<<"   Fit Probablility = "<<func_exp_Al->GetProb()<<" *****"<<endl;
  //func_exp_Al->Draw("same");
  //func_exp_Al->SetLineColor(5);

  //Fit the Gaussian to the elastic peak.
  TF1 *func_gaus_Al = new TF1("func_gaus_Al",fit_gaus,xmin,xmax,3);
  func_gaus_Al->SetLineColor(3);
  func_gaus_Al->SetParameter(0,1.);
  func_gaus_Al->SetParameter(1,1.);
  func_gaus_Al->SetParameter(2,1.);
  htot->Fit("func_gaus_Al","R 0 same M");
  //h1->Fit("g1","Rsame");
  //h1->GetFunction("g1")->SetLineColor(3);
  cout<<"***** Gaussian Fit: Chi^2 = "<<func_gaus_Al->GetChisquare()<<"   nDOF = "<<func_gaus_Al->GetNDF()<<"   Fit Probablility = "<<func_gaus_Al->GetProb()<<" *****"<<endl;

  //Fit the combined background exponential and Gaussian elastic peak.
  TF1 *func_total_Al = new TF1("func_total_Al",fit_total,2.5,fitmax,5);//2.5,3.25
  func_total_Al->SetLineColor(2);
  func_total_Al->SetNpx(1000);
  func_total_Al->SetParameter(0,func_exp_Al->GetParameter(0));
  //func_total_Al->SetParLimits(0,14.,18.);
  func_total_Al->SetParLimits(0,func_exp_Al->GetParameter(0),func_exp_Al->GetParameter(0));
  func_total_Al->SetParameter(1,func_exp_Al->GetParameter(1));
  //func_total_Al->SetParLimits(1,-4.4,-3.6);
  func_total_Al->SetParLimits(1,func_exp_Al->GetParameter(1),func_exp_Al->GetParameter(1));
  func_total_Al->SetParameter(2,func_gaus_Al->GetParameter(0));
  func_total_Al->SetParameter(3,func_gaus_Al->GetParameter(1));
  func_total_Al->SetParameter(4,func_gaus_Al->GetParameter(2));
  gStyle->SetOptFit(1111);
  htot->Fit("func_total_Al","R same M");
  cout<<"***** Combined Fit: Chi^2 = "<<func_total_Al->GetChisquare()<<"   nDOF = "<<func_total_Al->GetNDF()<<"   Fit Probablility = "<<func_total_Al->GetProb()<<" *****"<<endl;
  //Draw the total combined fit.
  func_total_Al->Draw("same");

  if(use_lorentz == 1)  
    {
      //Fit a Lorentzian to the elastic peak region.
      TF1 *func_lorentz = new TF1("func_lorentz",fit_lorentz,xmin,xmax,3);
      func_lorentz->SetParameter(0,2.37378);
      func_lorentz->SetParLimits(0,0.,100.);
      func_lorentz->SetParameter(1,3.02072);
      func_lorentz->SetParLimits(1,3.,3.04);
      func_lorentz->SetParameter(2,0.0217304);
      func_lorentz->SetParLimits(2,0.018,0.026);
      func_lorentz->SetLineColor(2);
      func_lorentz->SetNpx(1000);
      //func_lorentz->Draw("same");
      htot->Fit("func_lorentz","same");
      cout<<"***** Lorentz Fit: Chi^2 = "<<func_lorentz->GetChisquare()<<"   nDOF = "<<func_lorentz->GetNDF()<<"   Fit Probablility = "<<func_lorentz->GetProb()<<" *****"<<endl;
      
      //Fit a Lorentzian and an exponetial to the total fit.
      TF1 *func_total_lorentz = new TF1("func_total_lorentz",fit_total_lorentz,2.6,3.1,5);
      func_total_lorentz->SetParameter(0,func_exp_Al->GetParameter(0));
      func_total_lorentz->SetParLimits(0,18.,14.);
      func_total_lorentz->SetParameter(1,func_exp_Al->GetParameter(1));
      func_total_lorentz->SetParLimits(1,-3.6,-4.4);
      func_total_lorentz->SetParameter(2,1.5);//2.37378
      func_total_lorentz->SetParLimits(2,1.,2.);
      func_total_lorentz->SetParameter(3,3.02072);
      func_total_lorentz->SetParLimits(3,3.,3.04);
      func_total_lorentz->SetParameter(4,0.0217304);
      func_total_lorentz->SetParLimits(4,0.018,0.026);
      func_total_lorentz->SetLineColor(2);
      func_total_lorentz->SetNpx(1000);
      //func_total_lorentz->Draw("same");
      htot->Fit("func_total_lorentz","R same M");
      cout<<"***** Lorentz + exponential Fit: Chi^2 = "<<func_total_lorentz->GetChisquare()<<"   nDOF = "<<func_total_lorentz->GetNDF()<<"   Fit Probablility = "<<func_total_lorentz->GetProb()<<" *****"<<endl;
      func_total_lorentz->Draw("same");
    }
  
  /*
  cout<<"*****Fit of line after elastic peak with Al subtracted*****"<<endl;
  TF1 *func_line_Al_after = new TF1("func_line_Al_after",fit_line,xmax,fitmax,2);
  func_line_Al_after->SetLineColor(3);
  func_line_Al_after->SetParameter(0,40.);
  func_line_Al_after->SetParameter(1,-10.);
  htot->Fit("func_line_Al_after","R M 0");
  func_line_Al_after->Draw("same");
  */
  /*
    cout<<"*****Fit of quadratic after elastic peak with Al subtracted*****"<<endl;
    TF1 *func_quad_Al_after = new TF1("func_quad_Al_after",fit_quad,xmax,xbmax,3);
    func_quad_Al_after->SetLineColor(3);
    func_quad_Al_after->SetParameter(0,24.);
    func_quad_Al_after->SetParameter(1,-10.);
    func_quad_Al_after->SetParameter(2,1.);
    htot->Fit("func_quad_Al_after","R M 0");
    func_quad_Al_after->Draw("same");
  */
  
  /*
  //Fit with an exponential before and after the elastic peak plus a Gaussian for the peak.
  cout<<"*****Fit of exponential before elastic peak with Al subtracted*****"<<endl;
  //Fit the two exponentials and the Gaussian function to their respective regions.
  TF1 *func_exp_Al_before = new TF1("func_exp_Al_before",fit_exp,fitmin,xmin,2);
  func_exp_Al_before->SetLineColor(3);
  func_exp_Al_before->SetParameter(0,12.);
  func_exp_Al_before->SetParameter(1,-3.);
  htot->Fit("func_exp_Al_before","R M 0");
  func_exp_Al_before->Draw("same");
  
  cout<<"*****Fit of exponential after elastic peak with Al subtracted*****"<<endl;
  TF1 *func_exp_Al_after = new TF1("func_exp_Al_after",fit_exp_skip,xmax,xbmax,2);
  func_exp_Al_after->SetLineColor(3);
  func_exp_Al_after->SetParameter(0,12.);
  func_exp_Al_after->SetParameter(1,-3.);
  htot->Fit("func_exp_Al_after","R M 0");
  func_exp_Al_after->Draw("same");
  
  
  //Fit the combined background exponential and Gaussian elastic peak.
  TF1 *func_total_Al_new = new TF1("func_total_Al_new",fit_total_new,fitmin,fitmax,7);
  func_total_Al_new->SetLineColor(4);
  func_total_Al_new->SetNpx(1000);
  func_total_Al_new->SetParameter(0,func_exp_Al_before->GetParameter(0));
  func_total_Al_new->SetParameter(1,func_exp_Al_before->GetParameter(1));
  func_total_Al_new->SetParameter(2,func_exp_Al_after->GetParameter(0));
  func_total_Al_new->SetParameter(3,func_exp_Al_after->GetParameter(1));
  func_total_Al_new->SetParameter(4,func_gaus_Al->GetParameter(0));
  func_total_Al_new->SetParameter(5,func_gaus_Al->GetParameter(1));
  func_total_Al_new->SetParameter(6,func_gaus_Al->GetParameter(2));
  func_total_Al_new->SetParLimits(4,func_gaus_Al->GetParameter(0)-50.,func_gaus_Al->GetParameter(0)+50.);
  func_total_Al_new->SetParLimits(5,func_gaus_Al->GetParameter(1)-1.,func_gaus_Al->GetParameter(1)+1.);
  func_total_Al_new->SetParLimits(6,func_gaus_Al->GetParameter(2)-0.05,func_gaus_Al->GetParameter(2)+0.05);
  //func_total_Al_new->SetParLimits(2,func_line_Al_after->GetParameter(0)-0.0,func_line_Al_after->GetParameter(0)+0.0);
  //func_total_Al_new->SetParLimits(3,func_line_Al_after->GetParameter(1)-0.0,func_line_Al_after->GetParameter(1)+0.0);
  func_total_Al_new->SetParLimits(0,func_exp_Al_before->GetParameter(0)-10.0,func_exp_Al_before->GetParameter(0)+10.0);
  func_total_Al_new->SetParLimits(1,func_exp_Al_before->GetParameter(1)-3.0,func_exp_Al_before->GetParameter(1)+3.0);
  func_total_Al_new->SetParLimits(2,func_exp_Al_after->GetParameter(0)-0.0,func_exp_Al_after->GetParameter(0)+0.0);
  func_total_Al_new->SetParLimits(3,func_exp_Al_after->GetParameter(1)-0.0,func_exp_Al_after->GetParameter(1)+0.0);
  gStyle->SetOptFit(1111);
  htot->Fit("func_total_Al_new","R M 0");
  func_total_Al_new->Draw("same");
  */



  //Plot the full exponential fit including the elstic peak if it was skipped over.
  TF1 *func_exp_full_Al = new TF1("fit_exp_full_Al",fit_exp_full,fitmin,fitmax,2);
  if(use_histo_method==1)//This is immediately overwritten by the next line. Remove?
    {
      func_exp_full_Al->SetParameter(0,func_exp_Al->GetParameter(0));
      func_exp_full_Al->SetParameter(1,func_exp_Al->GetParameter(1));
      //func_exp_full->Draw("same");
      //cout<<"par[0] = "<<func_exp->GetParameter(0)<<"   par[1] = "<<func_exp->GetParameter(1)<<endl;
    }

  //Draw the exponential from the summed exponential background fit and Gaussian elastic peak fit.
  func_exp_full_Al->SetParameter(0,func_total_Al->GetParameter(0));
  func_exp_full_Al->SetParameter(1,func_total_Al->GetParameter(1));
  //func_exp_full->Draw("same");

  //Draw the Gaussian from the summed exponential background fit and Gaussian elastic peak fit.
  func_gaus_Al->SetParameter(0,func_total_Al->GetParameter(2));
  func_gaus_Al->SetParameter(1,func_total_Al->GetParameter(3));
  func_gaus_Al->SetParameter(2,func_total_Al->GetParameter(4));
  //func_gaus_Al->Draw("same");

  //Fit the exponential background without the elastic peak or main radiative tail to be summed with SIMC elastic results.
  TF1 *func_exp_Al_no_elastics = new TF1("func_exp_Al_no_elastics",fit_exp,2.4,2.7,2);
  func_exp_Al_no_elastics->SetLineColor(1);
  func_exp_Al_no_elastics->SetParameter(0,12.);
  func_exp_Al_no_elastics->SetParameter(1,-3.);
  htot->Fit("func_exp_Al_no_elastics","R M");
  cout<<"***** Exponential Fit Al Sub No Elastics: Chi^2 = "<<func_exp_Al_no_elastics->GetChisquare()<<"   nDOF = "<<func_exp_Al_no_elastics->GetNDF()<<"   Fit Probablility = "<<func_exp_Al_no_elastics->GetProb()<<" *****"<<endl;

  //Bin a histogram to the no elastics exponential fit that can be added to the SIMC elastic data to match the experimental result.
  TAxis *axis = htot->GetXaxis();
  Int_t bmin_no_elastics = axis->FindBin(fitmin); 
  Int_t bmax_no_elastics = axis->FindBin(fitmax);
  Double_t nbins_no_elastics = bmax_no_elastics-bmin_no_elastics;
  Double_t bin_width_no_elastics = (xmax_no_elastics-xmin_no_elastics)/nbins_no_elastics;
  cout<<"*********************************************"<<endl;
  cout<<"bmin_no_elastics = "<<bmin_no_elastics<<"   bmax_no_elastics = "<<bmax_no_elastics<<"   nbins_no_elastics = "<<nbins_no_elastics<<"   bin_width_no_elastics = "<<bin_width_no_elastics<<endl;
  TH1D *h_no_elastics = new TH1D("h_no_elastics","Fit Background No Elastics" , xb_nbins, xbmin, xbmax);

  //Remap the exponential fit skipping over a region to a full exponential for binning a histogram to the excluded region.
  TF1 *func_exp_full_Al_no_elastics = new TF1("fit_exp_full_Al_no_elastics",fit_exp_full,fitmin,fitmax,2);
  func_exp_full_Al_no_elastics->SetParameter(0,func_exp_Al_no_elastics->GetParameter(0));
  func_exp_full_Al_no_elastics->SetParameter(1,func_exp_Al_no_elastics->GetParameter(1));

  //Fill the fit histogram.
  for(Int_t i=0;i<xb_nbins;i++)
    {
      //cout<<"xbj = "<<htot->GetXaxis()->GetBinCenter(i)<<endl;
      if(subtract_SIMC == 0)
	{
	  if(i>bmin_no_elastics && i<bmax_no_elastics)
	    {
	      h_no_elastics->SetBinContent(i, func_exp_full_Al_no_elastics->Eval( htot->GetXaxis()->GetBinCenter(i) ) );
	    }
	  else
	    {
	      h_no_elastics->SetBinContent(i,0);
	    }
	}
      //If subtracting out the SIMC elastic results from the background first.
      if(subtract_SIMC == 1)
	{
	  if(i>bmin_no_elastics && i<axis->FindBin(xmin))
	    {
	      //h_no_elastics->SetBinContent(i, func_exp_full_Al_no_elastics->Eval( htot->GetXaxis()->GetBinCenter(i) ) );
	      //If subtracting out the SIMC elastic results from the background first.
	      h_no_elastics->SetBinContent(i, func_exp_full_Al_no_elastics->Eval( htot->GetXaxis()->GetBinCenter(i) ) - hSIMC->GetBinContent(i));
	      cout<<"Xbj = "<<htot->GetXaxis()->GetBinCenter(i)<<"   h_no_elastics = "<<func_exp_full_Al_no_elastics->Eval( htot->GetXaxis()->GetBinCenter(i) )<<"   hSIMC = "<<hSIMC->GetBinContent(i)<<"   h_no_elastics-SIMC = "<<func_exp_full_Al_no_elastics->Eval( htot->GetXaxis()->GetBinCenter(i) ) - hSIMC->GetBinContent(htot->GetXaxis()->GetBinCenter(i))<<endl;
	    }
	  else if(i>axis->FindBin(xmax) && i<bmax_no_elastics)
	    {
	      h_no_elastics->SetBinContent(i, func_exp_full_Al_no_elastics->Eval( htot->GetXaxis()->GetBinCenter(i) ) );
	    }
	  else
	    {
	      h_no_elastics->SetBinContent(i,0);
	    }
	}
    }
  if(show_histos==1)
    {
      h_no_elastics->Draw("same");
    }
  h_no_elastics->SetLineColor(1);

  //Now plot the sum of the SIMC results and the Al subtracted background so that it can be compared to the experimental data.
  TH1D *h_summed_SIMC = new TH1D("h_summed_SIMC","Fit Background No Elastics Summed with SIMC Elastic Results" , xb_nbins, xbmin, xbmax);
  h_summed_SIMC->Add(hSIMC,h_no_elastics,scale_SIMC,1.);
  h_summed_SIMC->Draw("same");
  h_summed_SIMC->SetLineColor(6);

  //Fit h_summed_SIMC histo with the same fit as the real data. Comparison of the magnitude of the Gaussians will allow use to calibrate the XS in the Monte Carlo.
  cout<<"*****Fit of SIMC plus one exponential fitting background in region skipping most elastics.*****"<<endl;
  TF1 *func_total_SIMC = new TF1("func_total_SIMC",fit_total,fitmin,fitmax,5);
  func_total_SIMC->SetLineColor(6);
  func_total_SIMC->SetNpx(1000);
  func_total_SIMC->SetParameter(0,func_exp_Al->GetParameter(0));
  func_total_SIMC->SetParameter(1,func_exp_Al->GetParameter(1));
  //func_total_SIMC->SetParLimits(0,10.,20.);
  //func_total_SIMC->SetParLimits(1,-2.,-7.);
  func_total_SIMC->SetParameter(2,0.5*func_gaus_Al->GetParameter(0));
  //func_total_SIMC->SetParLimits(2,200.,600.);
  func_total_SIMC->SetParameter(3,func_gaus_Al->GetParameter(1));
  //func_total_SIMC->SetParLimits(3,2.9,3.1);
  func_total_SIMC->SetParameter(4,func_gaus_Al->GetParameter(2));
  h_summed_SIMC->Fit("func_total_SIMC","R same M");
  func_total_SIMC->Draw("same");
  cout<<"***** Total Fit Al Sub Summed SIMC: Chi^2 = "<<func_total_SIMC->GetChisquare()<<"   nDOF = "<<func_total_SIMC->GetNDF()<<"   Fit Probablility = "<<func_total_SIMC->GetProb()<<" *****"<<endl;

  //Scale the SIMC results until the Gaussian SIMC peak nearly matches the data elastic peak height.
  TF1 *clear = new TF1("clear","0.",xbmin,xbmax);  //Function defined to clear histos. (multiply them by zero).
  if(match_data==1)
    {
      if(gaus_or_total==0)
	{
	  while(func_total_SIMC->GetParameter(2)>(func_total_Al->GetParameter(2)+0.5) || func_total_SIMC->GetParameter(2)<(func_total_Al->GetParameter(2)-0.5))
	    {
	      if(func_total_SIMC->GetParameter(2)>(func_total_Al->GetParameter(2)+0.5))
		{
		  cout<<"Scale factor for SIMC data+ = "<<scale_SIMC<<",   Gaussian data par[2] = "<<func_total_Al->GetParameter(2)<<",   Gaussian SIMC par[2] = "<<func_total_SIMC->GetParameter(2)<<endl;
		  if(fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2))>100)
		    {
		      scale_SIMC = scale_SIMC - 0.1;
		    }
		  if(fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2))<100 && fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2))>10)
		    {
		      //cout<<"fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2)) = "<<fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2))<<endl;
		      scale_SIMC = scale_SIMC - 0.01;
		    }
		  if(fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2))<10)
		    {
		      scale_SIMC = scale_SIMC - 0.001;
		    }
		  cout<<"Updated scale_SIMC = "<<scale_SIMC<<endl;
		  h_summed_SIMC->Multiply(clear,1.);    //Clear h_summed_SIMC.
		  h_summed_SIMC->Add(hSIMC,h_no_elastics,scale_SIMC,1.);    //Refill h_summed_SIMC using a new scale factor for the SIMC data.
		  h_summed_SIMC->Fit("func_total_SIMC","R same M q");
		  cout<<"Scale factor for SIMC data+ = "<<scale_SIMC<<",   Gaussian data par[2] = "<<func_total_Al->GetParameter(2)<<",   Gaussian SIMC par[2] = "<<func_total_SIMC->GetParameter(2)<<endl;
		  cout<<"*************************************************************************"<<endl;
		}
	      else if(func_total_SIMC->GetParameter(2)<(func_total_Al->GetParameter(2)-0.5))
		{
		  //cout<<"Scale factor for SIMC data- = "<<scale_SIMC<<",   Gaussian data par[2] = "<<func_total_Al->GetParameter(2)<<",   Gaussian SIMC par[2] = "<<func_total_SIMC->GetParameter(2)<<endl;
		  if(fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2))>100)
		    {
		      scale_SIMC = scale_SIMC + 0.1;
		    }
		  if(fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2))<100 && fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2))>10)
		    {
		      scale_SIMC = scale_SIMC + 0.01;
		    }
		  if(fabs(func_total_Al->GetParameter(2)-func_total_SIMC->GetParameter(2))<10)
		    {
		      scale_SIMC = scale_SIMC + 0.001;
		    }
		  //cout<<"Updated scale_SIMC = "<<scale_SIMC<<endl;
		  h_summed_SIMC->Multiply(clear,1.);    //Clear h_summed_SIMC.
		  h_summed_SIMC->Add(hSIMC,h_no_elastics,scale_SIMC,1.);    //Refill h_summed_SIMC using a new scale factor for the SIMC data.
		  h_summed_SIMC->Fit("func_total_SIMC","R same M q");
		  //cout<<"*************************************************************************"<<endl;
		}
	    }
	}

      //Using total fit peak matching. (Not just Gaussian portion of the total fit.
      if(gaus_or_total==1)
	{
	  while(func_total_SIMC->GetMaximum(xmin,xmax)>(func_total_Al->GetMaximum(xmin,xmax)+0.5) || func_total_SIMC->GetMaximum(xmin,xmax)<(func_total_Al->GetMaximum(xmin,xmax)-0.5))
	    {
	      if(func_total_SIMC->GetMaximum(xmin,xmax)>(func_total_Al->GetMaximum(xmin,xmax)+0.5))
		{
		  cout<<"Scale factor for SIMC data+ = "<<scale_SIMC<<",   total fit height = "<<func_total_Al->GetMaximum(xmin,xmax)<<",   total fit height  = "<<func_total_SIMC->GetMaximum(xmin,xmax)<<endl;
		  if(fabs(func_total_Al->GetMaximum(xmin,xmax)-func_total_SIMC->GetMaximum(xmin,xmax))>100)
		    {
		      scale_SIMC = scale_SIMC - 0.1;
		    }
		  if(fabs(func_total_Al->GetMaximum(xmin,xmax)-func_total_SIMC->GetMaximum(xmin,xmax))<100 && fabs(func_total_Al->GetMaximum(xmin,xmax)-func_total_SIMC->GetMaximum(xmin,xmax))>10)
		    {
		      scale_SIMC = scale_SIMC - 0.01;
		    }
		  if(fabs(func_total_Al->GetMaximum(xmin,xmax)-func_total_SIMC->GetMaximum(xmin,xmax))<10)
		    {
		      scale_SIMC = scale_SIMC - 0.001;
		    }
		  cout<<"Updated scale_SIMC = "<<scale_SIMC<<endl;
		  h_summed_SIMC->Multiply(clear,1.);    //Clear h_summed_SIMC.
		  h_summed_SIMC->Add(hSIMC,h_no_elastics,scale_SIMC,1.);    //Refill h_summed_SIMC using a new scale factor for the SIMC data.
		  /*
  func_total_SIMC->SetParameter(0,func_total_SIMC->GetParameter(0));
  func_total_SIMC->SetParameter(1,func_total_SIMC->GetParameter(1));
  //func_total_SIMC->SetParLimits(0,10.,20.);
  //func_total_SIMC->SetParLimits(1,-2.,-7.);
  func_total_SIMC->SetParameter(2,func_total_SIMC->GetParameter(2));
  //func_total_SIMC->SetParLimits(2,200.,600.);
  func_total_SIMC->SetParameter(3,func_total_SIMC->GetParameter(3));
  //func_total_SIMC->SetParLimits(3,2.9,3.1);
  func_total_SIMC->SetParameter(4,func_total_SIMC->GetParameter(4));
		  */
		  h_summed_SIMC->Fit("func_total_SIMC","R same M q");
		  cout<<"*************************************************************************"<<endl;
		}	      else if(func_total_SIMC->GetMaximum(xmin,xmax)<(func_total_Al->GetMaximum(xmin,xmax)-0.5))
		{
		  //cout<<"Scale factor for SIMC data- = "<<scale_SIMC<<",   total fit height = "<<func_total_Al->GetMaximum(xmin,xmax)<<",   total fit height = "<<func_total_SIMC->GetMaximum(xmin,xmax)<<endl;
		  if(fabs(func_total_Al->GetMaximum(xmin,xmax)-func_total_SIMC->GetMaximum(xmin,xmax))>100)
		    {
		      scale_SIMC = scale_SIMC + 0.1;
		    }
		  if(fabs(func_total_Al->GetMaximum(xmin,xmax)-func_total_SIMC->GetMaximum(xmin,xmax))<100 && fabs(func_total_Al->GetMaximum(xmin,xmax)-func_total_SIMC->GetMaximum(xmin,xmax))>10)
		    {
		      scale_SIMC = scale_SIMC + 0.01;
		    }
		  if(fabs(func_total_Al->GetMaximum(xmin,xmax)-func_total_SIMC->GetMaximum(xmin,xmax))<10)
		    {
		      scale_SIMC = scale_SIMC + 0.001;
		    }
		  //cout<<"Updated scale_SIMC = "<<scale_SIMC<<endl;
		  h_summed_SIMC->Multiply(clear,1.);    //Clear h_summed_SIMC.
		  h_summed_SIMC->Add(hSIMC,h_no_elastics,scale_SIMC,1.);    //Refill h_summed_SIMC using a new scale factor for the SIMC data.
		  /*
  func_total_SIMC->SetParameter(0,func_total_SIMC->GetParameter(0));
  func_total_SIMC->SetParameter(1,func_total_SIMC->GetParameter(1));
  //func_total_SIMC->SetParLimits(0,10.,20.);
  //func_total_SIMC->SetParLimits(1,-2.,-7.);
  func_total_SIMC->SetParameter(2,func_total_SIMC->GetParameter(2));
  //func_total_SIMC->SetParLimits(2,200.,600.);
  func_total_SIMC->SetParameter(3,func_total_SIMC->GetParameter(3));
  //func_total_SIMC->SetParLimits(3,2.9,3.1);
  func_total_SIMC->SetParameter(4,func_total_SIMC->GetParameter(4));
		  */
		  h_summed_SIMC->Fit("func_total_SIMC","R same M q");
		  //cout<<"*************************************************************************"<<endl;
		}
	    }
	}


      TF1 *func_total_gaus_data = new TF1("func_total_gaus_data",fit_gaus,xmin,xmax,3);
      func_total_gaus_data->SetParameter(0,func_total_Al->GetParameter(2));
      func_total_gaus_data->SetParameter(1,func_total_Al->GetParameter(3));
      func_total_gaus_data->SetParameter(2,func_total_Al->GetParameter(4));
      TF1 *func_total_gaus_SIMC = new TF1("func_total_gaus_SIMC",fit_gaus,xmin,xmax,3);
      func_total_gaus_SIMC->SetParameter(0,func_total_SIMC->GetParameter(2));
      func_total_gaus_SIMC->SetParameter(1,func_total_SIMC->GetParameter(3));
      func_total_gaus_SIMC->SetParameter(2,func_total_SIMC->GetParameter(4));
      Double_t area_data = (func_total_gaus_data->Integral(xmin,xmax)/xb_binwidth) * DT_correction * GC_eff;
      Double_t area_SIMC = func_total_gaus_SIMC->Integral(xmin,xmax)/xb_binwidth;

      //Match the areas of the Gaussian portions of both the data and SIMC total fits.
      if(gaus_or_total==2)
	{
	  cout<<"Area of data Gaussian before scaling = "<<area_data<<"   Area of SIMC Gaussian before scaling = "<<area_SIMC<<endl;

	  while(area_SIMC>(area_data+0.5) || area_SIMC<(area_data-0.5))
	    {
	      if(area_SIMC>(area_data+0.5))
		{
		  cout<<"Scale factor for SIMC+ = "<<scale_SIMC<<",   total data fit Gaussian area = "<<area_data<<",   total SIMC fit Gaussian area  = "<<area_SIMC<<endl;
		  if(fabs(area_data-area_SIMC)>100)
		    {
		      scale_SIMC = scale_SIMC - 0.1;
		    }
		  if(fabs(area_data-area_SIMC)<=100 && fabs(area_data-area_SIMC)>=10)
		    {
		      scale_SIMC = scale_SIMC - 0.01;
		    }
		  if(fabs(area_data-area_SIMC)<10)
		    {
		      scale_SIMC = scale_SIMC - 0.001;
		    }
		  cout<<"Updated scale_SIMC = "<<scale_SIMC<<endl;
		  h_summed_SIMC->Multiply(clear,1.);    //Clear h_summed_SIMC.
		  h_summed_SIMC->Add(hSIMC,h_no_elastics,scale_SIMC,1.);    //Refill h_summed_SIMC using a new scale factor for the SIMC data.
		  h_summed_SIMC->Fit("func_total_SIMC","R same M q");
		  //Redefine the Gaussian part of the new SIMC fit.
		  func_total_gaus_SIMC->SetParameter(0,func_total_SIMC->GetParameter(2));
		  func_total_gaus_SIMC->SetParameter(1,func_total_SIMC->GetParameter(3));
		  func_total_gaus_SIMC->SetParameter(2,func_total_SIMC->GetParameter(4));
		  area_SIMC = func_total_gaus_SIMC->Integral(xmin,xmax)/xb_binwidth;

		  cout<<"*************************************************************************"<<endl;
		}
	      else if(area_SIMC<(area_data-0.5))
		{
		  cout<<"Scale factor for SIMC+ = "<<scale_SIMC<<",   total data fit Gaussian area = "<<area_data<<",   total SIMC fit Gaussian area  = "<<area_SIMC<<endl;
		  if(fabs(area_data-area_SIMC)>100)
		    {
		      scale_SIMC = scale_SIMC + 0.1;
		    }
		  if(fabs(area_data-area_SIMC)<=100 && fabs(area_data-area_SIMC)>=10)
		    {
		      scale_SIMC = scale_SIMC + 0.01;
		    }
		  if(fabs(area_data-area_SIMC)<10)
		    {
		      scale_SIMC = scale_SIMC + 0.001;
		    }
		  cout<<"Updated scale_SIMC = "<<scale_SIMC<<endl;
		  h_summed_SIMC->Multiply(clear,1.);    //Clear h_summed_SIMC.
		  h_summed_SIMC->Add(hSIMC,h_no_elastics,scale_SIMC,1.);    //Refill h_summed_SIMC using a new scale factor for the SIMC data.
		  h_summed_SIMC->Fit("func_total_SIMC","R same M q");
		  //Redefine the Gaussian part of the new SIMC fit.
		  func_total_gaus_SIMC->SetParameter(0,func_total_SIMC->GetParameter(2));
		  func_total_gaus_SIMC->SetParameter(1,func_total_SIMC->GetParameter(3));
		  func_total_gaus_SIMC->SetParameter(2,func_total_SIMC->GetParameter(4));
		  area_SIMC = func_total_gaus_SIMC->Integral(xmin,xmax)/xb_binwidth;
		  cout<<"*************************************************************************"<<endl;
		}
	    }//End while loop.
	}//End Gaussian area matching function.

      //cout<<"Max value of elastic peak for data total fit = "<<func_total_Al->GetMaximum(xmin,xmax)<<endl;
      //cout<<"Max value of elastic peak for SIMC total fit = "<<func_total_SIMC->GetMaximum(xmin,xmax)<<endl;

      //Within reasonable range of target height or area. Print the fit values using the final value of scale_SIMC.
      cout<<"*************************************************************************"<<endl;
      h_summed_SIMC->Multiply(clear,1.);    //Clear h_summed_SIMC.
      h_summed_SIMC->Add(hSIMC,h_no_elastics,scale_SIMC,1.);    //Refill h_summed_SIMC using a new scale factor for the SIMC data.
      h_summed_SIMC->Fit("func_total_SIMC","R same M");
      cout<<"***** Fit after SIMC data scaled to match experimental data.: Chi^2 = "<<func_total_SIMC->GetChisquare()<<"   nDOF = "<<func_total_SIMC->GetNDF()<<"   Fit Probablility = "<<func_total_SIMC->GetProb()<<" *****"<<endl;
    }//End if match data.
  
  if(gaus_or_total == 1)
    {
      cout<<"Scale factor for SIMC data = "<<scale_SIMC<<".   Multiplying the scale factor for the elastic peak by the GC efficiency correction of "<<GC_eff<<" yields a scale factor of "<<scale_SIMC*GC_eff<<". Multiplying this value by the DT correction of "<<DT_correction<<" yields a final scale factor of "<<scale_SIMC*GC_eff*DT_correction<<"."<<endl;
      cout<<"Gaussian part of total data fit par[2] = "<<func_total_Al->GetParameter(2)<<",   Gaussian part of total SIMC fit par[2] = "<<func_total_SIMC->GetParameter(2)<<endl;
      cout<<"Max height of total data fit in elastic peak region = "<<func_total_Al->GetMaximum(xmin,xmax)<<"   Max height of total SIMC fit in elastic peak region = "<<func_total_SIMC->GetMaximum(xmin,xmax)<<endl;
    }

  if(gaus_or_total == 2)
    {
      cout<<"Scale factor for SIMC data after factoring in DT and GC efficiency corrections = "<<scale_SIMC<<endl;
    }

  //Calculate number of elactrons in elastic peak region for various fits.
  //Set exponential function to the exponential from func_total_Al.
  func_exp_Al->SetParameter(0,func_total_Al->GetParameter(0));
  func_exp_Al->SetParameter(1,func_total_Al->GetParameter(1));
  TAxis *axis = htot->GetXaxis();
  Int_t bmin = axis->FindBin(xmin); 
  Int_t bmax = axis->FindBin(xmax);
  Double_t nbins = bmax-bmin;
  Double_t bin_width = (xmax-xmin)/nbins;

  cout<<"********************************"<<endl;
  cout<<"********************************"<<endl;
  cout<<"The beginning of the elastic peak is bin "<<bmin<<" and the end of the elastic peak is "<<bmax<<"."<<endl;
  //Number of electrons in peak region based on func_total_Al which uses fit_total (one exponetial and a Gaussian).
  //cout<<"The number of counts under the exponential background fit after the peak extrapolated back in to the peak region is "<<func_exp_Al_after->Integral(xmin,xmax)/bin_width<<"."<<endl;
  //cout<<"The number of counts under the linear background fit after the peak extrapolated back in to the peak region is "<<func_line_Al_after->Integral(xmin,xmax)/bin_width<<"."<<endl;
  cout<<"The number of counts under the exponential fit to the total background in the combined fit is "<<func_exp_Al->Integral(xmin,xmax)/xb_binwidth<<"."<<endl;
  cout<<"The total counts under the elastic peak using the total fit of data ignoring all backgrounds is "<<func_total_Al->Integral(xmin,xmax)/xb_binwidth<<"."<<endl;
  cout<<"Data Gaussian: height = "<<func_total_Al->GetParameter(2)<<"   center = "<<func_total_Al->GetParameter(3)<<"   standard deviation = "<<func_total_Al->GetParameter(4)<<endl;
  cout<<"The total counts under the elastic peak using the total fit of SIMC ignoring all backgrounds is "<<func_total_SIMC->Integral(xmin,xmax)/xb_binwidth<<"."<<endl;
  cout<<"SIMC Gaussian: height = "<<func_total_SIMC->GetParameter(2)<<"   center = "<<func_total_SIMC->GetParameter(3)<<"   standard deviation = "<<func_total_SIMC->GetParameter(4)<<endl;
  cout<<"********************************"<<endl;
  cout<<"********************************"<<endl;

  
  //Calculate number of elctrons with the histgram method.
  if(use_histo_method==1)
    {
      //Create another histogram with the same bin sizes in the region of the elastic peak with the same bin widths as the Xbj histo. This new histo can then be integrated and compare with the Xbj histogram.
      //Calculate bin width in the region of the elastic peak.
      TAxis *axis = htot->GetXaxis();
      Int_t bmin = axis->FindBin(xmin); 
      Int_t bmax = axis->FindBin(xmax);
      Double_t nbins = bmax-bmin;
      Double_t bin_width = (xmax-xmin)/nbins;
      cout<<"*********************************************"<<endl;
      cout<<"bmin = "<<bmin<<"   bmax = "<<bmax<<"   nbins = "<<nbins<<"   bin_width = "<<bin_width<<endl;
      TH1D *h4 = new TH1D("h4","Fit Xbj" , nbins, xmin, xmax);

      //Fill the fit histogram.
      for(Int_t i=0;i<nbins;i++)
	{
	  //h2->Fill(func_exp_full->Eval(h1->GetXaxis()->GetBinCenter(h1->GetXaxis()->FindBin(bmin+i))));
	  h4->SetBinContent(i, func_exp_full_Al->Eval( htot->GetXaxis()->GetBinCenter(bmin+i) ) );
	  //cout<<"bin "<<bmin+i<<" has x value "<<h1->GetXaxis()->GetBinCenter(bmin+i)<<" and func_exp_full->Eval("<<h1->GetXaxis()->GetBinCenter(bmin+i)<<") = "<<func_exp_full->Eval( h1->GetXaxis()->GetBinCenter(bmin+i) )<<endl;
	}
      h4->Draw("same");
      h4->SetLineColor(kRed);
      
      //Find the integrals of the two histograms and subtract the fit histo integral from the Xbj histo to get the approximate number of events in the elastic peak.
      Double_t int_hist_Al = htot->Integral(bmin,bmax);
      Double_t int_fit_Al = h4->Integral(0,nbins-1);
      
      cout<<"The integral of the elastic peak region of the histogram is "<<int_hist_Al<<". The integral of the fit function in the region of the elastic peak is "<<int_fit_Al<<". Thus there are approximately "<<int_hist_Al-int_fit_Al<<" events in the elastic peak."<<endl;
      cout<<"*********************************************"<<endl;
    }
  
  //Print the integral of the Gaussian fitting the elastic peak to find the number of elastic electrons. Scale to number of bins.
  //xb_binwidth = (xbmax-xbmin)/xb_nbins;
  cout<<"xb_binwidth = "<<xb_binwidth<<endl;
  //Find the number of electrons in the elastic peak from the total fit (func_total_Al).
  cout<<"There are "<<func_gaus_Al->Integral(xmin,xmax)/xb_binwidth<<" elastic electrons in just the Gaussian part of the total data fit."<<endl;
  //cout<<"There are "<<func_gaus_Al->Integral(xmin,xmax)/bin_width<<" elastic electrons in just the Gaussian part of the total data fit."<<endl;
  cout<<"Gaussian paramerters data: par[0] = "<<func_gaus_Al->GetParameter(0)<<"   par[1] = "<<func_gaus_Al->GetParameter(1)<<"   par[2] = "<<func_gaus_Al->GetParameter(2)<<endl;
  TF1 *func_gaus_SIMC = new TF1("func_gaus_SIMC",fit_gaus,xmin,xmax,3);
  func_gaus_SIMC->SetParameter(0,func_total_SIMC->GetParameter(2));
  func_gaus_SIMC->SetParameter(1,func_total_SIMC->GetParameter(3));
  func_gaus_SIMC->SetParameter(2,func_total_SIMC->GetParameter(4));
  cout<<"There are "<<func_gaus_SIMC->Integral(xmin,xmax)/xb_binwidth<<" elastic electrons in just the Gaussian part of the total SIMC fit."<<endl;
  cout<<"Gaussian paramerters SIMC: par[0] = "<<func_gaus_SIMC->GetParameter(0)<<"   par[1] = "<<func_gaus_SIMC->GetParameter(1)<<"   par[2] = "<<func_gaus_SIMC->GetParameter(2)<<endl;
  //cout<<"There are "<<func_total_Al->Integral(xmin,xmax)*xb_binwidth<<" elastic electrons in this total fit."<<endl;


  //Now draw the production runs with no Al subraction on top of the scaled Al background for comparison.
  TCanvas* c4=new TCanvas("c4");
  c4->SetGrid();
  c4->SetLogy();

  TH1D *hAlScaled = new TH1D("hAlScaled","Xbj and Scaled Al Background" , xb_nbins, xbmin, xbmax);
  //h1->Add(hAl,-1.);
  h1->Draw();
  hAlScaled->Add(hAl,1.*charge*thickness*rc_dummy/rc_walls);
  hAlScaled->Draw("same");
  hAlScaled->SetLineColor(kRed);
  cout<<"*********************************************"<<endl;
  cout<<"Total events in unsubtracted Xbj = "<<h1->GetEntries()<<".   Total events in scaled Al background = "<<hAlScaled->GetEntries()<<"."<<endl;
  cout<<"The number of events in the unsubtracted Xbj peak = "<<h1->Integral(h1->FindBin(xmin),h1->FindBin(xmax))<<". The number of events in the elastic peak region of the scaled Al background = "<<hAlScaled->Integral(hAlScaled->FindBin(xmin),hAlScaled->FindBin(xmax))<<"."<<endl;
  //htot->Add(h1,hAl,1.,-1.*charge*thickness*rc_dummy/rc_walls);
  //htot->Draw();


  st->Stop();
  cout<<"*********************************************"<<endl;
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
