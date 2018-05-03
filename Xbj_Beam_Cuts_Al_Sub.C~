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
Int_t runlist[6] = {3892, 3893, 3894, 4073, 4074, 4075};
//Int_t runlist[1] = {4075};

Double_t xb_nbins = 1000.;
Double_t xbmin = 0., xbmax = 4.;
Double_t xb_binwidth = 0.;
//Double_t xmin = 2.92, xmax = 3.15;
Double_t xmin = 2.95, xmax = 3.10;
//Double_t fitmin = 2.8, fitmax = 3.2;
Double_t fitmin = 2.3, fitmax = 3.7;
//Double_t ymin = -0.028, ymax = 0.028;
Double_t ymin = -0.03, ymax = 0.03;      //I think this should be changed to a L.tr.vz cut instead.
Double_t thmin = -0.04, thmax = 0.055;
Double_t phmin = -0.03, phmax = 0.03;


void Xbj_Beam_Cuts() 
{
  //Define a new stopwatch.
  TStopwatch *st=new TStopwatch();
  st->Start(kTRUE);

  TH1D *h1 = new TH1D("h1","Xbj" , xb_nbins, xbmin, xbmax);

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
      //T->SetBranchStatus("L.tr.p",1);
      T->SetBranchStatus("DBB.evtypebits",1);
      T->SetBranchStatus("right_bcm_u1c",1);
      T->SetBranchStatus("right_bcm_d1c",1);
      T->SetBranchStatus("right_clkcount",1);
      
      Double_t x_bj=0.,Ext_x_bj=0.,L_tr_tg_y[21],L_tg_th[21],L_tg_ph[21],L_cer=0.,L_prl1_e=0.,L_prl2_e=0.,L_tr_p=0.,L_tr_n=0.,evtypebits;
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

      TCanvas* c1=new TCanvas("c1");
      c1->SetGrid();
      c1->SetLogy();
      
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
	  if(L_tr_tg_y[0]>ymin && L_tr_tg_y[0]<ymax && L_tr_n==1 && (evtypebits&1<<3)==1<<3 && L_prl1_e>(-L_prl2_e+2000) && L_tg_ph[0]>phmin && L_tg_ph[0]<phmax && L_tg_th[0]>thmin && L_tg_th[0]<thmax)
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

  //Create Gaussian fit for the elastic peak.
  Double_t fit_gaus(Double_t *x,Double_t *par) 
  {
    Double_t fitval = par[0]*TMath::Exp(-0.5*pow(((x[0]-par[1])/par[2]),2));
    return fitval;
  }
  //g1 = new TF1("g1","gaus",xmin,xmax);

  //Create a function that is the sum of the exponential background fit and the Gaussian elastic peak fit.
  Double_t fit_total(Double_t *x,Double_t *par) 
  {
    return fit_exp(x,par) + fit_gaus(x,&par[2]);
  }

  //Create a function that draws the full fit including over the elastic peak region.
  Double_t fit_exp_full(Double_t *x,Double_t *par) 
  {
    Double_t fitval = TMath::Exp(par[0]+par[1]*x[0]);
    return fitval;
  }

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
  h1->Fit("func_gaus","R 0 same M");
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
  h1->Fit("func_total","R SAME M");
  cout<<"***** Combined Fit: Chi^2 = "<<func_total->GetChisquare()<<"   nDOF = "<<func_total->GetNDF()<<"   Fit Probablility = "<<func_total->GetProb()<<" *****"<<endl;

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

  st->Stop();
  cout<<"CPU time = "<<st->CpuTime()<<" s = "<<st->CpuTime()/60.<<" min   Real time = "<<st->RealTime()<<" s = "<<st->RealTime()/60.<<" min"<<endl;
}
