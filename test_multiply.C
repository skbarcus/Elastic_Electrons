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

Int_t use_histo_method = 1;
Int_t show_histos = 1;
Int_t runlist[6] = {3892, 3893, 3894, 4073, 4074, 4075};
//Int_t runlist[1] = {4075};
  
Double_t Normfac1 = 0.278306e10, Normfac2 = 0.188633e10;
Int_t nevts_SIMC = 100000;
  
Double_t charge = 21.2708;     //Scale the Al background to the charge of the production runs.
Double_t thickness = 0.1979;   //Scale the Al background down to the thickness of 3He cell.
Double_t rc_dummy = 1./0.548506, rc_walls = 1./0.681787;  //Scale Al background by the RC for dummy and 3He cells.
Double_t xb_nbins = 200.;
Double_t xbmin = 0., xbmax = 4.;
Double_t xb_binwidth = abs(xbmax-xbmin)/xb_nbins;;
//Double_t xmin = 2.92, xmax = 3.15;
Double_t xmin = 2.95, xmax = 3.10;
Double_t xmin_no_elastics = 2.8, xmax_no_elastics = 3.15;
//Double_t fitmin = 2.8, fitmax = 3.2;
Double_t fitmin = 2.3, fitmax = 3.15;
//Double_t ymin = -0.028, ymax = 0.028;
Double_t ymin = -0.03, ymax = 0.028;      //I think this should be changed to a L.tr.vz cut instead.
Double_t thmin = -0.042, thmax = 0.049;
Double_t phmin = -0.03, phmax = 0.03;
Double_t dpmin = -0.02, dpmax = 0.03;
Double_t ymin_SIMC = ymin*100., ymax_SIMC = ymax*100.;
Double_t thmin_SIMC = TMath::ATan(thmin), thmax_SIMC = TMath::ATan(thmax);
Double_t phmin_SIMC = TMath::ATan(phmin), phmax_SIMC = TMath::ATan(phmax);
Double_t dpmin_SIMC = dpmin*100., dpmax_SIMC = dpmax*100.;
Double_t density_cut = -0.01469; //Y target position where the density profile changes.
Double_t GC_eff = 203./196.;     //GC efficiency correction to the elastic electron peak height. 
Double_t scale_SIMC = 1.0;       //Scale factor for SIMC histogram.

void test_multiply() 
{
  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  //c1->SetLogy();

  //Also read in SIMC elastic results to compare to the experimental elastic peak. Split into two ROOT files with different 3He densities to represent the boiling effects. Spliced together with a Y target cut.
  //Add first SIMC ROOT file with initial density.
  TChain *SIMC1 = new TChain("h666");
  SIMC1->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_y50_dp6_x70_rho0.0345.root");
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
  SIMC2->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/3he_elastic_y50_dp6_x70_rho0.0233.root");
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
  //hSIMC1->GetYaxis()->SetRange(0.,300.);
  //hSIMC1->GetYaxis()->SetRangeUser(0.,300000.);
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
  cout<<"Max hSIMC1 = "<<hSIMC1->GetMaximum()<<endl;
  cout<<"Max hSIMC2 = "<<hSIMC2->GetMaximum()<<endl;

  //Now sum the two different density ROOT files together.
  TH1D *hSIMC = new TH1D("hSIMC","SIMC Combined Density Xbj" , xb_nbins, xbmin, xbmax);
  hSIMC->GetYaxis()->SetRange(0,300);
  hSIMC->SetLineColor(kBlack);
  hSIMC->Add(hSIMC1,hSIMC2,1.,1.);
  hSIMC->Draw("same");
  cout<<"Max before multiplication = "<<hSIMC->GetMaximum()<<endl;

  hSIMC->Multiply(hSIMC);    //Multiplies hSIMC by hSIMC.
  hSIMC->Draw("same");
  cout<<"Max after multiplication = "<<hSIMC->GetMaximum()<<endl;

  hSIMC->Multiply(hSIMC1,hSIMC2,3.,1.);    //Replaces hSIMC with hSIMC1*3 * hSIMC2*1.
  hSIMC->Draw("same");
  cout<<"Max after multiplication second method = "<<hSIMC->GetMaximum()<<endl;

  TF1 *f1 = new TF1("f1","0.",xbmin,xbmax);
  hSIMC->Multiply(f1,1.);    //Multiplies hSIMC by f1*1.
  hSIMC->Draw("same");
  cout<<"Max after multiplication third method = "<<hSIMC->GetMaximum()<<endl;
  cout<<"nevts = "<<hSIMC->GetEntries()<<endl;
}
