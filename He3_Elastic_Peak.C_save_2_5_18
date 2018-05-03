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

Double_t PI = 3.14159265359;
Double_t deg2rad = PI/180.;

Double_t M = 3.01603*0.931;//3.0160293;         //Mass of Helium 3 GeV.

Double_t E = 3.356;             //Incident beam energy GeV.
Double_t Ep = 0.;               //Energy of deflected electron GeV.
Double_t dE = 0.155;              //Momentum acceptance of HRS.
Double_t dtheta = 0.03*1/deg2rad;    //Horizontal angular acceptance of HRS.

Double_t xmin = 0.;
Double_t xmax = E+0.2;
Double_t npdraw = 1000.;        //Number of points used to draw function.


//TChain chain("T");
//chain.Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4074.root");
//chain.Process("He3_Elastic_Peak.C","fillList");
//chain.SetBranchStatus("*",1);


//TFile *_file0 = TFile::Open("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4074.root");
//TFile *_file1 = TFile::Open("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4075.root");


Double_t Eprime(Double_t *Ep, Double_t *par)
{
  //Ep = E / ( 1 + (2*E/M)*pow(sin(angle[0]/2.),2.) );
  //return Ep;
  
  //theta = TMath::ASin(pow( (M/(2*E))*(E/Ep[0] - 1) ,0.5))*2 *1/deg2rad;
  theta = TMath::ACos(M/E+1-M/Ep[0]) * 1/deg2rad;
  return theta;
}



void He3_Elastic_Peak() 
{
  TChain *T = new TChain("T");
  //T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3892.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3893.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_3894.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4073.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4074*.root");
  T->Add("/home/skbarcus/Tritium/Analysis/He3/Rootfiles/e08014_4075*.root");

  //Double_t Eprime(Double_t *angle, Double_t *par)

  TF1 *fEprime = new TF1("fEprime", Eprime, xmin, xmax+.0, 1);

  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();
  //c1->SetLogy();

  fEprime->SetNpx(npdraw);
  fEprime->Draw();
  c1->SetTitle("Theta vs. E'");
  fEprime->GetHistogram()->SetTitle("Theta vs. E'");
  fEprime->GetHistogram()->GetYaxis()->SetTitle("Theta (deg.)");
  fEprime->GetHistogram()->GetXaxis()->SetTitle("E' (GeV)");
  //fEprime->GetXaxis()->SetRange(0.,3000000.5);
  //fEprime->GetXaxis()->SetRangeUser(0.,3.5);
  //fEprime->SetAxisRange(0.,3.5,"X");

  //Set Ep and theta again to draw the accpetance rectangle.
  Ep = 3.055;
  theta = 21;
  //ymax = fEprime->GetHistogram()->GetYmax();
  //c1->Update();
  Double_t ypadmax1 = gPad->GetUymax();
  Double_t xpadmax1 = gPad->GetUxmax();
 
  //Line representing the lower momentum acceptance.
  TLine *line1 = new TLine(Ep-dE, 0, Ep-dE, ypadmax1);
  line1->SetLineColor(kBlack);
  line1->SetLineWidth(2);
  line1->Draw("SAME");
  
  //Line representing the upper momentum acceptance.
  TLine *line2 = new TLine(Ep+dE, 0, Ep+dE, ypadmax1);
  line2->SetLineColor(kBlack);
  line2->SetLineWidth(2);
  line2->Draw("SAME");

  //Line representing the lower angular acceptance.
  TLine *line3 = new TLine(0, theta-dtheta, xpadmax1, theta-dtheta);
  line3->SetLineColor(kBlack);
  line3->SetLineWidth(2);
  line3->Draw("SAME");
  
  //Line representing the upper angular acceptance.
  TLine *line4 = new TLine(0, theta+dtheta, xpadmax1, theta+dtheta);
  line4->SetLineColor(kBlack);
  line4->SetLineWidth(2);
  line4->Draw("SAME");

Double_t theta1 = 21.;           //Scattering angle of incident electron degrees.

  //Make a second canvas for plotting theta (phi) as a function of E' (momentum).
  TCanvas* c2=new TCanvas("c2");
  c2->SetGrid();
  //L.tr.ph vs. L.tr.p
  //T->Draw("(L.tr.ph*180/3.14159265359+21):L.tr.p>>h2(1000,2.,3.5,1000,17.,24.)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1 && EKL.x_bj>2.98 && EKL.x_bj<3.11");

  //L.gold.ph vs. L.tr.p
  //T->Draw("(L.gold.ph*180/3.14159265359+21):(L.tr.p)>>h2(1000,2.,3.5,1000,17.,24.)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1 && EKL.x_bj>2.98 && EKL.x_bj<3.11");

  //L.gold.ph vs. 3.356 gold E'
  //T->Draw("(L.gold.ph*180./3.14159265359)+21:(  (3.356+L.gold.dp)/( 1.+(2.*(3.356+L.gold.dp)/3.0160293)*pow(sin((L.gold.ph+21*(3.14159265359/180.))/2.),2.) )  )>>h2(1000,2.,3.5,1000,17.,24.)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1 && EKL.x_bj>2.98 && EKL.x_bj<3.11");

  //L.gold.ph vs. 3.356 gold E' radians m=2.8084*0.931 GeV. (NIST = 3.0163)
  //T->Draw("TMath::ACos(  (cos(21.*3.14159265359/180.)-L.gold.ph*sin(21.*3.14159265359/180.))/pow(1+pow(L.gold.th,2.)+pow(L.gold.ph,2.),0.5)  )*180./3.14159265359:(  (3.356+L.gold.dp)/( 1.+(2.*(3.356+L.gold.dp)/(2.8084*0.931))*pow(sin((L.gold.ph+21*(3.14159265359/180.))/2.),2.) )  )>>h2(1000,2.,3.5,1000,17.,24.)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1 && EKL.x_bj>2.98 && EKL.x_bj<3.11");


  //L.gold.ph vs. 3.356 gold E' radians m=2.8084*0.931 GeV. (NIST = 3.0163)
  //T->Draw("TMath::ACos(  (cos(21.*3.14159265359/180.)-L.gold.ph*sin(21.*3.14159265359/180.))/pow(1+pow(L.gold.th,2.)+pow(L.gold.ph,2.),0.5)  )*180./3.14159265359:(  (3.356+L.gold.dp)/( 1.+(2.*(3.356+L.gold.dp)/(3.01603*0.931))*pow(sin(( TMath::ACos(  (cos(21.*3.14159265359/180.)-L.gold.ph*sin(21.*3.14159265359/180.))/pow(1+pow(L.gold.th,2.)+pow(L.gold.ph,2.),0.5)  )  )/2.),2.) )  )>>h2(1000,2.,3.5,1000,17.,24.)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1 && EKL.x_bj>2.98 && EKL.x_bj<3.11");

  //Using Nilanga
  T->Draw("TMath::ACos(  (cos(21.*3.14159265359/180.)-L.gold.ph*sin(21.*3.14159265359/180.))/pow(1+pow(L.gold.th,2.)+pow(L.gold.ph,2.),0.5)  )*180./3.14159265359:3.05666452*(1+L.gold.dp)>>h2(1000,2.,3.5,1000,17.,24.)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1 && EKL.x_bj>2.98 && EKL.x_bj<3.11");


  //L.tr.tg_ph vs. 3.356 gold E' radians m=2.8084*0.931 GeV. (NIST = 3.0163)
  //T->Draw("TMath::ACos(  (cos(21.*3.14159265359/180.)-L.tr.tg_ph*sin(21.*3.14159265359/180.))/pow(1+pow(L.tr.tg_th,2.)+pow(L.tr.tg_ph,2.),0.5)  )*180./3.14159265359:(  (3.356+L.tr.tg_dp)/( 1.+(2.*(3.356+L.tr.tg_dp)/(2.8084*0.931))*pow(sin(( TMath::ACos(  (cos(21.*3.14159265359/180.)-L.tr.tg_ph*sin(21.*3.14159265359/180.))/pow(1+pow(L.tr.tg_th,2.)+pow(L.tr.tg_ph,2.),0.5)  )  )/2.),2.) )  )>>h2(1000,2.,3.5,1000,17.,24.)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1 && EKL.x_bj>2.98 && EKL.x_bj<3.11");

  //L.gold.ph vs. 3.356 gold E' degrees.
  //T->Draw("0. + TMath::ACos(  ( (cos(21.*3.14159265359/180.)-L.gold.ph*(180./3.14159265359)*sin(21.*3.14159265359/180.))/pow(1+pow(L.gold.th*(180./3.14159265359),2.)+pow(L.gold.ph*(180./3.14159265359),2.),0.5) )*1  )*180./3.14159265359:(  (3.356+L.gold.dp)/( 1.+(2.*(3.356+L.gold.dp)/3.0160293)*pow(sin((L.gold.ph+21*(3.14159265359/180.))/2.),2.) )  )>>h2(1000,2.,3.5,1000,17.,24.)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1 && EKL.x_bj>2.98 && EKL.x_bj<3.11");

  //L.tr.tg_ph vs. 3.104 L.tr.tg_dp E'
  //T->Draw("(L.tr.tg_ph*180./3.14159265359+21):(  (3.104+L.tr.tg_dp)/( 1.+(2.*(3.104+L.tr.tg_dp)/3.0160293)*pow(sin((L.gold.ph+21*(3.14159265359/180.))/2.),2.) )  )>>h2(1000,2.,3.5,1000,17.,24.)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1 && EKL.x_bj>2.98 && EKL.x_bj<3.11");

  //h2->GetXaxis()->SetTitle("L.tr.p (momentum GeV)");
  //h2->GetYaxis()->SetTitle("L.tr.ph (angle radians)");

  Double_t ypadmax2 = gPad->GetUymax();
  Double_t xpadmax2 = gPad->GetUxmax();
  Double_t ypadmin2 = gPad->GetUymin();
  Double_t xpadmin2 = gPad->GetUxmin();

  //Line representing the lower momentum acceptance.
  TLine *line5 = new TLine(Ep-dE, ypadmin2, Ep-dE, ypadmax2);
  line5->SetLineColor(kBlue);
  line5->SetLineWidth(2);
  line5->Draw("SAME");
  
  //Line representing the upper momentum acceptance.
  TLine *line6 = new TLine(Ep+dE, ypadmin2, Ep+dE, ypadmax2);
  line6->SetLineColor(kBlue);
  line6->SetLineWidth(2);
  line6->Draw("SAME");

  //Line representing the lower angular acceptance.
  TLine *line7 = new TLine(xpadmin2,21.-dtheta,xpadmax2,21.-dtheta);
  line7->SetLineColor(kBlue);
  line7->SetLineWidth(2);
  line7->Draw("SAME");
  
  //Line representing the upper angular acceptance.
  TLine *line8 = new TLine(xpadmin2,21.+dtheta,xpadmax2,21.+dtheta);
  line8->SetLineColor(kBlue);
  line8->SetLineWidth(2);
  line8->Draw("SAME");

  fEprime->Draw("SAME");

  //Make a third canvas for plotting theta (phi) as a function of E' (momentum).
  TCanvas* c3=new TCanvas("c3");
  c3->SetGrid();
  c3->SetLogy();

  //T->Draw("EKL.x_bj>>h2(1000,0.5,3.5)","L.cer.asum_c>60 && (L.prl1.e+L.prl2.e)/L.tr.p/1000>0.7 && L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1");

  T->Draw("EKL.x_bj>>h2(1000,0.5,3.5)","L.tr.tg_y>-0.028 && L.tr.tg_y<0.028 && L.tr.n==1 && (DBB.evtypebits>>3)&1");

  //T->Draw("EKL.x_bj>>h2(1000,0.5,3.5)","L.tr.n==1 && (DBB.evtypebits>>3)&1");
}
