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

Double_t PI = 3.14159265359;
Double_t deg2rad = PI/180.;

Double_t M = 3.0160293;         //Mass of Helium 3 amu.
Double_t theta = 21.;           //Scattering angle of incident electron degrees.
//Double_t E = 3.104;             //Incident beam energy GeV.
Double_t E = 3.356;
Double_t Ep = 0.;               //Energy of deflected electron GeV.
Double_t dE = 0.155;              //Momentum acceptance of HRS.
Double_t dtheta = 0.03*1/deg2rad;    //Horizontal angular acceptance of HRS.

Double_t xmin = 0.;
Double_t xmax = E+0.2;
Double_t npdraw = 1000.;        //Number of points used to draw function.

void E_Prime_vs_Theta_Elastic() 
{

  //Double_t Eprime(Double_t *angle, Double_t *par)
  Double_t Eprime(Double_t *Ep, Double_t *par)
  {
    //Ep = E / ( 1 + (2*E/M)*pow(sin(angle[0]/2.),2.) );
    //return Ep;

    //theta = TMath::ASin(pow( (M/(2*E))*(E/Ep[0] - 1) ,0.5))*2 *1/deg2rad;
    theta = TMath::ACos(M/E+1-M/Ep[0]) * 1/deg2rad;
    return theta;
  }

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
  Double_t ypadmax = gPad->GetUymax();
  Double_t xpadmax = gPad->GetUxmax();
 
  //Line representing the lower momentum acceptance.
  TLine *line1 = new TLine(Ep-dE, 0, Ep-dE, ypadmax);
  line1->SetLineColor(kBlack);
  line1->SetLineWidth(2);
  line1->Draw("SAME");
  
  //Line representing the upper momentum acceptance.
  TLine *line2 = new TLine(Ep+dE, 0, Ep+dE, ypadmax);
  line2->SetLineColor(kBlack);
  line2->SetLineWidth(2);
  line2->Draw("SAME");

  //Line representing the lower angular acceptance.
  TLine *line3 = new TLine(0, theta-dtheta, xpadmax, theta-dtheta);
  line3->SetLineColor(kBlack);
  line3->SetLineWidth(2);
  line3->Draw("SAME");
  
  //Line representing the upper angular acceptance.
  TLine *line4 = new TLine(0, theta+dtheta, xpadmax, theta+dtheta);
  line4->SetLineColor(kBlack);
  line4->SetLineWidth(2);
  line4->Draw("SAME");
}
