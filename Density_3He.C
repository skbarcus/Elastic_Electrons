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

void Density_3He() 
{
  Double_t min1=-0.079, max1=-0.05;
  Double_t min2=-0.032, max2=0.071;

  FILE *fp1 = fopen("/home/skbarcus/Tritium/Analysis/He3/Elastic_Electrons/He3.rho","r");
  Int_t skip1 = 1;                          //Gives number of lines to skip at top of data file. 
  Int_t nlines1 = 0;                        //Counts number of lines in the data file. 
  Int_t ncols1;                             //Set how many columns of data we have in the data file.
  char* string1[1000];                          //Variable to read lines of the data file.
  Float_t zreact_temp=0.,rho0_temp=0.,rho20_temp=0.,rho60_temp=0.,rho120_temp=0.;
  Float_t zreact[60]={},rho[60]={};
  Int_t entries = 0;
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
	  ncols1 = fscanf(fp1,"%f %f %f %f %f",&zreact_temp,&rho0_temp,&rho20_temp,&rho60_temp,&rho120_temp);
	  if (ncols1 < 0) break;  
	  zreact[nlines1-skip1] = zreact_temp;
	  rho[nlines1-skip1] = rho120_temp; 
	  nlines1++;
	  entries++;
	}
    }
  cout<<"There are "<<entries<<" entries in the table."<<endl;
  for(Int_t i=0;i<60;i++)
    {
      cout<<"zreact["<<i<<"] = "<<zreact[i]<<"   rho["<<i<<"] = "<<rho[i]<<endl;
    }

  TH2D *h1 = new TH2D("h1","Density Profile" , 60, -0.122, 0.114,100,0.,0.04);

  for(Int_t i=0;i<60;i++)
    {
      h1->Fill(zreact[i],rho[i]);
    }

  TCanvas* c1=new TCanvas("c1");
  c1->SetGrid();

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.08);
  h1->GetXaxis()->SetTitle("Zreact (m)");
  h1->GetYaxis()->SetTitle("Rho (g/cm^{3})");
  h1->GetYaxis()->SetTitleOffset(1.3);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);
  h1->SetMarkerColor(2);
  h1->Draw("");

  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetLabelSize(0.035);
  h1->GetYaxis()->SetTitleSize(0.06);
  h1->GetYaxis()->SetTitleOffset(0.75);
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetLabelSize(0.05);
  h1->GetXaxis()->SetTitleSize(0.06);
  h1->GetXaxis()->SetTitleOffset(0.75);

  Double_t disc=0.0;
  Double_t avg1=0.0, avg2=0.0, avg3=0.0;
  Double_t sum1=0.0, sum2=0.0, sum3-0.0;
  Int_t n1=0,n2=0,n3=0;

  for(Int_t i=0;i<60;i++)
    {
      //Plateau 1.
      if(zreact[i]>min1 && zreact[i]<max1)
	{
	  sum1=sum1+rho[i];
	  n1++;
	}
      //Plateau 2.
      if(zreact[i]>min2 && zreact[i]<max2)
	{
	  sum2=sum2+rho[i];
	  n2++;
	}
    }
  avg1=sum1/n1;
  avg2=sum2/n2;
  disc=(max1+min2)/2.;
  cout<<"The average heights of the two plateau regions are: "<<avg1<<" for "<<min1<<"<zreact<"<<max1<<"  and "<<avg2<<" for "<<min2<<"<zreact<"<<max2<<"."<<endl;
  cout<<"The discontinuity for the step funtion is at zreact = "<<disc<<"."<<endl;

  TLine *line1 = new TLine(min1,avg1-0.0025,min1,avg1+0.0025);
  TLine *line2 = new TLine(max1,avg1-0.0025,max1,avg1+0.0025);
  TLine *line3 = new TLine(min2,avg2-0.0025,min2,avg2+0.0025);
  TLine *line4 = new TLine(max2,avg2-0.0025,max2,avg2+0.0025);
  TLine *line5 = new TLine(-0.1,avg1,disc,avg1);
  TLine *line6 = new TLine(disc,avg2,0.1,avg2);
  TLine *line7 = new TLine(disc,avg2,disc,avg1);
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");
  line5->SetLineColor(4);
  line5->SetLineWidth(2);
  line5->Draw("same");
  line6->SetLineColor(4);
  line6->SetLineWidth(2);
  line6->Draw("same");
  line7->SetLineColor(4);
  line7->SetLineWidth(2);
  line7->Draw("same");
}
