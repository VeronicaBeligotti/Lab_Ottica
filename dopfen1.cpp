#include <iostream>
#include "TMath.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TFitResult.h"
#include "TROOT.h"

Double_t funzione(Double_t *x, Double_t *par){

    Double_t xx= TMath::Abs(x[0]- par[1]);
    Double_t k = 2*TMath::Pi()/par[4];
    //Double_t D = 0.30E-3;
    //Double_t L = 1.135;
    //Double_t d = 0.10E-3;
    Double_t sint = xx/sqrt(xx*xx+par[5]*par[5]);
    Double_t argc = k*par[2]*sint*0.5;
    Double_t args = k*par[3]*sint*0.5;
    Double_t sen = sin(args)/args;
    Double_t val = par[0]*pow(cos(argc),2)*pow(sen,2);
    return val;
}

void diffrazione(){

gStyle->SetOptStat(2210);
gStyle->SetOptFit(111);
gStyle->SetFitFormat("4.3g");

TGraphErrors *graph = new TGraphErrors("b010errore.txt", "%lg %lg %*lg %lg"); 

TCanvas *c1 = new TCanvas("c1","Prima combinazione di dati");

Double_t x, y;
    graph->GetPoint(0, x, y);         //Prende la x e la y del punto i
    Double_t max_x = x, max_y = y;
    for(int i = 1; i < graph->GetN(); i++)    //Finchè i è minore del numero totale di punti il for continua a girare
    { graph->GetPoint(i, x, y);
        if(y > max_y) {
           max_x = x;
           max_y = y;
        }
    }
    std::cout<<max_x<<"     "<<max_y<<'\n';

TF1 *f = new TF1("f", funzione, 0.045, 0.07 ,6);
f->SetParameter(0,max_y);
f->SetParLimits(0,max_y-1,max_y+1);
f->SetParName(0,"Imax");
f->SetParameter(1,max_x);
f->SetParLimits(1,max_x-0.001,max_x+0.001);
f->SetParameter(2,0.30E-3);
f->SetParLimits(2,0.29E-3,0.31E-3);
f->SetParameter(3,0.10E-3);
f->SetParLimits(3,0.09E-3,0.11E-3);
f->SetParameter(4,632.8E-9);
f->SetParLimits(4,622.8E-9,642.8E-9);
f->SetParameter(5,1.135);
f->SetParLimits(5,1.130,1.140);      
f->SetLineColor(kRed);
f->SetNpx(1000);
f->SetParName(1,"x");
f->SetParName(2,"D");
f->SetParName(3,"d");
f->SetParName(4,"Lambda");
f->SetParName(5,"L");

graph->SetMarkerStyle(7);
graph->SetMarkerColor(kBlue);
graph->GetXaxis()->SetTitle("Posizione (m)");
graph->GetXaxis()->SetTitleOffset(1.05);
graph->GetYaxis()->SetTitle("Intensita' luminosa (u arbitraria)");
graph->GetYaxis()->SetTitleOffset(1.05);
graph->SetTitle("Prima combinazione di dati");


graph->GetXaxis()->SetRangeUser(0.045, 0.07);
graph->Fit("f","RLL");
graph->Draw("APE");

TLegend *leg = new TLegend (.1, .7, .3, .9, "Legenda");
leg->SetFillColor(0);
leg->AddEntry(graph,"Punti");
leg->AddEntry(f,"Fit Globale");
leg->Draw("SAME");

}