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
    Double_t k = 2*TMath::Pi()/632.8E-9;
    Double_t D = 0.30E-3;
    Double_t L = 1.135;
    Double_t d = 0.20E-3;
    Double_t sint = xx/sqrt(xx*xx+L*L);
    Double_t argc = k*D*sint*0.5;
    Double_t args = k*d*sint*0.5;
    Double_t sen = sin(args)/args;
    Double_t val = par[0]*pow(cos(argc),2)*pow(sen,2);
    return val;
}

void diffrazione(){

gStyle->SetOptStat(2210);
gStyle->SetOptFit(1111);
gStyle->SetFitFormat("7.6g");

TGraphErrors *graph = new TGraphErrors("b020errore.txt", "%lg %lg %*lg %lg"); 

TCanvas *c1 = new TCanvas("c1","Doppia Fenditura");

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

TF1 *f = new TF1("f", funzione, 0.045, 0.075 ,2);
f->SetParameter(0,max_y);
f->SetParLimits(0,max_y-1,max_y+1);
f->SetParameter(1,max_x);
f->SetParLimits(1,max_x-0.001,max_x+0.001);  
f->SetLineColor(kRed);

graph->SetMarkerStyle(7);
graph->SetMarkerColor(kBlue);
graph->GetXaxis()->SetTitle("Posizione (m)");
graph->GetYaxis()->SetTitle("Intensita' luminosa (u arbitraria)");
graph->SetTitle("Doppia fenditura 3");


graph->GetXaxis()->SetRangeUser(0.045, 0.075);
graph->Fit("f","R");
graph->Draw("APE");

}