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
    //Double_t d = 0.15E-3;
    Double_t sint = xx/sqrt(xx*xx+par[5]*par[5]);
    Double_t argc = k*par[2]*sint*0.5;
    Double_t args = k*par[3]*sint*0.5;
    Double_t sen = sin(args)/args;
    Double_t val = par[0]*pow(cos(argc),2)*pow(sen,2);
    return val;
}

void diffrazione(){

gStyle->SetOptStat(2210);
gStyle->SetOptFit(1111);
gStyle->SetFitFormat("5.4g");

TGraphErrors *graph = new TGraphErrors("b015errore.txt", "%lg %lg %*lg %lg"); 

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

TF1 *f = new TF1("f", funzione, 0.045, 0.07 ,6);
f->SetParameter(0,max_y);
f->SetParLimits(0,max_y-1,max_y+1);
f->SetParameter(1,max_x);
f->SetParLimits(1,max_x-0.001,max_x+0.001);
f->SetParameter(2,0.30E-3);
f->SetParLimits(2,0.29E-3,0.31E-3);
f->SetParameter(3,0.15E-3);
f->SetParLimits(3,0.14E-3,0.16E-3);
f->SetParameter(4,632.8E-9);
f->SetParLimits(4,622.8E-9,642.8E-9);
f->SetParameter(5,1.135);
f->SetParLimits(5,1.125,1.145);  
f->SetNpx(1000);

graph->SetMarkerStyle(7);
graph->SetMarkerColor(kBlue);
graph->GetXaxis()->SetTitle("Posizione (m)");
graph->GetXaxis()->SetTitleOffset(1.05);
graph->GetYaxis()->SetTitle("Intensita' luminosa (u arbitraria)");
graph->GetYaxis()->SetTitleOffset(1.05);
graph->SetTitle("Doppia fenditura 2");

graph->GetXaxis()->SetRangeUser(0.045, 0.07);
graph->Fit("f","RL");
graph->Draw("APE");

double sumy = 0;
   for (int j = 1; j < graph->GetN(); j++)
    {
        sumy += graph->GetPointY(j);
    }
std::cout << sumy << endl;

double media = sumy/graph->GetN();
std::cout << media << endl;

double den = 0;
for (int k = 1; k < graph->GetN(); k++) {
    den += pow(graph->GetPointY(k)-media, 2);
}
std::cout << den << endl;

double numi = 0;
for (int t = 1; t < 240; t++) 
{
    numi += pow(((graph->GetPointY(t)) - (f->Eval(graph->GetPointX(t)))), 2);
}

//Devo evitare la posizione centrale altrimenti il risultato sballa

double numf = 0;
for (int r = 241; r < graph->GetN(); r++) 
{
    numf += pow(((graph->GetPointY(r)) - (f->Eval(graph->GetPointX(r)))), 2);
}
double num = numi + numf;
double rapp = num/den;
double R = 1 - rapp;

std::cout << "Parametro R^2: " << R << endl;    //Si chiama residuo?

}