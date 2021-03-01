#include "TFile.h"
#include "TTree.h"
#include "TInterpreter.h"
#include <iostream>
#include <fstream>

{
////////////////////////////////////////////////////////////////////////////
//Ce code permet de convertir des spectres geant4 et MCNP en spectre avec
// Un FWHM. Cela permet de connaitre le comportement d'un détecteur. Pour 
// Rappel le FWHM peut se calculer comme : 
// FWHM = a + b (E + c E**2)**1/2
// Avec a, b et c qui dépendent du détecteur. 
// Le FWHM est lié au sigma par : FWHM ~ 2.35 sigma
///////////////////////////////////////////////////////////////////////////// 
TFile f("/path/file.root");
TCanvas* c1 = new TCanvas("c1", "My title",800,800);

TRandom *MT=  new TRandom3();
///////Facteur FWHM
double fa=0.0;
double fb=0.0;
double fc=0.0;
double FWHM=0.0;
double edepG4=0.0;
double edepMCNP=0.0;
/////////////////////////

//Read G4Spectrum
TDirectory *Dir_1=(TDirectory*)f.Get("histo");
TH1D * h1 = (TH1D*)Dir_1->Get("h1.5");

//Read MCNP Spectrum
std::ifstream myReadFile;
TNtuple* a= new TNtuple("a","a","ener:num:err");
a -> ReadFile("/path/file.dat");

// Definition des histogrames
TH1D *h2 = new TH1D("Geant4","Geant4",2000,0,1.6);
TH1D *h3 = new TH1D("MCNP","MCNP",2000,0,1.6);
//Spectre Geant4 
for (Int_t i=0;i<h1->GetSize();i++) {
    for(Int_t j=0;j<h1->GetBinContent(i);j++){
    edepG4+=h1->GetBinCenter(i);
    FWHM=fa+fb*sqrt(h1->GetBinCenter(i)+fc*h1->GetBinCenter(i)*h1->GetBinCenter(i));
    double sigma=FWHM/2.35;
    double r= MT ->Gaus();
    double ener= h1->GetBinCenter(i) + sigma* r;
    h2->Fill(ener);
    }
}

//Spectre MCNP
Int_t n=a->GetEntries();
Float_t num,cps,ener;
Int_t numpart=5000000;
a->SetBranchStatus("err*", 0);
//a->SetBranchStatus("ener*", 0);

for (Int_t k=0; k<n ;k++){
a->SetBranchAddress("num", &num);
a->GetEntry(k);
cps=num*numpart;
for(Int_t l=0; l<(int)cps ;l++){
    a->SetBranchAddress("ener", &ener);
    a->GetEntry(k);
    edepMCNP+=ener;
    FWHM=fa+fb*sqrt(ener+fc*ener*ener);
    double sigma=FWHM/2.35;
    double r= MT ->Gaus();
    double ene= ener + sigma* r;
    h3->Fill(ene);
    }
/////////////////////////////////////////////
}
h2->SetStats(0); //ne display pas la stat box
h2->SetMaximum(200);// Max y
h3->SetLineColor(kBlue);
h2->SetTitle("Spectre du Xenon 135:");
h2->GetXaxis()->SetTitle("Energie [MeV]");
h2->GetYaxis()->SetTitle("Coups");
h2->Draw();
gPad->Update();
h3->SetLineColor(kRed);
h3->Draw("same");
gPad->Update();

auto legend = new TLegend(0.9,0.9,0.68,0.8);
//legend->SetHeader("","C"); // option "C" allows to center the header
legend->AddEntry(h2,"Geant4","l");
legend->AddEntry(h3,"MCNP","l");
legend->Draw();
gPad->Update();
cout<<"Dose Geant4= "<<edepG4<<endl;
cout<<"Dose MCNP= "<<edepMCNP<<endl;
c1->SaveAs("noFWH.root");
}


