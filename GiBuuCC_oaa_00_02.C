#define GiBuuCC_oaa_00_02_cxx
#include "GiBuuCC_oaa_00_02.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TVector3.h"
#include <vector>
#include  "TMath.h"
#include <iostream>
#include <fstream> 
using namespace std;

/// define some common parameters that may be used later

#define PI 3.14159265
#define mumass 0.1056583755 // GeV/c2
#define protonmass 0.93827208816 //GeV/c2
#define pionmass 0.1349768  //GeV/c2

/////////////////////////////// defining output files. one ROOT file (TFile format) and one text file (txt file) ///////////////////////////

ofstream listOAASBND_00_02CC_GiBUU;
TFile *outOAASBND_00_02CC_GiBUU;

//////////////////////////////////////////////////////////////////// DECLARING HISTOGRAMS /////////////////////////////////////////////////////////////////

///// 1-D histograms, TH1D ///////////////////////////////

// true particle multiplicity multiplicity (only visibles)

TH1D *trueNproton_CC_00_02;
TH1D *trueNpion_CC_00_02;
TH1D *trueNpi0_CC_00_02;
// true neutrino energy
TH1D *trueEnu_CC_00_02;
// true energy transfer
TH1D *trueOmega_CC_00_02;

TH1D *trueEnu_CCQE_00_02;
TH1D *trueEnu_CCMEC_00_02;
TH1D *trueEnu_CCRES_00_02;
TH1D *trueEnu_CCDIS_00_02;
TH1D *trueEnu_other_00_02;

///////////////////////////////////////////////////////////////////// END of DECLARING HISTROGRAMS ///////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////  SOME COUNTERS //////////////////////////////////////////////////////////////////

int totentries =0;/// number of total numu CC events selected
int nCCQE = 0;/// number of total numu CCQE events selected
int nCCMEC =0;/// number of total numu CCMEC events selected
int nCCRes =0;/// number of total numu CCRes events selected
int nCCDIS =0;/// number of total numu CCDIS events selected
int nother =0;/// just for debbuging

///////////////////////////////////////////////////////////////////////// FUNCTIONS //////////////////////////////////////////////////////////

double cosine(double px, double py , double pz, double nupx, double nupy , double nupz){

  double norm_P = sqrt(pow(px,2) + pow(py,2)+ pow(pz,2));
  double norm_nuP = sqrt(pow(nupx,2) + pow(nupy,2)+ pow(nupz,2));
  
  TVector3 pdir;
  pdir[0] = px;
  pdir[1] = py;
  pdir[2] = pz;

  TVector3 beamdir;
  beamdir[0] = nupx;
  beamdir[1] = nupy;
  beamdir[2] = nupz;

  double cos; // scalar product
  cos = pdir.Dot(beamdir);
  cos = cos/(norm_P*norm_nuP);

  return cos;

}

//// calculating the energy transfer from the neutrino to the nucleus. We could use as input reconstructred neutrino energy also to test resolutions
double w(double emu, double enu){ ///// energy transfer

  return (enu - emu);
  
}

void GiBuuCC_oaa_00_02::Loop()
{
//   In a ROOT session, you can do:
//      root> .L GiBuuCC_oaa_00_02.C
//      root> GiBuuCC_oaa_00_02 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   bool numuCC = false;

   Long64_t nentries = fChain->GetEntriesFast();

   //////////////// output files location and histogram definitions

   listOAASBND_00_02CC_GiBUU.open("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/info_SBND_CC_00_02.txt");
   outOAASBND_00_02CC_GiBUU = new TFile("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/analysis_SBND_CC_00_02.root","RECREATE");

   trueNproton_CC_00_02 = new TH1D("trueNproton_CC_00_02", "trueNproton_CC_00_02", 15, 0, 15);
   trueNpion_CC_00_02 = new TH1D("trueNpion_CC_00_02", "trueNpion_CC_00_02", 6, 0, 6);
   trueNpi0_CC_00_02 = new TH1D("trueNpi0_CC_00_02", "trueNpi0_CC_00_02", 6, 0, 6);

   trueEnu_CC_00_02 = new TH1D("trueEnu_CC_00_02", "trueEnu_CC_00_02", 60, 0, 3);
   trueOmega_CC_00_02 = new TH1D("trueOmega_CC_00_02", "trueOmega_CC_00_02", 60, 0, 3);
   trueEnu_CCQE_00_02 = new TH1D("trueEnu_CCQE_00_02", "trueEnu_CCQE_00_02", 60, 0, 3);
   trueEnu_CCMEC_00_02 = new TH1D("trueEnu_CCMEC_00_02", "trueEnu_CCMEC_00_02", 60, 0, 3);
   trueEnu_CCRES_00_02 = new TH1D("trueEnu_CCRES_00_02", "trueEnu_CCRES_00_02", 60, 0, 3);
   trueEnu_CCDIS_00_02 = new TH1D("trueEnu_CCDIS_00_02", "trueEnu_CCDIS_00_02", 60, 0, 3);
   trueEnu_other_00_02 = new TH1D("trueEnu_other_00_02", "trueEnu_other_00_02", 60, 0, 3);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      double pmu = sqrt(lepOut_Px*lepOut_Px + lepOut_Py*lepOut_Py + lepOut_Pz*lepOut_Pz); /// true muon (absolute) momentum

      if(lepIn_E>0 && pmu >0) numuCC = true;;/// this is the condition to select a numu CC event in true variables
      if (!numuCC) continue;
      
      totentries++; /// counting number of events selected as numu CC

      //couting number of particles int eh final state by type.
      int nneutron = 0; /// number of neutrons
      int npion = 0; /// number of charged pions
      int nproton = 0; /// number of protons
      int npi0 = 0;/// number of neutral pions
      int nkaon = 0;/// number of kaons (we don't care if negative or positive)
      int neta = 0;/// number of etas (we don't care if negative or positive)
      int nrho = 0;/// number of rhos (we don't care if negative or positive)

      /////// some true kinematics for all events //////

      double emu = lepOut_E; /// true muon energy
      double omega = w(emu, lepIn_E); /// true energy transfer
      double protonmom =0; //// momentum of the proton
      double pionmom = 0.;
      double Epionplus = 0.;
      double Epionneg = 0.;
      double Epion = 0.;
      double pi0mom = 0.;

      double greatest_proton=-1;
      double greatest_pion=-1;
      double greatest_pi0=-1;

      /////////////////////////
      for (int k=0;k<barcode->size();k++){
	if ((*barcode)[k]==2212){
	  greatest_proton = sqrt((*Px)[k]*(*Px)[k] + (*Py)[k]*(*Py)[k] + (*Pz)[k]*(*Pz)[k]);
	  if(greatest_proton < 0.256) continue; //// ************************** proton detection threshold ********************
	  nproton++;
	  if(greatest_proton > protonmom) protonmom = greatest_proton;
         }
	if (abs((*barcode)[k])==211 ){
	  greatest_pion = sqrt((*Px)[k]*(*Px)[k] + (*Py)[k]*(*Py)[k] + (*Pz)[k]*(*Pz)[k]);
	  if(greatest_pion < 0.053) continue; //// ************************** proton detection threshold ********************
	  npion++;
	  if(greatest_pion > pionmom) pionmom = greatest_pion;
         }
	if (abs((*barcode)[k])==111 ){
	  greatest_pi0 = sqrt((*Px)[k]*(*Px)[k] + (*Py)[k]*(*Py)[k] + (*Pz)[k]*(*Pz)[k]);
	  if(greatest_pi0 < 0.053) continue; //// ************************** proton detection threshold ********************
	  npi0++;
	  if(greatest_pi0 > pi0mom) pi0mom = greatest_pi0;
         }

	trueNproton_CC_00_02->Fill(nproton, weight);
	trueNpion_CC_00_02->Fill(npion, weight);
	trueNpi0_CC_00_02->Fill(npi0, weight);
	trueEnu_CC_00_02->Fill(lepIn_E, weight);
	trueOmega_CC_00_02->Fill(omega, weight);

	if (evType==1){
	  trueEnu_CCQE_00_02->Fill(lepIn_E, weight);
	  nCCQE++;
	}
	else if ((evType>=2)&&(evType<=31)){ ///// non-strange baryon resonances
      	  trueEnu_CCRES_00_02->Fill(lepIn_E, weight);
	  nCCRes++;
	}
	else if (evType==32){ /// pi-neutron bkg (nu+ neutron -> muon + piplus + neutron)
	  trueEnu_CCDIS_00_02->Fill(lepIn_E, weight);
	  nCCDIS++;
	}
	else if (evType==33){ ///// pi-proton bkg (nu + neutron -> muon + pi0 + proton)
	  trueEnu_CCDIS_00_02->Fill(lepIn_E, weight);
	  nCCDIS++;
	}
	else if (evType==34){ ///// DIS (the most common definition accross experiments)
	   trueEnu_CCDIS_00_02->Fill(lepIn_E, weight);
	   nCCDIS++;
	}
	else if (evType==35){ //// 2p2h QE
	  trueEnu_CCMEC_00_02->Fill(lepIn_E, weight);
	  nCCMEC++;
	}
	else if (evType==36){ //// 2p2h Delta
	  trueEnu_CCMEC_00_02->Fill(lepIn_E, weight);
	  nCCMEC++;
	}
	else if (evType==37){ ///// 2pi bkg
	  trueEnu_CCDIS_00_02->Fill(lepIn_E, weight);
	  nCCDIS++;
	}
	else {
	  trueEnu_other_00_02->Fill(lepIn_E, weight);
	  nother++;
	}
      }
      
   }
   /// scaling down to the multiplicity added per run
   trueNproton_CC_00_02->Scale(100);
   trueNpion_CC_00_02->Scale(100);
   trueNpi0_CC_00_02->Scale(100);
   trueEnu_CC_00_02->Scale(100);
   trueOmega_CC_00_02->Scale(100);
   trueEnu_CCQE_00_02->Scale(100);
   trueEnu_CCMEC_00_02->Scale(100);
   trueEnu_CCRES_00_02->Scale(100);
   trueEnu_CCDIS_00_02->Scale(100);
   trueEnu_other_00_02->Scale(100);
   
   ////////////////////// PLOTTING //////////////////////////////////////
   ///// write few information in the txt file, and write the ROOT file before closing
   
   listOAASBND_00_02CC_GiBUU<<" number of CC :  "<<totentries<<endl;
   listOAASBND_00_02CC_GiBUU<<" number of CCQE events in CC :  "<<nCCQE<<endl;
   listOAASBND_00_02CC_GiBUU<<" number of CCMEC events in CC :  "<<nCCMEC<<endl;
   listOAASBND_00_02CC_GiBUU<<" number of CCRes events in CC :  "<<nCCRes<<endl;
   listOAASBND_00_02CC_GiBUU<<" number of CCDIS events in CC :  "<<nCCDIS<<endl;
   listOAASBND_00_02CC_GiBUU<<" number of other events in CC :  "<<nother<<endl;

   listOAASBND_00_02CC_GiBUU<< "Output file written" << std::endl;

   outOAASBND_00_02CC_GiBUU->cd();
   outOAASBND_00_02CC_GiBUU->Write();

   ////////////////////////////////////////////////////////// DRAWING //////////////////////////

   /////// preliminary drawing style ////

   gStyle->SetOptStat(0000);
   gStyle->SetOptFit(1111);
   gStyle->SetOptTitle(0);
   gStyle->SetPadColor(kWhite);
   gStyle->SetStatY(0.90);
   gStyle->SetStatX(0.90);
   gStyle->SetStatW(0.15);
   gStyle->SetStatH(0.2);
   gStyle->SetLabelSize(0.04,"X");
   gStyle->SetLabelFont(62,"X");
   gStyle->SetTitleSize(0.04,"X");
   gStyle->SetTitleFont(62,"X");
   gStyle->SetTitleOffset(0.85,"X");

  //  gStyle->SetLabelOffset(0.015,"Y");
   gStyle->SetLabelSize(0.03,"Y");
   gStyle->SetLabelFont(62,"Y");
   gStyle->SetTitleSize(0.04,"Y");
   gStyle->SetTitleFont(62,"Y");
   gStyle->SetTitleOffset(1.3,"Y");
   gStyle->SetTitleX(0.22);
   gStyle->SetTitleY(0.98);
   gStyle->SetTitleSize(0.04,"t");
   gStyle->SetTitleBorderSize(0);
   gStyle->SetCanvasBorderSize(0);

   //////////////////// PLOTS //////
   //// each plot goes in a canvas
   /// most 1-D histograms are grouped into TStack
   /// each plot is stored as pdf, eps png and also in .C (sometimes helps when one cannot re-run the entire code but needs to change any style)
   /// some legends are misplaced, need to check one-by-one
   
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
  
  trueEnu_CCQE_00_02-> SetFillColor(kRed);
  trueEnu_CCMEC_00_02-> SetFillColor(kGreen);
  trueEnu_CCRES_00_02-> SetFillColor(kYellow);
  trueEnu_CCDIS_00_02-> SetFillColor(kBlue);
  trueEnu_other_00_02-> SetFillColor(kMagenta);
 
  THStack *trueEnu_00_02 = new THStack("trueEnu_00_02","");
  trueEnu_00_02-> Add(trueEnu_CCQE_00_02);
  trueEnu_00_02-> Add(trueEnu_CCMEC_00_02);
  trueEnu_00_02-> Add(trueEnu_CCRES_00_02);
  trueEnu_00_02-> Add(trueEnu_CCDIS_00_02);
  trueEnu_00_02-> Add(trueEnu_other_00_02);   
  trueEnu_00_02-> Draw("hist");
  trueEnu_00_02->GetXaxis()->SetTitle("E_{#nu} [GeV/c]");
  trueEnu_00_02->GetYaxis()->SetTitle("number of events");
  
  TLegend *l1 = new TLegend(0.5, 0.7, 0.9, 0.9);
  l1 -> AddEntry(trueEnu_CCQE_00_02, "CCQE", "f");
  l1 -> AddEntry(trueEnu_CCMEC_00_02, "CCMEC", "f");
  l1 -> AddEntry(trueEnu_CCRES_00_02, "CCRES", "f");
  l1 -> AddEntry(trueEnu_CCDIS_00_02, "CCDIS", "f");
  l1 -> Draw();
  
  c1->Update();
  c1->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueEnu_00_02.pdf");
  c1->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueEnu_00_02.eps");
  c1->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueEnu_00_02.png");
  c1->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueEnu_00_02.C");

  TCanvas *c2 = new TCanvas("c2", "c2", 900, 900);
  
  trueNproton_CC_00_02-> SetFillColor(kRed);
  trueNproton_CC_00_02-> Draw("hist");
  trueNproton_CC_00_02->GetXaxis()->SetTitle("proton multiplicity");
  trueNproton_CC_00_02->GetYaxis()->SetTitle("number of events");

  TLegend *l2 = new TLegend(0.5, 0.7, 0.9, 0.9);
  l2 -> AddEntry(trueNproton_CC_00_02, "protons above 35MeV KE", "f");
  l2 -> Draw();

  c2->Update();
  c2->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNproton_CC_00_02.pdf");
  c2->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNproton_CC_00_02.eps");
  c2->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNproton_CC_00_02.png");
  c2->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNproton_CC_00_02.C");  

  TCanvas *c3 = new TCanvas("c3", "c3", 900, 900);
  
  trueNpion_CC_00_02-> SetFillColor(kGreen);
  trueNpion_CC_00_02-> Draw("hist");
  trueNpion_CC_00_02->GetXaxis()->SetTitle("charged pion multiplicity");
  trueNpion_CC_00_02->GetYaxis()->SetTitle("number of events");

  TLegend *l3 = new TLegend(0.5, 0.7, 0.9, 0.9);
  l3 -> AddEntry(trueNpion_CC_00_02, "pions above 10MeV KE", "f");
  l3 -> Draw();

  c3->Update();
  c3->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNpion_CC_00_02.pdf");
  c3->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNpion_CC_00_02.eps");
  c3->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNpion_CC_00_02.png");
  c3->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNpion_CC_00_02.C");  

  TCanvas *c4 = new TCanvas("c4", "c4", 900, 900);
  
  trueNpi0_CC_00_02-> SetFillColor(kCyan);
  trueNpi0_CC_00_02-> Draw("hist");
  trueNpi0_CC_00_02->GetXaxis()->SetTitle("neutral pion multiplicity");
  trueNpi0_CC_00_02->GetYaxis()->SetTitle("number of events");

  TLegend *l4 = new TLegend(0.5, 0.7, 0.9, 0.9);
  l4 -> AddEntry(trueNpi0_CC_00_02, "pions above 10MeV KE", "f");
  l4 -> Draw();

  c4->Update();
  c4->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNpi0_CC_00_02.pdf");
  c4->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNpi0_CC_00_02.eps");
  c4->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNpi0_CC_00_02.png");
  c4->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueNpi0_CC_00_02.C");  
   
  TCanvas *c5 = new TCanvas("c5", "c5", 900, 900);
  
  trueOmega_CC_00_02->SetTitle("");
  trueOmega_CC_00_02->GetXaxis()->SetTitle("true energy transfer, #omega [GeV]");
  trueOmega_CC_00_02->GetYaxis()->SetTitle("number of events");
  trueOmega_CC_00_02->SetLineColor(kRed);
  trueOmega_CC_00_02->SetLineWidth(2);
  trueOmega_CC_00_02->SetMarkerColor(kRed);
  trueOmega_CC_00_02->SetMarkerStyle(20);
  trueOmega_CC_00_02->Draw("e1");
  
  c5->Update();
  c5->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueOmega_CC_00_02.pdf");
  c5->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueOmega_CC_00_02.eps");
  c5->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueOmega_CC_00_02.png");
  c5->Print("/Users/castillofernr/Documents/GiBUU2021/SBND/analysis/plots_00_02/trueOmega_CC_00_02.C");

}
