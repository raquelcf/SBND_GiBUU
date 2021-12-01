//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 30 11:09:02 2021 by ROOT version 6.24/06
// from TTree RootTuple/RootTuple
// found on file: ../simCC_oaa_00_02/simCC_oaa_00_02.root
//////////////////////////////////////////////////////////

#ifndef GiBuuCC_oaa_00_02_h
#define GiBuuCC_oaa_00_02_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class GiBuuCC_oaa_00_02 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        weight;
   vector<int>     *barcode;
   vector<double>  *Px;
   vector<double>  *Py;
   vector<double>  *Pz;
   vector<double>  *E;
   Int_t           evType;
   Double_t        lepIn_E;
   Double_t        lepIn_Px;
   Double_t        lepIn_Py;
   Double_t        lepIn_Pz;
   Double_t        lepOut_E;
   Double_t        lepOut_Px;
   Double_t        lepOut_Py;
   Double_t        lepOut_Pz;
   Double_t        nuc_E;
   Double_t        nuc_Px;
   Double_t        nuc_Py;
   Double_t        nuc_Pz;

   // List of branches
   TBranch        *b_weight;   //!
   TBranch        *b_barcode;   //!
   TBranch        *b_Px;   //!
   TBranch        *b_Py;   //!
   TBranch        *b_Pz;   //!
   TBranch        *b_E;   //!
   TBranch        *b_evType;   //!
   TBranch        *b_lepIn_E;   //!
   TBranch        *b_lepIn_Px;   //!
   TBranch        *b_lepIn_Py;   //!
   TBranch        *b_lepIn_Pz;   //!
   TBranch        *b_lepOut_E;   //!
   TBranch        *b_lepOut_Px;   //!
   TBranch        *b_lepOut_Py;   //!
   TBranch        *b_lepOut_Pz;   //!
   TBranch        *b_nuc_E;   //!
   TBranch        *b_nuc_Px;   //!
   TBranch        *b_nuc_Py;   //!
   TBranch        *b_nuc_Pz;   //!

   GiBuuCC_oaa_00_02(TTree *tree=0);
   virtual ~GiBuuCC_oaa_00_02();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef GiBuuCC_oaa_00_02_cxx
GiBuuCC_oaa_00_02::GiBuuCC_oaa_00_02(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../simCC_oaa_0.0_0.2/simCC_oaa_0.0_0.2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../simCC_oaa_0.0_0.2/simCC_oaa_0.0_0.2.root");
      }
      f->GetObject("RootTuple",tree);

   }
   Init(tree);
}

GiBuuCC_oaa_00_02::~GiBuuCC_oaa_00_02()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GiBuuCC_oaa_00_02::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GiBuuCC_oaa_00_02::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void GiBuuCC_oaa_00_02::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   barcode = 0;
   Px = 0;
   Py = 0;
   Pz = 0;
   E = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("barcode", &barcode, &b_barcode);
   fChain->SetBranchAddress("Px", &Px, &b_Px);
   fChain->SetBranchAddress("Py", &Py, &b_Py);
   fChain->SetBranchAddress("Pz", &Pz, &b_Pz);
   fChain->SetBranchAddress("E", &E, &b_E);
   fChain->SetBranchAddress("evType", &evType, &b_evType);
   fChain->SetBranchAddress("lepIn_E", &lepIn_E, &b_lepIn_E);
   fChain->SetBranchAddress("lepIn_Px", &lepIn_Px, &b_lepIn_Px);
   fChain->SetBranchAddress("lepIn_Py", &lepIn_Py, &b_lepIn_Py);
   fChain->SetBranchAddress("lepIn_Pz", &lepIn_Pz, &b_lepIn_Pz);
   fChain->SetBranchAddress("lepOut_E", &lepOut_E, &b_lepOut_E);
   fChain->SetBranchAddress("lepOut_Px", &lepOut_Px, &b_lepOut_Px);
   fChain->SetBranchAddress("lepOut_Py", &lepOut_Py, &b_lepOut_Py);
   fChain->SetBranchAddress("lepOut_Pz", &lepOut_Pz, &b_lepOut_Pz);
   fChain->SetBranchAddress("nuc_E", &nuc_E, &b_nuc_E);
   fChain->SetBranchAddress("nuc_Px", &nuc_Px, &b_nuc_Px);
   fChain->SetBranchAddress("nuc_Py", &nuc_Py, &b_nuc_Py);
   fChain->SetBranchAddress("nuc_Pz", &nuc_Pz, &b_nuc_Pz);
   Notify();
}

Bool_t GiBuuCC_oaa_00_02::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GiBuuCC_oaa_00_02::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t GiBuuCC_oaa_00_02::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef GiBuuCC_oaa_00_02_cxx
