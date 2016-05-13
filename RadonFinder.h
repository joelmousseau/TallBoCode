//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jan  7 12:13:45 2016 by ROOT version 5.34/32
// from TChain wavedata/
//////////////////////////////////////////////////////////

#ifndef RadonFinder_h
#define RadonFinder_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

static const Int_t N_SAMPLES = 50000;
static const Double_t TICK_TO_NS = 8.0;

class RadonFinder {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           WaveNo;
   Int_t           PulseNo;
   Int_t           StartTick;
   Int_t           EndTick;
   Long64_t        Date;
   Double_t        Amplitude;
   Double_t        Area;

   // List of branches
   TBranch        *b_WaveNo;   //!
   TBranch        *b_PulseNo;   //!
   TBranch        *b_StartTick;   //!
   TBranch        *b_EndTick;   //!
   TBranch        *b_Date;   //!
   TBranch        *b_Amplitude;   //!
   TBranch        *b_Area;   //!

   RadonFinder(TTree *tree=0);
   virtual ~RadonFinder();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     MakeHists(std::string fillRate = "Unkown");
   virtual void     MakePlots(std::string fillRate = "Unkown");
   virtual void     MakeComparisonPlots();
   virtual Bool_t  Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual double    CalcLogError(double val, double err);
   virtual void     SetRedHeatPalette();
};

#endif

#ifdef RadonFinder_cxx
RadonFinder::RadonFinder(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("wavedata",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("wavedata","");
      chain->Add("FastFill2015NewPulse.root/wavedata");
      chain->Add("FastFill2016NewPulse.root/wavedata");
      //chain->Add("SlowFill50mvDiv.root/wavedata");
      //chain->Add("SlowFill10mvDiv.root/wavedata");
      //chain->Add("SlowFillAfterSplit.root/wavedata");
      //chain->Add("SlowFillAfterSplitNePulse.root/wavedata");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

RadonFinder::~RadonFinder()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RadonFinder::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RadonFinder::LoadTree(Long64_t entry)
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

void RadonFinder::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("WaveNo", &WaveNo, &b_WaveNo);
   fChain->SetBranchAddress("PulseNo", &PulseNo, &b_PulseNo);
   fChain->SetBranchAddress("StartTick", &StartTick, &b_StartTick);
   fChain->SetBranchAddress("EndTick", &EndTick, &b_EndTick);
   fChain->SetBranchAddress("Date", &Date, &b_Date);
   fChain->SetBranchAddress("Amplitude", &Amplitude, &b_Amplitude);
   fChain->SetBranchAddress("Area", &Area, &b_Area);
   Notify();
}

Bool_t RadonFinder::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RadonFinder::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RadonFinder::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RadonFinder_cxx
