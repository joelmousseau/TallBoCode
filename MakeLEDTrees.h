//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 18 11:21:46 2015 by ROOT version 5.34/32
// from TChain waveformdata/
//////////////////////////////////////////////////////////

#ifndef MakeLEDTrees_h
#define MakeLEDTrees_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

static const Int_t N_SAMPLES = 2500;
static const Double_t TICK_TO_NS = 0.8;

class MakeLEDTrees {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   
   

   // Declaration of leaf types
   Double_t        ch1wfms[N_SAMPLES];
   Double_t        ch2wfms[N_SAMPLES];
   Double_t        ch3wfms[N_SAMPLES];
   Double_t        ch4wfms[N_SAMPLES];

   // List of branches
   TBranch        *b_ch1wfms;   //!
   TBranch        *b_ch2wfms;   //!
   TBranch        *b_ch3wfms;   //!
   TBranch        *b_ch4wfms;   //!

   MakeLEDTrees(TTree *tree=0);
   virtual ~MakeLEDTrees();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual void     ViewWave();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual double    getAmplitude (std::vector<double> wave, int startTick, int endTick);
   virtual int      FindPulseEndTime(std::vector<double> ADCs, int MPulseStartTime, double Pedestal);
   virtual int      FindPulseStartTime(std::vector<double> ADCs, int MPulseStartTime, double Threshold);
   virtual double    getArea(std::vector<double> wave, int startTick, int endTick);
   virtual std::string getFileName();
};

#endif

#ifdef MakeLEDTrees_cxx
MakeLEDTrees::MakeLEDTrees(TTree *tree) : fChain(0) 
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
      f->GetObject("waveformdata",tree);


#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("waveformdata","");
            chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-19_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-19_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-20_00_02.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-20_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-20_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-20_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-21_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-21_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-21_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-21_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-22_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-22_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-22_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-22_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-23_00_02.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-23_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-23_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-27-23_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-00_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-00_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-00_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-00_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-01_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-01_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-01_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-01_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-02_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-02_15_02.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-02_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-02_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-03_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-03_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-03_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-03_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-04_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-04_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-04_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-04_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-05_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-05_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-05_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-05_45_02.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-06_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-06_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-06_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-06_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-07_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-07_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-07_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-07_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-08_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-08_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-08_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-08_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-09_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-09_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-09_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-09_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-10_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-10_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-10_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-10_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-11_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-11_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-11_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-11_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-12_00_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-12_15_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-12_30_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-12_45_01.root/waveformdata");
      chain->Add("/uboone/data/users/joelam/RadonTests/TallBo_Rn_2016-01-28-13_00_01.root/waveformdata");
      
    
      
      
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

MakeLEDTrees::~MakeLEDTrees()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MakeLEDTrees::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MakeLEDTrees::LoadTree(Long64_t entry)
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

std::string MakeLEDTrees::getFileName(){
    return ( fChain->GetCurrentFile() )->GetName();
}

void MakeLEDTrees::Init(TTree *tree)
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

   fChain->SetBranchAddress("ch1wfms", ch1wfms, &b_ch1wfms);
   fChain->SetBranchAddress("ch2wfms", ch2wfms, &b_ch2wfms);
   fChain->SetBranchAddress("ch3wfms", ch3wfms, &b_ch3wfms);
   fChain->SetBranchAddress("ch4wfms", ch4wfms, &b_ch4wfms);
   Notify();
}

Bool_t MakeLEDTrees::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MakeLEDTrees::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MakeLEDTrees::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MakeLEDTrees_cxx
