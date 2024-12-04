//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jul  6 14:19:33 2024 by ROOT version 6.32.00
// from TTree apv_raw/apv_raw
// found on file: ./run2124_xtalk1_xtalk2.root
//////////////////////////////////////////////////////////

#ifndef apv_raw_h
#define apv_raw_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TCanvas.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TF1.h>

// Headers needed by this particular selector
#include <vector>
#include <string>

class apv_raw : public TSelector
{
public:
   TTreeReader fReader; //! the tree reader
   TTree *fChain = 0;   //! pointer to the analyzed TTree or TChain

   TFile *out_file = NULL;
   TTree *fOutputTree = NULL;

   TString run;
   Int_t fProcessedEntries = 0;

   void SetRunNumber(TString num) { run = num; }

   std::vector<float> fitted_times;
   std::vector<float> fitted_times_err;

   TH1F *signalhisto;
   TH1F *times;
   // TF1 *fermidirac;
   TF1 *fermidirac = new TF1("fermidirac", "([0]/(1+exp(-(x-[1])/[2])))+[3]", 0, 675);

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<ULong64_t> evt = {fReader, "evt"};
   TTreeReaderValue<UInt_t> error = {fReader, "error"};
   TTreeReaderValue<Int_t> daqTimeSec = {fReader, "daqTimeSec"};
   TTreeReaderValue<Int_t> daqTimeMicroSec = {fReader, "daqTimeMicroSec"};
   TTreeReaderValue<Int_t> srsTimeStamp = {fReader, "srsTimeStamp"};
   TTreeReaderValue<UInt_t> srsTrigger = {fReader, "srsTrigger"};
   TTreeReaderArray<unsigned int> srsFec = {fReader, "srsFec"};
   TTreeReaderArray<unsigned int> srsChip = {fReader, "srsChip"};
   TTreeReaderArray<unsigned int> srsChan = {fReader, "srsChan"};
   TTreeReaderArray<std::string> mmChamber = {fReader, "mmChamber"};
   TTreeReaderArray<int> mmLayer = {fReader, "mmLayer"};
   TTreeReaderArray<char> mmReadout = {fReader, "mmReadout"};
   TTreeReaderArray<int> mmStrip = {fReader, "mmStrip"};
   TTreeReaderArray<std::vector<short>> raw_q = {fReader, "raw_q"};

   apv_raw(TTree * /*tree*/ = 0) {}
   ~apv_raw() override {}
   Int_t Version() const override { return 2; }
   void Begin(TTree *tree) override;
   void SlaveBegin(TTree *tree) override;
   void Init(TTree *tree) override;
   bool Notify() override;
   bool Process(Long64_t entry) override;
   Int_t GetEntry(Long64_t entry, Int_t getall = 0) override { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   void SetOption(const char *option) override { fOption = option; }
   void SetObject(TObject *obj) override { fObject = obj; }
   void SetInputList(TList *input) override { fInput = input; }
   TList *GetOutputList() const override { return fOutput; }
   void SlaveTerminate() override;
   void Terminate() override;

   ClassDefOverride(apv_raw, 0);
};

#endif

#ifdef apv_raw_cxx
void apv_raw::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

bool apv_raw::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

#endif // #ifdef apv_raw_cxx
