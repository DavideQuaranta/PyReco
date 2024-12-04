#define apv_raw_cxx
// The class definition in apv_raw.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("apv_raw.C")
// root> T->Process("apv_raw.C","some options")
// root> T->Process("apv_raw.C+")
//

#include "apv_raw.h"
#include <TH2.h>
#include <TStyle.h>

void apv_raw::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   // TString outfilename("run" + run + "_timefit.root");
   // out_file = new TFile(outfilename, "RECREATE");

   fOutputTree = new TTree("output_tree", "apv_raw_with_fitted_times");
   fOutputTree->Branch("fitted_time", &fitted_times);
   fOutputTree->Branch("fitted_time_err", &fitted_times_err);

   signalhisto = new TH1F("signal_histo", "signal_histo", 27, 0, 675);
   // times = new TH1F("times", "times", 200, -200, 675);
   // fermidirac = new TF1("fermidirac", "([0]/(1+exp(-(x-[1])/[2])))+[3]", 0, 675);
}

void apv_raw::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

bool apv_raw::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

   fitted_times.clear();
   fitted_times_err.clear();

   fProcessedEntries++;

   // commment this if to process everything
   // if (fProcessedEntries >= 1000)
   // {
   //    // out_file->cd();
   //    // out_file->WriteTObject(fOutputTree);
   //    // out_file->Close();
   //    return kFALSE; // Return false to stop further processing
   // }

   // Example: Fill signalhisto with data
   for (size_t i = 0; i < mmStrip.GetSize(); ++i)
   {
      for (size_t l = 0; l < raw_q[i].size(); ++l)
      {
         signalhisto->SetBinContent(l, raw_q[i][l]);
         signalhisto->SetBinError(l, 10);
      }

      // Taken from TBReco MMHit, all rights reserved
      // Example: Fit signalhisto with fermidirac function
      double maxval = signalhisto->GetMaximum();
      // double rms= signalhisto->GetRMS() ;
      int maxBin = signalhisto->GetMaximumBin();
      double timeAtMaxBin = signalhisto->GetBinCenter(maxBin);
      double timeAtMaxBinPlus1 = signalhisto->GetBinCenter(maxBin + 1);
      // checking baseline values
      double baselineMin(999);
      double baselineMax(-999);
      double baseEstimate[4];
      double baselineMean(0);

      baseEstimate[0] = signalhisto->GetBinContent(1);
      baseEstimate[1] = signalhisto->GetBinContent(signalhisto->GetNbinsX() - 2);
      baseEstimate[2] = signalhisto->GetBinContent(signalhisto->GetNbinsX() - 1);
      baseEstimate[3] = signalhisto->GetBinContent(signalhisto->GetNbinsX());

      for (int k = 0; k < 4; k++)
      {
         if (baseEstimate[k] < baselineMin)
            baselineMin = baseEstimate[k];
         if (baseEstimate[k] > baselineMax)
            baselineMax = baseEstimate[k];
         baselineMean += baseEstimate[k];
      }
      baselineMean = baselineMean / 4.;

      fermidirac->SetParameters(maxval, timeAtMaxBin / 2., timeAtMaxBin / 4.5, (baselineMin + baselineMax) / 2.);
      fermidirac->SetParLimits(1, -500, 2000);
      fermidirac->SetParLimits(3, baselineMin, baselineMax);

      if (timeAtMaxBinPlus1 < 120.)
         fermidirac->FixParameter(3, baselineMean);

      // // DQ: if the time is too short the baseline its affected by the peak
      // if (timeAtMaxBinPlus1<120.) {
      //    double baselineMean2 = (baseEstimate[1]+baseEstimate[2]+baseEstimate[3])/3.;
      //    fermidirac->FixParameter(3,baselineMean2);
      // }
      // // DQ:end of my addition to TBReco

      signalhisto->Fit(fermidirac, "Q", " ", 0., timeAtMaxBinPlus1);

      // out_file->WriteTObject(signalhisto);

      // TCanvas histo_canvas("canvas", "canvas", 100,100, 800, 600);
      // signalhisto->Draw();
      // histo_canvas.Draw();

      // Store results or perform further analysis as needed
      // Example: Store fit parameters in a TTree or print to console
      double fit_time = fermidirac->GetParameter(1);
      double fit_time_err = fermidirac->GetParError(1);

      fitted_times.push_back(fit_time);
      fitted_times_err.push_back(fit_time_err);
   }

   fOutputTree->Fill();
   // std::cout << "Entry: " << entry << "Event : " << *evt << std::endl;
   if (entry % 100 == 0)
      printf("Event: %llu\n", *evt);

   return true;
}

void apv_raw::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void apv_raw::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   // fOutputTree->SetDirectory(out_file);
   // fOutputTree->Write();
   // out_file->cd();
   out_file->WriteTObject(fOutputTree);
   out_file->Close();

   // TCanvas *c1 = new TCanvas("", "", 100,100,800,600);
   // times->Draw();
   // c1->Draw();
}