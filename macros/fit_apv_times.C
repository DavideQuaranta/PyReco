// process_files.C

#include <iostream>
#include <string>
#include <vector>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
// #include "apv_raw.h"
#include <vector>

void fit_apv_times()
{

    // Specify the directory containing your ROOT files
    // TString directory("../TestBeamNovember/raw/run"); // Replace with your directory path
    // std::vector<TString> runs = {"2082", "2116", "2117", "2121", "2124", "2125", "2126", "2127", "2128", "2129", "2130", "2131", "2132", "2133", "2134", "2135", "2136", "2137"};
    std::vector<TString> directory = {"../data/raw/run"};
    std::vector<TString> runs = {"2124"};
    TString extension("_xtalk1_xtalk2.root");

    for (int i = 0; i < runs.size(); i++)
    {
        // Create a TChain
        TChain chain("apv_raw"); // Replace "treeName" with the actual name of the tree in your files
        std::cout << "analyzing run" << runs[i].Data() << std::endl;
        TString inputFile(directory[i] + runs[i] + extension);
        chain.Add(inputFile);

        // Create the output file and tree names
        TString outputFileName("../data/timefit/run" + runs[i] + "_timefit.root");
        TFile *outputFile = TFile::Open(outputFileName, "RECREATE");
        // TTree *outputTree = new TTree("output_tree", "output_tree");

        // Create an instance of the selector
        apv_raw *selector = new apv_raw();
        selector->out_file = outputFile;
        // selector->fOutputTree = outputTree;

        // Process the chain with the selector
        chain.Process(selector);

        // Clean up
        delete selector;
    }
}
