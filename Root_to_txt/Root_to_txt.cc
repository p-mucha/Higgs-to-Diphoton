//libraries
#include <fstream>
#include "TROOT.h" 
#include "iostream"
#include <iomanip>
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "string.h"
#include "TChain.h"
#include "TTree.h"


// Declare void function
void Root_to_txt(){

  // Input file directory
  string inDir = "./";  
  // Output file directory
  string outDir = "./";

  // Name of the input file
  string inFilehist1 = "Input_File.root";

  // Access data from tree branch: demo/branch_name
  TFile *file = new TFile((inDir + inFilehist1).c_str());
  TH1D *histogram = (TH1D*) file->Get("demo/PHOTONS_lower_range")->Clone();

  // Declare output file in txt format
  std::ofstream output("Output_Histogram.txt");

  // Loop over histogram bins
  // Data will be in format: bin centres, bin entries
  for (int i = 1; i <= histogram->GetNbinsX(); i++) {
    double bin_center = histogram->GetBinCenter(i);
    double bin_entry = histogram->GetBinContent(i);
    output << bin_center << ", " << bin_entry << std::endl;
  }

  output.close();
  file->Close();

}
