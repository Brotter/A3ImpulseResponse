#include <iostream>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
//root
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TProfile2D.h"
//anita
#include "RawAnitaHeader.h"
#include "RawAnitaEvent.h"
#include "UsefulAnitaEvent.h"
#include "AnitaConventions.h"

#include "antScopeParser.h"

#include "FFTtools.h"
#include "FFTWComplex.h"

using namespace std;



TGraph* antennaParseAndAverage() {

  //set where the antenna response scope waveform is
  stringstream inputFile;
  inputFile << "/Volumes/ANITA3Data/Rooftop_Seevey_Antenna/Antenna_Impulse_Testing/";
  inputFile << "Attempt_3_23July2012/NewSeavey/Vpol/VPulse_BoresightToBoresight3_10dBatt_23July2012.dat";

  ifstream inData;
  inData.open(inputFile.str());

  //get the first comment out of the way
  string dataLine;
  getline(inData,dataLine);

  //set up for correlating and averaging
  const int numEntries = 248;
  TGraph* waveGraphs[numEntries];

  //get all the graphs out of the scope readout 
  for (int entry=0; entry<numEntries; entry++) {
    waveGraphs[entry] = getScopeGraph(inData);
  }
    
  inData.close();

  //correlate and average graphs
  TGraph *avgWaveGraph = FFTtools::correlateAndAverage(numEntries,waveGraphs);
  avgWaveGraph->SetName("avgAntGraph");


  //memory management!
  for (int entry=0; entry<numEntries; entry++) {
    delete(waveGraphs[entry]);
  }
      
  return avgWaveGraph;
  
}


TGraph* FIDParseAndAverage() {

  //set where the antenna response scope waveform is
  stringstream inputFile;
  inputFile << "/Volumes/ANITA3Data/Rooftop_Seevey_Antenna/Antenna_Impulse_Testing/";
  inputFile << "Attempt_1_28June2012/Cables/Orange_20dB_30dB_20dB_2.dat";

  ifstream inData;
  inData.open(inputFile.str());

  //get the first comment out of the way
  string dataLine;
  getline(inData,dataLine);

  //set up for correlating and averaging
  const int numEntries = 248;
  TGraph* waveGraphs[numEntries];

  //get all the graphs out of the scope readout 
  for (int entry=0; entry<numEntries; entry++) {
    waveGraphs[entry] = getScopeGraph(inData);
  }
    
  inData.close();

  //correlate and average graphs
  TGraph *avgWaveGraph = FFTtools::correlateAndAverage(numEntries,waveGraphs);
  avgWaveGraph->SetName("avgFIDGraph");


  //memory management!
  for (int entry=0; entry<numEntries; entry++) {
    delete(waveGraphs[entry]);
  }
      
  return avgWaveGraph;
  
}



int main(int argc, char** argv) {
  
  cout << "Starting!  :D" << endl;

  //Get the two pulses
  TGraph *waveGraphFID = FIDParseAndAverage();
  TGraph *waveGraphAnt = antennaParseAndAverage();

  //save them to a .txt file
  ofstream outStream;
  //FID
  outStream.open("avgFID.txt");
  for (int pt=0; pt<waveGraphFID->GetN(); pt++) {
    outStream << waveGraphFID->GetX()[pt] << " " << waveGraphFID->GetY()[pt] << endl;
  }
  outStream.close();

  //antenna
  outStream.open("avgAnt.txt");
  for (int pt=0; pt<waveGraphAnt->GetN(); pt++) {
    outStream << waveGraphAnt->GetX()[pt] << " " << waveGraphAnt->GetY()[pt] << endl;
  }
  outStream.close();

  //save them to a .root file
  TFile *rootFile = TFile::Open("antennaMeasurementAverages.root","recreate");
  waveGraphFID->Write();
  waveGraphAnt->Write();
  rootFile->Close();

  cout << "Done! :)" << endl;

  


  return -1;
}

