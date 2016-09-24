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
#include "AnitaGeomTool.h"

#include "scopeParser.h"

#include "FFTtools.h"
#include "brotterTools.h"

using namespace std;



TGraph* scopeParseAndAverage(int run) {

  //set where the antenna response scope waveform is
  stringstream inputFile,inputDir;

  inputDir << "/Volumes/ANITA3Data/antarctica14/scopeTraces/";
  inputDir << "run" << run << "/";
    

  //set up for correlating and averaging
  const int numEntries = 248;
  TGraph* waveGraphs[numEntries];


  string dataLine;

  //get all the graphs out of the scope readout 
  for (int entry=0; entry<numEntries; entry++) {  
    inputFile.str("");
    inputFile << inputDir.str() << "wave_" << entry+1;
    
    waveGraphs[entry] = getScopeGraph(inputFile.str(),0);
    
  }
  cout << "scopeParseAndAverage():got all the graphs I wanted" << endl;

  //correlate and average graphs
  TGraph *avgWaveGraph = FFTtools::correlateAndAverage(numEntries,waveGraphs);
  avgWaveGraph->SetName("avgAntGraph");


  //memory management!
  for (int entry=0; entry<numEntries; entry++) {
    delete(waveGraphs[entry]);
  }
      
  cout << avgWaveGraph->GetN() << endl;
  return avgWaveGraph;
  
}


void getRingPolPhiFromName( string name, AnitaRing::AnitaRing_t &ring, AnitaPol::AnitaPol_t &pol, int &phi ) {

  char num[2];
  num[0] = name[0];
  num[1] = name[1];

  phi = atoi( num );

  if (name[2] == 'M' || name[2] =='m') {
    ring = AnitaRing::kMiddleRing; }
  else if (name[2] == 'B' || name[2] == 'b') {
    ring = AnitaRing::kBottomRing; }
  else if (name[2] == 'T' || name[2] == 't') {
    ring = AnitaRing::kTopRing; }
  else {
    cout << "Couldn't figure out which ring you meant!" << endl;
    return; }
  
  if (name[3] =='V' || name[3] == 'v') {
    pol = AnitaPol::kVertical; }
  else if (name[3] =='H' || name[3] =='h') {
    pol = AnitaPol::kHorizontal; }
  else {
    cout << "Couldn't figure out which polarization you meant!" << endl;
    return; }

  return;
}


int findCalRun(string name) {

  int runNum = -1;
  ifstream runLogFile;
  runLogFile.open("runLog66dB.txt");
  if (!runLogFile) cout << "didn't open file?" << endl;

  int pos = 0;
  string subLine;
  string nameStr;
    

  string line;
  int found = 0;
  while ( getline(runLogFile,line) ) {
    if (line.front() == '#') continue;

    pos = line.find(name);
    runNum = stoi(line.substr(0,pos));

    if (pos != -1) break;
  }

  cout << name << " associated " << nameStr << " with run " << runNum << endl;
  return runNum;
}



TGraph* surfParseAndAverage(string antName) {
  /*
    Raw: 256ish samples @2.6GS/s (.0385ns)
    Plan: Upsample to 10GS/s (0.1ns), zero pad to 1024 samples
  */

  
  //for exiting after an error
  TGraph *nullGr = new TGraph();

  int calRun = -1;
  calRun = findCalRun(antName);
    


  //figure out what that antName string means in useful types
  AnitaPol::AnitaPol_t pol = AnitaPol::kNotAPol;
  AnitaRing::AnitaRing_t ring = AnitaRing::kNotARing;
  int phi = -1;
  getRingPolPhiFromName(antName,ring,pol,phi);
  if (phi == -1 || pol == AnitaPol::kNotAPol || ring == AnitaRing::kNotARing) {
    cout << "Didn't unwrap name string correctly :( I quit!" << endl;
    return nullGr; }

  //I need to know the channel number too for determining which labrador
  AnitaGeomTool *geom = new AnitaGeomTool();
  int chanIndex = geom->getChanIndexFromRingPhiPol(ring,phi-1,pol);
  delete geom;

  cout << "phi:" << phi << " ring:" << AnitaRing::ringAsChar(ring) << " pol:" << AnitaPol::polAsChar(pol) << " chanIndex:" << chanIndex << endl;

  //set up where the data is
  stringstream name;
  stringstream dirName;
  dirName << "/Volumes/ANITA3Data/antarctica14/root/";

//  const char* anitaDataDir= getenv("ANITA3_CALDATA");
//  if (anitaDataDir == NULL) {
//    cout << "ANITA3_CALDATA bash environment variable not set!  I don't know where to look for data :(" << endl;
//    return nullGr; }
//  dirName << anitaDataDir << "/root/";
//  //  dirName << "/storagec/antarctica14/root/";
  //TEST - REMOVE LATER (I don't have any runs to look at here on this plane :'( )
 
  

  //Importing Events!
  TChain *eventTree = new TChain("eventTree","eventTree");
  name.str("");
  name << dirName.str() << "run" << calRun << "/eventFile" << calRun << ".root";
  eventTree->Add(name.str().c_str());
  if (eventTree->GetEntries() == 0) {
    cout << "You messed up and there aren't any events, I quit." << endl;
    return nullGr; }  
  //set address
  RawAnitaEvent *event = NULL;
  eventTree->SetBranchAddress("event",&event);


  //Importing Headers!
  TChain *headTree = new TChain("headTree","headTree");
  name.str("");
  name << dirName.str() << "run" << calRun << "/headFile" << calRun << ".root";
  headTree->Add(name.str().c_str());
  if (eventTree->GetEntries() == 0) {
    cout << "You messed up and there aren't any headers, I quit." << endl;
    return nullGr; }
  //set address
  RawAnitaHeader *head = NULL;
  headTree->SetBranchAddress("header",&head);


  //Just check that they have the same number of entires for fun
  int lenEntries = -1;
  if (headTree->GetEntries() == eventTree->GetEntries()) {
    lenEntries = headTree->GetEntries();
    cout << "Number of entries the same for both eventTree and headTree :D (" << lenEntries << ")" << endl; }
  else {
    lenEntries = headTree->GetEntries();
    cout << "headTree and eventTree have a different number of entries..." << endl;
    cout << "headTree:" << headTree->GetEntries() << " eventTree:" << eventTree->GetEntries();  }

  
  //define how many waveforms we are going to average together
  const int numEntries = 1024;

  //set up the pointer array for all those TGraphs
  TGraph *waveGraphs[numEntries];

  TGraph *averagedFFT = new TGraph();
  averagedFFT->SetName("averagedFFT");

  //**Loop Through Those Entries!!!!**
  int entryAveraged = 0;
  for (int entry=0; entry<headTree->GetEntries(); entry++) {
    //just in case there aren't any entries, I'd rather not segfault
    if (entryAveraged == numEntries) break;
    //move through entries
    cout << entryAveraged << "(" << entry << ") /" << numEntries << "\r";
    fflush(stdout);
    headTree->GetEntry(entry);
    eventTree->GetEntry(entry);

    
    //I only want to look at one lab chip since they can all be different!
    // continue will make it skip everything after this point if it isn't satisfied
    if (event->getLabChip(chanIndex) != 0) continue;

    //get calibrated event
    UsefulAnitaEvent *useful = new UsefulAnitaEvent(event,WaveCalType::kOnlyTiming,head);

    //Get calibrated waveform for the channel of interest
    TGraph *calibGraph = useful->getGraph(ring,phi-1,pol);
    delete(useful);
    //upsample to 10GS/s (0.1ns bins)
    TGraph *interpGraph = FFTtools::getInterpolatedGraph(calibGraph,0.1);
    delete(calibGraph);
    //zero pad to 1024 bins
    TGraph *finalGraph = brotterTools::zeroPadToLength(interpGraph,1024);
    delete(interpGraph);


    //store that in the array to average together later
    waveGraphs[entryAveraged] = new TGraph(*finalGraph);
    delete(finalGraph);
    
    //set the names and titles of those graphs for fun (also I can save them later if I want)
    name.str("");
    name << "wave_" << antName << "_" << entryAveraged;
    waveGraphs[entryAveraged]->SetName(name.str().c_str());

    name.str("");
    name << "ev" << head->eventNumber;
    waveGraphs[entryAveraged]->SetTitle(name.str().c_str());


    //incriment the number of graphs averaged, since a bunch of events aren't valid since they are the wrong lab
    entryAveraged++;
  }
  cout << "Saved " << entryAveraged << " waveforms to memory, next step: averaging" << endl;

//  //Correlate and Average all those graphs together
//  TGraph *avgWaveGraph = new TGraph(*waveGraphs[entryAveraged-1]);
  TGraph *avgWaveGraph = FFTtools::correlateAndAverage(entryAveraged,waveGraphs);
  name.str("");
  name << "avgSurfGraph";
  avgWaveGraph->SetName(name.str().c_str());
  avgWaveGraph->SetLineColor(kRed);
  avgWaveGraph->SetMarkerColor(kRed);
  




  // Save a bunch of interesting things about generating the surf waveform (since it messes up so often)
  name.str("");
  name << "waveforms/" << antName << "_surfInfo.root";
  TFile *testFile = TFile::Open(name.str().c_str(),"recreate");

  //Save a copy of the histogram
  TH1D *hist = brotterTools::correlationDistribution(entryAveraged,waveGraphs);
  TGraph *corrPattern = brotterTools::correlationPattern(entryAveraged,waveGraphs);
  corrPattern->SetName("corrPattern1");
  corrPattern->SetTitle("corrPattern1");

  TGraph *corrPattern2 = brotterTools::correlationPattern(entryAveraged,waveGraphs);
  corrPattern2->SetName("corrPattern2");
  corrPattern2->SetTitle("corrPattern2");

  TGraph *corrPattern3 = brotterTools::correlationPattern(entryAveraged,waveGraphs);
  corrPattern3->SetName("corrPattern3");
  corrPattern3->SetTitle("corrPattern3");

  TGraph *corrPattern4 = brotterTools::correlationPattern(entryAveraged,waveGraphs);
  corrPattern4->SetName("corrPattern4");
  corrPattern4->SetTitle("corrPattern4");

  //write
  for (int i=2; i<entryAveraged; i++) {
    cout << i << endl;
    TGraph *tempGraph = FFTtools::correlateAndAverage(i,waveGraphs);
    name.str("");
    name << "averagedGraph" << i;
    tempGraph->SetName(name.str().c_str());
    tempGraph->Write();
    delete(tempGraph);
  }
  hist->Write();
  corrPattern->Write();
  corrPattern2->Write();
  corrPattern3->Write();
  corrPattern4->Write();
  //  averagedFFT->Write();

  testFile->Close();

  
  //memory management!
  for (int i=0; i<entryAveraged; i++) {
    delete(waveGraphs[i]); 
  }

  
  cout << "done with making surf averaged graph" << endl;
  //Done!
  return avgWaveGraph;

}



int main(int argc, char** argv) {
  
  cout << "Starting!  :D" << endl;

  string antName;

  if (argc == 1) {
    cout << "Using default channel of 09BH" << endl;
    antName = "09BH";
  }
  else if (argc == 2) {
    cout << "Using " << argv[1] << " as channel" << endl;
    antName = argv[1];
  }
  else {
    cout << "Usage: " << argv[0] << " channelName (eg:09BH)" << endl;
  }

  int calRun =findCalRun(antName);

  //Get the two pulses
  TGraph *waveGraphScope = scopeParseAndAverage(calRun);
  cout << "Got Scope Waveforms!" << endl;
  TGraph *waveGraphSurf = surfParseAndAverage(antName);
  cout << "Got Surf Waveforms!" << endl;

  //save them to a .txt file
  ofstream outStream;
  stringstream fileName;
  //Surf
  fileName.str("");
  fileName << "waveforms/" << antName << "_avgSurfWaveform.txt";
  outStream.open(fileName.str());
  for (int pt=0; pt<waveGraphSurf->GetN(); pt++) {
    outStream << waveGraphSurf->GetX()[pt] << " " << waveGraphSurf->GetY()[pt] << endl;
  }
  outStream.close();

  //Scope 
  fileName.str("");
  fileName << "waveforms/" << antName << "_avgScopeWaveform.txt";
  outStream.open(fileName.str());
  for (int pt=0; pt<waveGraphScope->GetN(); pt++) {
    outStream << waveGraphScope->GetX()[pt] << " " << waveGraphScope->GetY()[pt] << endl;
  }
  outStream.close();

  //save them to a .root file
  fileName.str("");
  fileName << "waveforms/" << antName << "_avgWaveforms.root";
  TFile *rootFile = TFile::Open(fileName.str().c_str(),"recreate");
  waveGraphSurf->Write();
  waveGraphScope->Write();

  cout << "Done! :)" << endl;

  

  return -1;
}

