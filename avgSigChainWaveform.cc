#include <iostream>
#include <math.h>
#include <string>
#include <sstream>
#include <fstream>
#include <sys/stat.h>

//root
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"

//anita
#include "RawAnitaHeader.h"
#include "RawAnitaEvent.h"
#include "UsefulAnitaEvent.h"
#include "AnitaConventions.h"
#include "AnitaGeomTool.h"
#include "AnitaVersion.h"
#include "AnalysisWaveform.h"
#include "RFInterpolate.h"

#include "scopeParser.h"

#include "FFTtools.h"
#include "brotterTools.h"

using namespace std;


string outDir = "waveforms_new/";


void printDT(TGraph* inGraph,string note = "") {

  cout << "printDT (" << note << "): " <<  inGraph->GetX()[1] - inGraph->GetX()[0] << endl;

  return;

}



TGraph* scopeParseAndAverage(int run,int chan,int numEntries=248) {

  //set where the antenna response scope waveform is
  stringstream inputFile,inputDir;

  char* calDir = getenv("ANITA3_CALDATA");
  inputDir << calDir << "/antarctica14/scopeTraces/";
  inputDir << "run" << run << "/";
    
  struct stat sb;
  if (!(stat(inputDir.str().c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))) {
    cout << "run" << run << " scope traces not found, defaulting to 10600" << endl;
    inputDir.str("");
    inputDir << calDir << "/antarctica14/scopeTraces/run10600/";
  }

    


  //set up for correlating and averaging
  TGraph* waveGraphs[numEntries];


  string dataLine;

  //get all the graphs out of the scope readout 
  for (int entry=0; entry<numEntries; entry++) {  
    inputFile.str("");
    inputFile << inputDir.str() << "wave_" << entry+1;
    
    waveGraphs[entry] = getScopeGraph(inputFile.str(),chan);
    
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
  runLogFile.open("runLog56dB.txt");
  if (!runLogFile) cout << "didn't open runlog file? check if runLog56dB.txt exists" << endl;

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



void surfInfoMaker(string antName, int numGraphs,TGraph **waveGraphs) {

  stringstream name;

  // Save a bunch of interesting things about generating the surf waveform (since it messes up so often)
  name.str("");
  name << outDir << antName << "_surfInfo.root";
  TFile *testFile = TFile::Open(name.str().c_str(),"recreate");

  //Save a copy of the histogram
  TH1D *hist = brotterTools::correlationDistribution(numGraphs,waveGraphs);
  TGraph *corrPattern = brotterTools::correlationPattern(numGraphs,waveGraphs);
  corrPattern->SetName("corrPattern1");
  corrPattern->SetTitle("corrPattern1");

  TGraph *corrPattern2 = brotterTools::correlationPattern(numGraphs,waveGraphs);
  corrPattern2->SetName("corrPattern2");
  corrPattern2->SetTitle("corrPattern2");

  TGraph *corrPattern3 = brotterTools::correlationPattern(numGraphs,waveGraphs);
  corrPattern3->SetName("corrPattern3");
  corrPattern3->SetTitle("corrPattern3");

  TGraph *corrPattern4 = brotterTools::correlationPattern(numGraphs,waveGraphs);
  corrPattern4->SetName("corrPattern4");
  corrPattern4->SetTitle("corrPattern4");

  //write the movie of the correlation, this takes forever (and is super inefficient) so you should probably leave it commented most of the time
  //also try doing this correlation and averaging myself?
  TGraph *myCorr = new TGraph(*waveGraphs[0]);

  for (int i=2; i<numGraphs; i++) {
    cout << "surfInfoMaker: " << i << "/" << numGraphs << "\r";
    TGraph *tempGraph = FFTtools::correlateAndAverage(i,waveGraphs);
    name.str("");
    name << "averagedGraph" << i;
    tempGraph->SetName(name.str().c_str());
    tempGraph->Write();

    TGraph *tempGraph2 = brotterTools::rotateToMatch(myCorr,waveGraphs[i]);
    name.str("");
    name << "rotated" << i;
    tempGraph2->SetName(name.str().c_str());
    tempGraph2->Write();
        
    for (int pt=0; pt<myCorr->GetN(); pt++) {
      myCorr->GetY()[pt] += tempGraph2->GetY()[pt];
    }
    name.str("");
    name << "myCorr" << i;
    myCorr->SetName(name.str().c_str());
    myCorr->Write();



    delete(tempGraph);
    delete(tempGraph2);

  }



  hist->Write();
  corrPattern->Write();
  corrPattern2->Write();
  corrPattern3->Write();
  corrPattern4->Write();
  //  averagedFFT->Write();

  testFile->Close();

  return;
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
  AnitaGeomTool *geom = AnitaGeomTool::Instance();
  int chanIndex = geom->getChanIndexFromRingPhiPol(ring,phi-1,pol);
  delete geom;

  cout << "phi:" << phi << " ring:" << AnitaRing::ringAsChar(ring) << " pol:" << AnitaPol::polAsChar(pol) << " chanIndex:" << chanIndex << endl;

  //set up where the data is
  stringstream name;
  stringstream dirName;
  char* calDir = getenv("ANITA3_CALDATA");
  dirName << calDir << "/antarctica14/root/";

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
    UsefulAnitaEvent *useful = new UsefulAnitaEvent(event,WaveCalType::kFull,head);

    //Get calibrated waveform for the channel of interest
    TGraph *calibGraph = useful->getGraph(ring,phi-1,pol);
    delete(useful);


    //FFTtools has a freq domain interpolator!
    //    TGraph *interpGraphA = FFTtools::getInterpolatedGraph(calibGraph,1./2.6);
    //    


    //Cosmin wants me to try stuff in RFInterpolate.h
    TGraphErrors *cosminGraph = FFTtools::getInterpolatedGraphSparseInvert(calibGraph,1/2.6,260);
    TGraph *interpGraphA = new TGraph(cosminGraph->GetN(),cosminGraph->GetX(),cosminGraph->GetY());
    TGraph *interpGraph = FFTtools::getInterpolatedGraphFreqDom(interpGraphA,0.1);

    //create  an Analysis Waveform object to do some other stuff
    //This also does an AKIMA spline to get even sampling (1/2.6GS/s = 
    //    AnalysisWaveform *aWave = new AnalysisWaveform(calibGraph->GetN(),calibGraph->GetX(),calibGraph->GetY());
    //    delete(calibGraph);
    //also time domain zero pad so they are the same length
    //    aWave->forceEvenSize(260);

    //need to ensure "beginning and ending of waveform go to zero
    



//    int nPts = aWave->Neven();
//    int nPtsFFT = nPts/2. - 1;
//    double dT = 1./2.6;
//    double new_dT = 1./10.;
//    int new_nPtsFFT = ((dT/new_dT)+1)*(nPtsFFT - 1);
//    int padding = new_nPtsFFT - nPtsFFT;
    //    cout << "nPts=" << nPts << " nPtsFFT=" << nPtsFFT << " new_nPtsFFT=" << new_nPtsFFT << " padding=" << padding << endl;
    

    //    aWave->padFreqAdd(padding);
    //    TGraph *interpGraph = new TGraph(aWave->Neven(),aWave->even()->GetX(),aWave->even()->GetY());
    //    delete(aWave);

    //upsample to 10GS/s (0.1ns bins)
    //    TGraph *interpGraph = FFTtools::getInterpolatedGraph(calibGraph,0.1);
    //    delete(calibGraph);
    //zero pad to 1024 bins
    TGraph *finalGraph = brotterTools::makeLength(interpGraph,1024);
    delete(interpGraph);
    //    printDT(finalGraph,"finalGraph");


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
  
  
  //if you want to save the surf info (about correlating and averaging mostly)
  //  surfInfoMaker(antName,entryAveraged,waveGraphs);


  

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

  AnitaVersion::set(3);

  string antName;

  string suffix = "";

  if (argc == 1) {
    cout << "Using default channel of 09BH" << endl;
    antName = "09BH";
  }
  else if (argc == 2) {
    cout << "Using " << argv[1] << " as channel" << endl;
    antName = argv[1];
  }
  else if (argc == 3) {
    cout << "Using " << argv[1] << " as channel" << endl;
    antName = argv[1];
    suffix = argv[2];
    cout << "Appending " << suffix << " to filename" << endl;
  }
  else {
    cout << "Usage: " << argv[0] << " [channelName (eg:09BH)] [opt: suffix]" << endl;
  }

  int calRun =findCalRun(antName);

  //Get the two pulses - and one more for the trigger!
  TGraph *waveGraphScope = scopeParseAndAverage(calRun,0);
  cout << "Got Scope Waveforms!" << endl;
  TGraph *waveGraphSurf = surfParseAndAverage(antName);
  cout << "Got Surf Waveforms!" << endl;
  TGraph *waveGraphTrig = scopeParseAndAverage(calRun,1);
  waveGraphTrig->SetName("waveGraphTrig");
  cout << "Got Trigger Waveforms!" << endl;

  //save them to a .txt file
  ofstream outStream;
  stringstream fileName;
  //Surf
  fileName.str("");
  fileName << outDir << antName << "_avgSurfWaveform";
  if (suffix != "") fileName << "_" << suffix;
  fileName << ".txt";
   outStream.open(fileName.str());
  for (int pt=0; pt<waveGraphSurf->GetN(); pt++) {
    outStream << waveGraphSurf->GetX()[pt] << " " << waveGraphSurf->GetY()[pt] << endl;
  }
  outStream.close();

  //Scope 
  fileName.str("");
  fileName << outDir << antName << "_avgScopeWaveform";
  if (suffix != "") fileName << "_" << suffix;
  fileName << ".txt";  
  outStream.open(fileName.str());
  for (int pt=0; pt<waveGraphScope->GetN(); pt++) {
    outStream << waveGraphScope->GetX()[pt] << " " << waveGraphScope->GetY()[pt] << " " << waveGraphTrig->GetY()[pt] << endl;
  }
  outStream.close();

  //save them to a .root file
  fileName.str("");
  fileName << outDir << antName << "_avgWaveforms";
  if (suffix != "") fileName << "_" << suffix;
  fileName << ".root";  
  TFile *rootFile = TFile::Open(fileName.str().c_str(),"recreate");
  waveGraphSurf->Write();
  waveGraphScope->Write();
  waveGraphTrig->Write();
  rootFile->Close();


  delete waveGraphScope;
  delete waveGraphSurf;

  cout << "Done! :)" << endl;

  

  return 1;
}




