#include "brotterTools.h"

using namespace std;

brotterTools::brotterTools()
{
}

brotterTools::~brotterTools()
{
}

TGraph *brotterTools::correlationPattern(Int_t numGraphs, TGraph **grPtrPtr)
{



  TGraph *corrPattern = new TGraph();
 
  //Assume they are all at same sampling rate
  if(numGraphs<2) return NULL;
  TGraph *grA = new TGraph(*grPtrPtr[0]);
  Int_t numPoints=grA->GetN();  

  double dT = grA->GetX()[1]-grA->GetX()[0];

  Double_t *timeVals= grA->GetX();
  Double_t *safeTimeVals = new Double_t[numPoints];
  Double_t *sumVolts = new Double_t [numPoints];
  for(int i=0;i<numPoints;i++) 
    safeTimeVals[i]=timeVals[i];  
  
  int countWaves=1;
  for(int graphNum=1;graphNum<numGraphs;graphNum++) {
    TGraph *grB = grPtrPtr[graphNum];
    if(grB->GetN()<numPoints)
      numPoints=grB->GetN();
    TGraph *grCorAB = FFTtools::getCorrelationGraph(grA,grB);

    Int_t peakBin = FFTtools::getPeakBin(grCorAB);
    corrPattern->SetPoint(corrPattern->GetN(),graphNum,peakBin*dT);
    TGraph *grComAB = new TGraph(numPoints,safeTimeVals,sumVolts);

    //    delete grB;
    delete grCorAB;
    if(graphNum>1)
      delete grA;
    grA=grComAB;
    countWaves++;

  }
  delete grA;
  delete [] safeTimeVals;
  delete [] sumVolts;


  return corrPattern;


}

TH1D *brotterTools::correlationDistribution(Int_t numGraphs, TGraph **grPtrPtr)
{
  TGraph *grA = grPtrPtr[0];

  Int_t numPoints=grA->GetN();  
  int corrLength=int(TMath::Power(2,int(TMath::Log2(numPoints))+2));
  double dT = grA->GetX()[1]-grA->GetX()[0];
  cout << "dT= " << dT << " corrLength= " << corrLength << endl;

  TGraph *corrPattern = correlationPattern(numGraphs,grPtrPtr);

  TH1D *outputHist = new TH1D("correlationDistribution","correlationDistribution",corrLength,0,corrLength*dT);
  for (int pt=0; pt<corrPattern->GetN(); pt++) {
    outputHist->Fill(corrPattern->GetY()[pt]);
  }

  return outputHist;
  
 }


TGraph *brotterTools::zeroPadToLength(const TGraph *inGraph, Int_t endLength) 
{
  /* 
     Zero pad a TGraph to have endLength points
     Creates and returns a new TGraph
  */

  //make output graph
  TGraph *outGraph = new TGraph();

  //check if you are actually trying to shorten something, if you are, just copy the initial graph and return it
  if (inGraph->GetN() >= endLength) {
    cout << "Warning from zeroPadToLength: Graph is already longer than you want!! ";
    cout << "(" << inGraph->GetN() << " >= " << endLength << ")" << endl;
    delete(outGraph);
    outGraph = new TGraph(*inGraph);
    return outGraph;
  }
  
  //See how much longer you need to make it
  int needToAdd = endLength - inGraph->GetN();
  //if it is odd, note that (add the extra one on to the end)
  int odd = 0;
  if (needToAdd%2 != 0) odd = 1;

  //print out values as a check
  //  cout << "beginning length: " << inGraph->GetN() << " end length:" << endLength << " adding: " << needToAdd << endl;
  

  //find out what the binning will be, as well as inGraph start point
  double tStart = inGraph->GetX()[0];
  double dT = inGraph->GetX()[1] - tStart;

  //zero pad front
  for (int pt=0; pt<needToAdd/2; pt++) {
    outGraph->SetPoint(outGraph->GetN(),tStart-(needToAdd/2-pt)*dT,0);
  }
  //fill in inGraph
  for (int pt=0; pt<inGraph->GetN(); pt++) {
    outGraph->SetPoint(outGraph->GetN(),inGraph->GetX()[pt],inGraph->GetY()[pt]);
  }
  //zero pad end (remember odd value!)
  for (int pt=0; pt<needToAdd/2 + odd; pt++) {
    outGraph->SetPoint(outGraph->GetN(),outGraph->GetX()[outGraph->GetN()-1]+dT,0);
  }

  //done!
  return outGraph;

}
