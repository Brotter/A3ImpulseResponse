#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include <time.h>
#include "FFTtools.h"
#include "AnitaDataset.h"
#include "UsefulAnitaEvent.h"
#include "AnalysisWaveform.h"
#include "TH2F.h"
#include <unistd.h>
#include "AnitaGeomTool.h"
#include "TGraphErrors.h"
#include "RFInterpolate.h"

TGraph* cutToLength(const TGraph *inGraph, Int_t endLength)
{
  /* 
   *      Cuts off the final points of a TGraph so that it is endLength long
   *           Creates and returns a new TGraph
   *             */

  //make output graph
  TGraph *outGraph = new TGraph(*inGraph);

  //    //check if you are actually trying to shorten something, if you are, just copy the initial graph and return it
  if (inGraph->GetN() <= endLength) {
    cout << "Warning from zeroPadToLength: Graph is already shorter than you want!! ";
    cout << "(" << inGraph->GetN() << " >= " << endLength << ")" << endl;
    delete(outGraph);
    outGraph = new TGraph(*inGraph);
    return outGraph;
  }

  for (int pt=0; pt<endLength; pt++) {
    outGraph->SetPoint(pt,inGraph->GetX()[pt],inGraph->GetY()[pt]);
  }

  //done!
  return outGraph;

}

TGraph* zeroPadToLength(const TGraph *inGraph, Int_t endLength)
{ 
  /* 
   *      Zero pad a TGraph to have endLength points
   *           Creates and returns a new TGraph
   *             */

  //make output graph
  TGraph *outGraph = new TGraph(*inGraph);
  //      
  //check if you are actually trying to shorten something, if you are, just copy the initial graph and return it
  if (inGraph->GetN() >= endLength) {
    cout << "Warning from zeroPadToLength: Graph is already longer than you want!! ";
    cout << "(" << inGraph->GetN() << " >= " << endLength << ")" << endl;
    delete(outGraph);
    outGraph = new TGraph(*inGraph);
    return outGraph;
  }

  //find out what the binning will be, as well as inGraph start point
  double tStart = inGraph->GetX()[0];
  double dT = inGraph->GetX()[1] - tStart;

  //zero pad end (remember odd value!)
  while (outGraph->GetN() < endLength) {
    int pt = outGraph->GetN();
    outGraph->SetPoint(pt,outGraph->GetX()[pt-1]+dT,0);
  }

  //done!
  return outGraph;
}

TGraph* makeLength(const TGraph *inGraph, Int_t endLength)
{
  /*
   *
   *     Makes a graph the exact length you want by cutting or zero padding the end
   *
   *         Returns a new TGraph
   *
   *            */

  TGraph *outGraph = NULL;
  if (inGraph->GetN() < endLength) {
    outGraph = zeroPadToLength(inGraph,endLength);
  }
  else if (inGraph->GetN() > endLength) {
    outGraph = cutToLength(inGraph,endLength);
  }
  else {
    outGraph = new TGraph(*inGraph);
  }

  return outGraph;

}

double bestTime(TGraph * corrGraph)
{
	double best;
	int index = TMath::LocMax(corrGraph->GetN(), corrGraph->GetY());
	best = corrGraph->GetX()[index];
	return best;
}

double VPP(TGraph * g)
{
	return (TMath::MaxElement(g->GetN(), g->GetY()) - TMath::MinElement(g->GetN(), g->GetY()));
}

void histos(int run, TString names[], AnitaRing::AnitaRing_t whichRing, int chans)
{
	long t0 = time(0);

	AnitaDataset * d = new AnitaDataset(run);
	d->setCalType(WaveCalType::kNoChannelToChannelDelays);
	const int numEntries = 150;//d->N();
	int startEvent = 0;
	UsefulAnitaEvent * uae = 0;

	TGraph* gWaves0[numEntries] = {NULL};
	TGraph* gWaves1[numEntries] = {NULL};
	TGraph* gWaves2[numEntries] = {NULL};
	TGraph* gWaves3[numEntries] = {NULL};
	TGraph* gWaves4[numEntries] = {NULL};
	TGraph* gWaves5[numEntries] = {NULL};
	TGraph* gWaves6[numEntries] = {NULL};
	TGraph* gWaves7[numEntries] = {NULL};
	TGraph* gWaves8[numEntries] = {NULL};
	TGraph* gWaves9[numEntries] = {NULL};
	TGraph* gWaves10[numEntries] = {NULL};
	TGraph* gWaves11[numEntries] = {NULL};
	TGraph* gWaves12[numEntries] = {NULL};
	TGraph* gWaves13[numEntries] = {NULL};
	TGraph* gWaves14[numEntries] = {NULL};
	TGraph* gWaves15[numEntries] = {NULL};

	AnitaRing::AnitaRing_t midRing = AnitaRing::kMiddleRing;
	AnitaRing::AnitaRing_t botRing = AnitaRing::kBottomRing;
	AnitaRing::AnitaRing_t topRing = AnitaRing::kTopRing;

	AnitaPol::AnitaPol_t HPOL = AnitaPol::kHorizontal;
	AnitaPol::AnitaPol_t VPOL = AnitaPol::kVertical;

	int nChan = 16;

	//This is the ring each of the 16 cables
	AnitaRing::AnitaRing_t rings[16];
	//this is the pol of each of the 16 cables
	AnitaPol::AnitaPol_t pols[16];
	//this is the channel of each of the 16 cables
	int chanNum[16];
  int finished[16] = {0};
	
	double refDelay = 81.3919;

	double tDelay = 65.9021;

	double delays[16] = {25.625, 22.827, 23.028, 25.624,
											 32.104, 22.988, 22.848, 25.64,
											 22.986, 22.824, 25.751, 22.932,
											 25.514, 25.722, 25.514, 25.547};

	TGraph* gAverage[16] = {NULL};
	for (int i = 0; i < 16; i++)
	{
		rings[i] = whichRing;
		chanNum[i] = (i)%8 + chans;
		if (i < 8) pols[i] = VPOL;
		if (i >= 8) pols[i] = HPOL;
		TString vppname = names[i] + "VPP";
		gAverage[i] = 0;
	}


	TGraph * gChan = 0;
	double bt = 0;
	
	AnitaGeomTool * agt;
	
	int surf;
	int chan;
	int ant;

  int entryUsed[16] = {0};

  for (int i = startEvent; i < d->N(); i++)
  {
    int nchandone = 0;
    for(int phi = 0; phi < nChan; phi++) nchandone += finished[phi];
    if(nchandone == 16) break;
    d->getEntry(i);
    uae = d->useful();
    if(i%10 == 0) printf("%d events finished!\n", i);
    for (int phi = 0; phi < nChan; phi++)
    {
      if(entryUsed[phi] == numEntries)
      {
        finished[phi] = 1;
        continue;
      }
      int chanIndex = agt->getChanIndexFromRingPhiPol(rings[phi], chanNum[phi], pols[phi]);
      if(uae->getLabChip(chanIndex) != 0) continue;
      gChan=uae->getGraph(rings[phi], chanNum[phi], pols[phi]);
      AnalysisWaveform wf(gChan->GetN(), gChan->GetX(), gChan->GetY(), .1);
      const TGraphAligned* temp = wf.even();
      TGraph* interpGraph = new TGraph(temp->GetN(), temp->GetX(), temp->GetY());
      TGraph* gAligned = makeLength(interpGraph, 1024);

      if(phi==0) gWaves0[entryUsed[0]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==1) gWaves1[entryUsed[1]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==2) gWaves2[entryUsed[2]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==3) gWaves3[entryUsed[3]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==4) gWaves4[entryUsed[4]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==5) gWaves5[entryUsed[5]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==6) gWaves6[entryUsed[6]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==7) gWaves7[entryUsed[7]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==8) gWaves8[entryUsed[8]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==9) gWaves9[entryUsed[9]]=new TGraph(gAligned->GetN(), gAligned->GetX(), gAligned->GetY());
      if(phi==10) gWaves10[entryUsed[10]]=new TGraph(gAligned->GetN(),gAligned->GetX(),gAligned->GetY());
      if(phi==11) gWaves11[entryUsed[11]]=new TGraph(gAligned->GetN(),gAligned->GetX(),gAligned->GetY());
      if(phi==12) gWaves12[entryUsed[12]]=new TGraph(gAligned->GetN(),gAligned->GetX(),gAligned->GetY());
      if(phi==13) gWaves13[entryUsed[13]]=new TGraph(gAligned->GetN(),gAligned->GetX(),gAligned->GetY());
      if(phi==14) gWaves14[entryUsed[14]]=new TGraph(gAligned->GetN(),gAligned->GetX(),gAligned->GetY());
      if(phi==15) gWaves15[entryUsed[15]]=new TGraph(gAligned->GetN(),gAligned->GetX(),gAligned->GetY());
      entryUsed[phi]++;
      delete interpGraph;
      delete gAligned;
    }
  }
  gAverage[0] = FFTtools::correlateAndAverage(entryUsed[0], gWaves0);
  gAverage[1] = FFTtools::correlateAndAverage(entryUsed[1], gWaves1);
  gAverage[2] = FFTtools::correlateAndAverage(entryUsed[2], gWaves2);
  gAverage[3] = FFTtools::correlateAndAverage(entryUsed[3], gWaves3);
  gAverage[4] = FFTtools::correlateAndAverage(entryUsed[4], gWaves4);
  gAverage[5] = FFTtools::correlateAndAverage(entryUsed[5], gWaves5);
  gAverage[6] = FFTtools::correlateAndAverage(entryUsed[6], gWaves6);
  gAverage[7] = FFTtools::correlateAndAverage(entryUsed[7], gWaves7);
  gAverage[8] = FFTtools::correlateAndAverage(entryUsed[8], gWaves8);
  gAverage[9] = FFTtools::correlateAndAverage(entryUsed[9], gWaves9);
  gAverage[10] = FFTtools::correlateAndAverage(entryUsed[10], gWaves10);
  gAverage[11] = FFTtools::correlateAndAverage(entryUsed[11], gWaves11);
  gAverage[12] = FFTtools::correlateAndAverage(entryUsed[12], gWaves12);
  gAverage[13] = FFTtools::correlateAndAverage(entryUsed[13], gWaves13);
  gAverage[14] = FFTtools::correlateAndAverage(entryUsed[14], gWaves14);
  gAverage[15] = FFTtools::correlateAndAverage(entryUsed[15], gWaves15);

  for(int i = 0; i < entryUsed[0]; i++) delete gWaves0[i];
  for(int i = 0; i < entryUsed[1]; i++) delete gWaves1[i];
  for(int i = 0; i < entryUsed[2]; i++) delete gWaves2[i];
  for(int i = 0; i < entryUsed[3]; i++) delete gWaves3[i];
  for(int i = 0; i < entryUsed[4]; i++) delete gWaves4[i];
  for(int i = 0; i < entryUsed[5]; i++) delete gWaves5[i];
  for(int i = 0; i < entryUsed[6]; i++) delete gWaves6[i];
  for(int i = 0; i < entryUsed[7]; i++) delete gWaves7[i];
  for(int i = 0; i < entryUsed[8]; i++) delete gWaves8[i];
  for(int i = 0; i < entryUsed[9]; i++) delete gWaves9[i];
  for(int i = 0; i < entryUsed[10]; i++) delete gWaves10[i];
  for(int i = 0; i < entryUsed[11]; i++) delete gWaves11[i];
  for(int i = 0; i < entryUsed[12]; i++) delete gWaves12[i];
  for(int i = 0; i < entryUsed[13]; i++) delete gWaves13[i];
  for(int i = 0; i < entryUsed[14]; i++) delete gWaves14[i];
  for(int i = 0; i < entryUsed[15]; i++) delete gWaves15[i];

	TString fgname = "averageGraphs/run" + TString::Itoa(run, 10) + "wfs.root";
	TFile fg(fgname.Data(), "RECREATE");
	fg.cd();
	for (int i = 0; i < nChan; i++)
	{
		TString avgtitle = names[i] + "wf";
		gAverage[i]->SetTitle(avgtitle.Data());
		gAverage[i]->Write(avgtitle.Data());
	}
	fg.Close();	
	for (int i = 0; i < nChan; i++)
	{
    delete gAverage[i];
  }
	
	printf("%ld\n", time(0) - t0);
}

void doAll()
{
/*	
*/	
	TString names1858[16] = { "MV1", "MV2", "MV3", "MV4",
                            "MV5", "MV6", "MV7", "MV8",
                            "MH1", "MH2", "MH3", "MH4",
                            "MH5", "MH6", "MH7", "MH8"};

	histos(1858, names1858, AnitaRing::kMiddleRing, 0); 
	
/*	
	TString names1879[16] = { "MV9", "MV10", "MV11", "MV12",
                            "MV13", "MV14", "MV15", "MV16",
                            "MH9", "MH10", "MH11", "MH12",
                            "MH13", "MH14", "MH15", "MH16"};

	histos(1879, names1879, AnitaRing::kMiddleRing, 8); 

	TString names1884[16] = { "BV9", "BV10", "BV11", "BV12",
                            "BV13", "BV14", "BV15", "BV16",
                            "BH9", "BH10", "BH11", "BH12",
                            "BH13", "BH14", "BH15", "BH16"};

	histos(1884, names1884, AnitaRing::kBottomRing, 8); 

	TString names1891[16] = { "BV1", "BV2", "BV3", "BV4",
                            "BV5", "BV6", "BV7", "BV8",
                            "BH1", "BH2", "BH3", "BH4",
                            "BH5", "BH6", "BH7", "BH8"};

	histos(1891, names1891, AnitaRing::kBottomRing, 0);

	
	TString names1909[16] = { "TV9", "TV10", "TV11", "TV12",
                            "TV13", "TV14", "TV15", "TV16",
                            "TH9", "TH10", "TH11", "TH12",
                            "TH13", "TH14", "TH15", "TH16"};

	histos(1909, names1909, AnitaRing::kTopRing, 8); 

	TString names1900[16] = { "TV1", "TV2", "TV3", "TV4",
                            "TV5", "TV6", "TV7", "TV8",
                            "TH1", "TH2", "TH3", "TH4",
                            "TH5", "TH6", "TH7", "TH8"};

	histos(1900, names1900, AnitaRing::kTopRing, 0);
*/	
}
