void plotSurfAll(string antName){

  stringstream name;
  name.str("");
  name << "waveforms/" << antName << "_surfInfo.root";
  TFile *inFile = TFile::Open(name.str().c_str());


  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);

  for (int i=2; i<500; i++) {
    name.str("");
    name << "averagedGraph" << i;
    TGraph *currGraph = (TGraph*)inFile->Get(name.str().c_str());
    currGraph->Draw("alp");
    c1->Update();
  }


}
