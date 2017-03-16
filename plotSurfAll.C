void plotSurfAll(string antName){

  stringstream name;
  name.str("");
  name << "waveforms_new/" << antName << "_surfInfo.root";
  TFile *inFile = TFile::Open(name.str().c_str());


  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);

  for (int i=2; i<500; i++) {
    name.str("");
    name << "myCorr" << i;
    TGraph *currGraph = (TGraph*)inFile->Get(name.str().c_str());
    for (int pt=0; pt<currGraph->GetN(); pt++) currGraph->GetY()[pt] /= i;
    name.str("");
    name << "rotated" << i;
    TGraph *rotGraph = (TGraph*)inFile->Get(name.str().c_str());
    rotGraph->SetLineColor(kRed);
      
    name.str("");
    name << "ev" << i;
    currGraph->SetTitle(name.str().c_str());
    currGraph->Draw("aC*");
    currGraph->GetHistogram()->SetMaximum(120);
    currGraph->GetHistogram()->SetMinimum(-120);
    currGraph->Draw("aC*");
    rotGraph->Draw("lSame");
    c1->Update();
    //    if (i==2) c1->SaveAs("08MH.gif");
    //    else c1->SaveAs("08MH.gif+");
  }


}
