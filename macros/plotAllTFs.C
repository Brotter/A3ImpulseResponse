{
//TString thename = "/home/abl/impulseResponseCode/A3ImpulseResponse/transferFunctions/palAntA4_";
  TString thename = "/home/abl/impulseResponseCode/A3ImpulseResponse/oldTransferFunctions/palAnt_";
  TString end = ".txt";
  TGraph* gH[48];
  TGraph* gV[48];
  int cH = 0;
  int cV = 0;
  for(int phi = 1; phi < 17; phi++)
  {
    for(int tmb = 0; tmb < 3; tmb++)
    {
      for(int pol = 0; pol < 2; pol++)
      {
        TString ring;
        if(tmb == 0) ring = "T";
        if(tmb == 1) ring = "M";
        if(tmb == 2) ring = "B";
        TString p = pol ? "V" : "H";
        TString num = "0";
        if(phi < 10)
        {
          TString tempNum = Form("%d", phi);
          num = num + tempNum;
        }
        if(phi > 9) num = Form("%d", phi);
        TString fname = thename + num + ring + p + end;
        if(pol == 0)
        {
          gH[cH] = new TGraph(fname.Data());
          cH++;
        }
        if(pol == 1)
        {
          gV[cV] = new TGraph(fname.Data());
          cV++;
        }
      }
    }
  }
  TCanvas * c1 = new TCanvas("c1","c1", 1000,400);
  c1->cd();
  gH[0]->GetYaxis()->SetRangeUser(-.09,.09);
  gH[0]->GetXaxis()->SetRangeUser(0,20);
  gH[0]->Draw("alp");
  for(int i = 1; i < 48; i++) gH[i]->Draw("lpsame");
  TCanvas * c2 = new TCanvas("c2","c2", 1000,400);
  c2->cd();
  gV[0]->Draw("alp");
  gV[0]->GetYaxis()->SetRangeUser(-.09,.09);
  gV[0]->GetXaxis()->SetRangeUser(0,20);
  for(int i = 1; i < 48; i++) gV[i]->Draw("lpsame");
}
