#ifndef _ANTSCOPEPARSER
#define _ANTSCOPEPARSER

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

///PARSING SCOPE DATA
void ParseData(string line,vector<double> &t, vector<double> &v);
TGraph* VectorToGraph(vector<double> t, vector<double> v);


TGraph *getScopeGraph(ifstream &inData) {
  int currWave = 0;
  int lineNum = 0;

  string dataLine;
  
  vector<double> dataVector;
  vector<double> timeVector;


  getline(inData,dataLine);

  //read lines until you see a '#'
  while (dataLine[0] != '#') {
    lineNum++;
    ParseData(dataLine,timeVector, dataVector);
    getline(inData,dataLine);

   }
  
  //  cout << "waveform was " << lineNum << " lines long" << endl;

  return VectorToGraph(timeVector,dataVector);
}

void ParseData(string line,vector<double> &t, vector<double>& v) {
  string sT,sV;
  double T,V;

  sT = line.substr(0,line.find(" "));
  T = stod(sT)*1e9; //make it in ns!!
  sV = line.substr(line.find(" "),line.npos);
  V = stod(sV);

  //  cout << T << " " << V << endl;

  v.push_back(V);
  t.push_back(T);

}


TGraph* VectorToGraph(vector<double> t,vector<double> v) {
  TGraph *gr = new TGraph();
  for (int i = 0; i < v.size(); i++) {
    gr->SetPoint(i,t[i],v[i]); }
  return gr;
}
  

#endif
