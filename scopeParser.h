#ifndef _SCOPEPARSER
#define _SCOPEPARSER

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

///PARSING SCOPE DATA
void ParseHeader(string line,double& dt, double& t0, double& yMulti, double& yOff, double& yZero);
void ParseData(string line,vector<double>& v,double yOff, double yMulti);
TGraph* VectorToGraph(vector<double> v,double dt);


TGraph *getScopeGraph(string fileNameBase,int channel) {
  ifstream inData,inHead;
  stringstream name;
  name.str("");
  name << fileNameBase << ".dat";
  inData.open(name.str());
  
  name.str("");
  name << fileNameBase << ".hdr";
  inHead.open(name.str());
    
  if (!inData.is_open() || !inHead.is_open()) {
    cout << fileNameBase << "Doesn't exist!!  Scope data or header files are missing!!!!" << endl;
  }
  
  double dt,t0,yMulti,yOff,yZero;

  string headerLine, dataLine;
  
  for (int i=0; i<=channel*2; i++) {
    getline(inHead,headerLine);
    getline(inData,dataLine);
  }
  ParseHeader(headerLine,dt,t0,yMulti,yOff,yZero);
  vector<double> dataVector;
  ParseData(dataLine,dataVector,yOff,yMulti);

  inData.close();
  inHead.close();
  
  return VectorToGraph(dataVector,dt);
}


void ParseHeader(string line,double& dt, double& t0, double& yMulti, double& yOff, double& yZero) {

  int pos;
  pos =line.find("XINCR");
  if (pos != string::npos) {
    int locpos = line.find(";XZE");
    stringstream ss(line.substr(pos +6,locpos - (pos + 6)));
    ss >> dt;
    //    cout << "dt: " << dt << endl;
  }
  
  pos =line.find("XZERO ");
  if (pos != string::npos) {
    int locpos = line.find(";PT_OF");
    stringstream ss(line.substr(pos +6,locpos - (pos + 6)));
    ss >> t0;
    //    cout << "t0: " << t0 << endl;
  }
  
  pos =line.find("YMULT");
  if (pos != string::npos) {
    int locpos = line.find(";YOFF");
    stringstream ss(line.substr(pos +6,locpos - (pos + 6)));
    ss >> yMulti;
    //    cout << "Ymultiplier: " << yMulti << endl;
  }
  
  pos =line.find("YOFF");
  if (pos != string::npos) {
    int locpos = line.find(";YZE");
    stringstream ss(line.substr(pos +6,locpos - (pos + 6)));
    ss >> yOff;
    //    cout << "yOff: " << yOff << endl;
  }
  
  pos =line.find("YZERO");
  if (pos != string::npos) {
    int locpos = line.find(";NR_FR");
    stringstream ss(line.substr(pos +6,locpos - (pos + 6)));
    ss >> yZero;
    //    cout << "yZero: " << yZero << endl;
  }
}

void ParseData(string line,vector<double>& v,double yOff, double yMulti) {


  for (int i = 0; i < line.size(); i++) {
    if (line[i] == ',') line[i] = ' '; }
  
  stringstream ss(line);
  double V;
  while (ss >> V) {
    v.push_back((V  + yOff) * yMulti);
  }
}


TGraph* VectorToGraph(vector<double> v,double dt) {
  TGraph *gr = new TGraph();
  for (int i = 0; i < v.size(); i++) {
    gr->SetPoint(i,i*dt,v[i]); }
  return gr;
}
  

#endif
