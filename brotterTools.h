
#ifndef BROTTERTOOLS_H
#define BROTTERTOOLS_H

#ifdef __CINT__
#ifdef FFTW_64_BIT // Hack for Hawaii install of FFTW
typedef struct {char a[16];} __float128; /* 16 chars have the same size as one __float128 */
#endif
#endif 

// c++ libraries thingies
#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>

// ROOT
#include "TGraph.h"
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"
#include "TH1D.h"

//My includes
#include "FFTWComplex.h"
#include "RFSignal.h"
#include "FFTtools.h"


// FFTW
#include <fftw3.h>

/* I shamelessly stole this first part from Ryan, and the Makefile from Ben Strutt, and I'm going to dump a bunch of the utilities I
   use all the time in every piece of code I write and try to make a library I can move around

*/


class brotterTools
{
 public:
  brotterTools(); ///< Constructor
  ~brotterTools(); ///< Destructor


  /*A shameless copy of Ryan's FFTtools::correlateAndAverage(), except it just gives you back a histogram of the correlation values*/
  static TH1D *correlationDistribution(Int_t numGraphs, TGraph **grPtrPtr);
  static TGraph *correlationPattern(Int_t numGraphs, TGraph **grPtrPtr);
  static TGraph *zeroPadToLength(const TGraph *inGraph, Int_t endLength);

};

#endif //BROTTERTOOLS_H

