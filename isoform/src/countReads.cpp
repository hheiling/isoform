/*****************************************************************************
  intersectMain.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licenced under the GNU General Public License 2.0 license.
 
  Modified by Wei Sun, Aug 7th, 2011
  Department of Biostatistics, 
  Univerisity of North Carolina, Chapel Hill
  weisun@email.unc.edu

  Licenced under the GNU General Public License 2.0 license.
 
******************************************************************************/
#include "intersectBed.h"

extern "C" {

using namespace std;

int countReads(char** RbedAFile, char** RbedBFile, char** RoutputFile, 
               double* RoverlapFraction, int* RforceStrand) {

  // input files
  string bedAFile;
  string bedBFile;

  bedAFile = RbedAFile[0];
  bedBFile = RbedBFile[0];

  // output files
  char* outputFile;
  outputFile = RoutputFile[0];
  
  // input arguments
  float overlapFraction = (float) RoverlapFraction[0];
  bool forceStrand      = false;
  if (RforceStrand[0]) { forceStrand = true; }

  BedIntersect *bi = new BedIntersect(bedAFile, bedBFile, outputFile, 
                                      overlapFraction, forceStrand);
  delete bi;
  return 0;
}

}// extern "C"
