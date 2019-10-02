/*****************************************************************************
  intersectBed.h

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
#ifndef INTERSECTBED_H
#define INTERSECTBED_H

#include "bedFile.h"
#include "BamReader.h"
#include "BamWriter.h"
#include "BamAncillary.h"
#include "BamAux.h"
using namespace BamTools;


#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;



class BedIntersect {

public:

    // constructor
    BedIntersect(string bedAFile, string bedBFile, const char* outputFile,
                 float overlapFraction, bool forceStrand);

    // destructor
    ~BedIntersect(void);

private:

    //------------------------------------------------
    // private attributes
    //------------------------------------------------
    string _bedAFile;
    string _bedBFile;
    string _outputFile;
  
    bool  _forceStrand;
    float _overlapFraction;

    // instance of a bed file class.
    BedFile *_bedA, *_bedB;

    //------------------------------------------------
    // private methods
    //------------------------------------------------

    void IntersectBam(string bamFile, const char* outputFile);

    bool FindOverlaps(const BED &a, vector<BED> &hits, set<string> &nms);

};

#endif /* INTERSECTBED_H */
