/*****************************************************************************
  intersectBed.cpp

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
#include "lineFileUtilities.h"
#include "intersectBed.h"


/*
    Constructor
*/
BedIntersect::BedIntersect(string bedAFile, string bedBFile, const char* outputFile, 
                           float overlapFraction, bool forceStrand) {

    _bedAFile        = bedAFile;
    _bedBFile        = bedBFile;
    _outputFile      = outputFile;
    _overlapFraction = overlapFraction;
    _forceStrand     = forceStrand;

    // create new BED file objects for A and B
    _bedA = new BedFile(bedAFile);
    _bedB = new BedFile(bedBFile);

    IntersectBam(bedAFile, outputFile);
}


/*
 Destructor
 */
BedIntersect::~BedIntersect(void) {
  if(_bedA != NULL){
    delete _bedA;
    _bedA= NULL;
  }
  if(_bedB != NULL){
    delete _bedB;
    _bedB= NULL;
  }
  
}


bool BedIntersect::FindOverlaps(const BED &a, vector<BED> &hits, set<string> &nms) {
  
  // how many overlaps are there b/w the bed and the set of hits?
  int s, e, overlapBases=0;
  bool hitsFound   = false;
  int aLength      = (a.end - a.start);   // the length of a in b.p.
  
  // collect and report the sufficient hits
  _bedB->FindOverlapsPerBin(a.chrom, a.start, a.end, a.strand, hits, _forceStrand);  
  
  // loop through the hits and report those that meet the user's criteria
  vector<BED>::const_iterator h       = hits.begin();
  vector<BED>::const_iterator hitsEnd = hits.end();
    
  for (; h != hitsEnd; ++h) {
    s            = max(a.start, h->start);
    e            = min(a.end, h->end);
    overlapBases = overlapBases + (e - s); // the number of overlapping bases b/w a and b    
    nms.insert(h->name);
  }
  
  // is there enough overlap relative to the user's request? (default ~ 1bp)
  // 'a' is a consecutive piece of DNA sequence. However it can overlap with 
  // more than one consecutive exons
  if ( ( (float) overlapBases / (float) aLength ) >= _overlapFraction ) {
    hitsFound = true;
  }
  
  return hitsFound;
}



void BedIntersect::IntersectBam(string bamFile, const char* outputFile) {

  bool foundIt = true;

  // load the "B" bed file into a map so
  // that we can easily compare "A" to it for overlaps
  _bedB->loadBedFileIntoMap();

  // open the BAM file
  BamReader reader;
  BamWriter writer;
  reader.Open(bamFile);

  // get header & reference information
  string header  = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  vector<BED> hits;
  // reserve some space
  hits.reserve(100);

  _bedA->bedType = 6;
  BamAlignment bam;

  // name of the previous read
  string name0 ("");
  string name1 ("");

  // names of the elements in hits, e.g., the features in fileB
  set<string> nms;
  set<string>::iterator it;

  // output file
  ofstream outfile;
  outfile.open (outputFile);

  // get each set of alignments for each pair.
  while (reader.GetNextAlignment(bam)) {

    if (! bam.IsMapped()) {
      continue;
    }
    
    name1 = bam.Name;
    name1 = name1.substr(0, name1.find_last_not_of('/')-1);

    if (name0.compare(name1) != 0) {
      
      if (! nms.empty()) {
        
        if(foundIt){
          outfile << name0 << "\t";
          for ( it=nms.begin() ; it != nms.end(); it++ ){
            outfile << *it << ";";
          }
          
          outfile << "\n";
        }
        
        nms.clear();
      }
      
      name0 = name1;
      foundIt = true;
    }
    
    // split the BAM alignment into discrete BED blocks and
    // look for overlaps only within each block.
    bedVector bedBlocks;  // vec to store the discrete BED "blocks" from a
    getBamBlocks(bam, refs, bedBlocks, false);
    
    vector<BED>::const_iterator bedItr  = bedBlocks.begin();
    vector<BED>::const_iterator bedEnd  = bedBlocks.end();
    
    for (; bedItr != bedEnd; ++bedItr) {
      foundIt = FindOverlaps(*bedItr, hits, nms);
      hits.clear();
      if (! foundIt) { break; }
    }
    
  }

  if(foundIt){
    outfile << name0 << "\t";
    for ( it=nms.begin() ; it != nms.end(); it++ ){
      outfile << *it << ";";
    }
    
    outfile << "\n";
  }
  
  nms.clear();
  
  // close the relevant BAM files.
  reader.Close();

  // close output file
  outfile.close();
}

