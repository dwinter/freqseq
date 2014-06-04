#ifndef utils_H
#define utils_H

#include "utils/bamtools_pileup_engine.h"

using namespace std;
using namespace BamTools;

//typededs and classes

union ReadData{
    uint16_t reads[4];
    uint64_t key;
};

typedef vector<ReadData> ReadDataVector;

typedef pair<uint16_t, uint64_t> SiteSummary;
typedef vector<SiteSummary> SiteSummaryVector;
typedef vector<string> SampleNames;


//helper functions for BAM parsing

uint16_t get_base_index(char b){
    switch(b){        
        case 'A':
        case 'a':    
            return 0;
         case 'T':
         case 't':
             return 3;
         case 'C':
         case 'c':
             return 1;
         case 'G':
         case 'g':
             return 2;
         case '-':
         case 'N':
             return -1 ;
         default: // Unknown base, alert in debug mode?
             return -1;
    }
}


uint32_t find_sample_index(string s, SampleNames sv){
    for (size_t i=0; i < sv.size(); i++){
        if(s.compare(sv[i])==0){
            return i;
        }
    }
    return(-1);
}

bool passes_QC(PileupAlignment pileup, uint16_t map_cut, uint16_t qual_cut){
    const BamAlignment *ali = &pileup.Alignment;
    if(ali->MapQuality > map_cut){
        uint16_t bqual = static_cast<short>(ali->Qualities[pileup.PositionInAlignment]) - 33;
        if(bqual > qual_cut){
            return(not (ali->IsDuplicate()) && not(ali->IsFailedQC()) && ali->IsPrimaryAlignment());
        }
    }
    return false;
}

SiteSummary summarize_site(ReadData sample){
    uint16_t most_common = distance(sample.reads, max_element(sample.reads, sample.reads+4));
    uint64_t minor_allele_count = 0;
    uint64_t depth = 0;
    for(size_t i = 0; i < 4; i++){
        if(i != most_common){
            minor_allele_count += sample.reads[i];
        }
        depth += sample.reads[i];
    }
    return SiteSummary(most_common, minor_allele_count/(double)depth);
}

bool is_poly (ReadDataVector sample_data){
    auto it = sample_data.begin();
    auto first_sample = summarize_site(*it);
    uint16_t major_allele = first_sample.first;
    if( first_sample.second > 0.05 ){
        return true;
    }
    for(; it != sample_data.end(); it++){
       SiteSummary summarized = summarize_site(*it);
        if( summarized.first != major_allele){
            return true;
        }
        if( summarized.second > 0.05){
            return true;
        }
    }
    return false;
}

#endif 








