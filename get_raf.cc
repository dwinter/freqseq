#include <iostream>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <sys/stat.h>
#include <unordered_map>

#include "utils.h"
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"


using namespace BamTools;
using namespace std;

typedef unordered_map<string, uint16_t> SampleMap;

uint16_t base_index(char b){                                                                        
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

class RafVisitor: public PileupVisitor{
    public:
        RafVisitor(const RefVector& bam_references,
                    BamAlignment& ali, 
                    const Fasta& ref_genome,
                    const SamHeader& header,
                    const SampleMap& samples,
                    ostream *out_stream
                    ):
            PileupVisitor(), m_samples(samples), 
                             m_ref_genome(ref_genome),
                             m_header(header),
                             m_bam_ref(bam_references),
                             m_ostream(out_stream)


            {             

                nsamp = m_samples.size(); 
            }
        ~RafVisitor(void) { }
    public:
        void Visit(const PileupPosition& pileupData){
            uint64_t pos  = pileupData.Position;                                                   
            m_ref_genome.GetBase(pileupData.RefId, pos, current_base);
            uint16_t ref_base_idx = base_index(current_base);
            vector<uint32_t> nref (nsamp, 0);
            vector<uint32_t> depth (nsamp, 0);
            for(auto it = begin(pileupData.PileupAlignments);
                     it !=end(pileupData.PileupAlignments);
                     ++it){
                if( passes_QC(*it, 30, 13) ){
                    
                    uint16_t b_index = get_base_index(it->Alignment.QueryBases[it->PositionInAlignment]);
                    if (b_index < 4){
                        string tag_id;
                        it->Alignment.GetTag("RG", tag_id);
                        uint32_t sindex = m_samples[tag_id];
                        depth[sindex] += 1;
                        if(b_index == ref_base_idx){
                            nref[sindex] += 1;
                        }
                    }
                }
            }
            for (size_t i = 0; i< nsamp; i++){
                if(depth[i] > 0){
                    *m_ostream << m_bam_ref[pileupData.RefId].RefName << '\t' 
                               << pos << '\t' << i << '\t' << nref[i]  <<'\t' << depth[i] 
                               << endl;
                }

            }
    }                          
    private:
        RefVector m_bam_ref;
        Fasta m_ref_genome;
        ostream* m_ostream;
        int nsamp;
        SampleMap m_samples;
        SamHeader m_header;
        int m_ref_pos;
        char current_base;
};
        
int main(int argc, char* argv[]){
    string bam_path = argv[1];
    string ref_file = argv[2];
    string out_prefix = argv[3];
    ofstream result_stream (out_prefix + "_bases.tsv");
    
    //Set up the Bam for parsing...
    BamReader bam;
    bam.Open(bam_path);
    bam.OpenIndex(bam_path + ".bai");
    RefVector references = bam.GetReferenceData(); 
    SamHeader header = bam.GetHeader();

    //Now the reference genome, creating index if there isn't already one
    Fasta reference_genome;
    struct stat file_info;
    string faidx_path = ref_file + ".fai";
    if (stat(faidx_path.c_str(), &file_info) !=0){
        reference_genome.Open(ref_file);
        reference_genome.CreateIndex(faidx_path);
    }
    reference_genome.Open(ref_file, faidx_path);

    //Map readgroups to samples
    SampleMap name_map; 
    uint16_t sindex = 0;
    for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
        if(it->HasSample()){
            cerr << it->Sample <<endl;
            auto s  = name_map.find(it->Sample);
            if( s == name_map.end()){ // not in there yet
                name_map[it->Sample] = sindex;
                sindex += 1;
            }
        }
    }
    // And now, go back over the read groups to map RG:sample index
    SampleMap samples;
    for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
        if(it->HasSample()){
            samples[it->ID] = name_map[it->Sample];  
        }
    }

    uint32_t nsamp = samples.size();

    BamAlignment ali;
    PileupEngine pileup;

    RafVisitor *f = new RafVisitor(references, 
                                    ali, reference_genome,
                                    header,
                                    samples,
                                    &result_stream);
    pileup.AddVisitor(f); 
    while(bam.GetNextAlignment(ali)){
        pileup.AddAlignment(ali);
    }
    pileup.Flush();
    
    return 0;
}      

