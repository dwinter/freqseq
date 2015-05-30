#include <iostream>
#include <stdint.h>
#include <string>
#include <algorithm>

#include "utils.h"
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"


using namespace BamTools;
using namespace std;

class RafVisitor: public PileupVisitor{
    public:
        RafVisitor(const RefVector& bam_references,
                    BamAlignment& ali, 
                    const Fasta& ref_genome,
                    const SamHeader& header,
                    const SampleMap samples,
                    vector<uint32_t>& nref,
                    vector<uint32_t>& depth,
                    ostream *out_stream
                    ):
            PileupVisitor(), m_samples(samples), 
                             m_ref_genome(ref_genome),
                             m_header(header),
                             m_bam_ref(bam_references),
                             m_ostream(out_stream),
                             m_nref(nref),
                             m_depth(depth)


            {             
                nsamp = m_samples.size(); 
            }
        ~RafVisitor(void) { }
    public:
        void Visit(const PileupPosition& pileupData){
            uint64_t pos  = pileupData.Position;                                                   
            m_ref_genome.GetBase(pileupData.RefId, pos, current_base);
            uint16_t ref_base_idx = base_index(current_base)
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
                *m_ostream << m_bam_ref[pileupData.RefId].RefName << '\t' 
                           << pos << '\t' << i << '\t' << nref[i] << depth[i] 
                           << endl;
            }
                              
    private:
        RefVector m_bam_ref;
        vector<uint32_t>& m_nref;
        vector<uint32_t>&m_depth;
        ostream* m_ostream;
        int nsamp;
        SampleMap m_samples;
        SamHeader m_header;
        int m_ref_pos;
        string m_initial_data;
        char current_base;
};
        
int main(int argc, char* argv[]){
    string bam_path = argv[1];
    string out_prefix = argv[2];
    ofstream result_stream (out_prefix + "_bases.tsv");
    

    BamReader bam;
    bam.Open(bam_path);
    bam.OpenIndex(bam_path + ".bai");
    SamHeader header = bam.GetHeader();
    SampleNames samples;
    RefVector references = bam.GetReferenceData(); 
    for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
        if(it->HasSample()){
                samples.push_back(it->Sample);
        }
    }
    uint32_t nsamp = samples.size();
    vector<uint32_t> minor_denom (nsamp, 0);
    vector<uint32_t> minor_count (nsamp, 0);

    BamAlignment ali;
    PileupEngine pileup;

    FreqVisitor *f = new FreqVisitor(references, 
                                    ali, 
                                    header,
                                    samples,
                                    minor_count,
                                    minor_denom,
                                    &result_stream);
    pileup.AddVisitor(f); 
    while(bam.GetNextAlignment(ali)){
        pileup.AddAlignment(ali);
    }
    pileup.Flush();
    
    ofstream summary (out_prefix + "_summary.tsv");
    for (size_t i=0; i < nsamp; i++){
        summary << samples[i] << '\t'
                << minor_count[i]/(double)minor_denom[i];
    }


    return 0;
}      
