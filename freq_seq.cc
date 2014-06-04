#include <iostream>
#include <stdint.h>
#include <string>
#include <algorithm>

#include "utils.h"
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"


using namespace BamTools;
using namespace std;




class FreqVisitor: public PileupVisitor{
    public:
        FreqVisitor(const RefVector& bam_references,
                    BamAlignment& ali, 
                    const SamHeader& header,
                    const SampleNames& samples
                    ):
            PileupVisitor(), m_samples(samples), 
                             m_header(header),
                             m_bam_ref(bam_references) {
                nsamp = m_samples.size(); 
            }
        ~FreqVisitor(void) { }
    public:
        void Visit(const PileupPosition& pileupData){
            ReadDataVector sample_reads;
            for(auto it = begin(pileupData.PileupAlignments);
                     it !=end(pileupData.PileupAlignments);
                     ++it){
                if( passes_QC(*it, 13, 13) ){
                    uint16_t b_index = get_base_index(it->Alignment.QueryBases[it->PositionInAlignment]);
                    if (b_index < 4){
                        string tag_id;
                        it->Alignment.GetTag("RG", tag_id);
                        string sm = m_header.ReadGroups[tag_id].Sample;
                        uint32_t sindex = find_sample_index(sm, m_samples);
                        sample_reads[sindex].reads[b_index] += 1;
                        }
                    }
                }
                if (is_poly(sample_reads)){
                    string chr = m_bam_ref[pileupData.RefId].RefName;
                    uint64_t pos  = pileupData.Position;
                    ReadData overall = {0,0,0,0};
                    //find major allele over all
                    for (auto summary_it =  begin(sample_reads); 
                              summary_it != end(sample_reads);
                              summary_it++){
                        for(size_t i = 0; i < 4; i++){
                            overall.reads[i] += summary_it->reads[i];
                        }
                    }
                    uint16_t major_allele = distance(overall.reads, 
                                             max_element(overall.reads, overall.reads+4));
    
                    cout << chr << '\t' <<  pos << '\t' << major_allele << '\t';
                    //find frequency of non-major alleles in all sites
                    for (auto summary_it =  begin(sample_reads); 
                              summary_it != end(sample_reads);
                              summary_it++){
                        uint64_t depth = 0;
                        uint64_t minor_allele_count = 0;
                        for (size_t i=0; i < 4; i++){
                            if (i != major_allele){
                                minor_allele_count += summary_it->reads[i];
                            }
                            depth += summary_it->reads[i];
                        }
                        cout << minor_allele_count / (double)depth << '\t';
                    }
                    cout << endl;
                }
            }
    private:
        RefVector m_bam_ref;
        int nsamp;
        SampleNames m_samples;
        SamHeader m_header;
        int m_ref_pos;
        string m_initial_data;
};
        
int main(int argc, char* argv[]){
    string bam_path = argv[1];
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
    BamAlignment ali;
    PileupEngine pileup;

    FreqVisitor *f = new FreqVisitor(references, 
                                    ali, 
                                    header,
                                    samples);
    pileup.AddVisitor(f); 
    while(bam.GetNextAlignment(ali)){
        pileup.AddAlignment(ali);
    }
    pileup.Flush();
    return 0;
}
      
