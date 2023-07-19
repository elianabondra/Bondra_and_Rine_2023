# code
code_for_Bondra-and-Rine_PNAS_2023 "Context dependent function of the transcriptional regulator Rap1 in gene silencing and activation in Saccharomyces cerevisiae"

The source codes in this folder are for the analysis and plotting used in the manuscript above and NCBI GEO Series GSE227763
Included in this folder are the custom genomes to which the raw sequencing data were aligned.

For data associated with figures 1-3, the general workflow is alignment_example >> rmdup_folder >> sort_index >> coverage >> normalize_bedgraphs >> example_bedgraph_plot


For data associated with figures 4-5, a similar workflow takes place through alignment (to hybrid_genome_forAA_fig4-5.fa). Then "AA_genome_stats" was used to calculate the read coverage assigned to the paradoxus genome. 

To define peak regions, "define_peaks_MACS.R" was called using a merged JRY15237_00_IP and JRY15238_00_IP as the IP, and a no-tag control JRY12923_IP as the Input. 


