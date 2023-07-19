# code
code_for_Bondra-and-Rine_PNAS_2023 "Context dependent function of the transcriptional regulator Rap1 in gene silencing and activation in Saccharomyces cerevisiae"

The source codes in this folder are for the analysis and plotting used in the manuscript above and NCBI GEO Series GSE227763
Included in this folder are the custom genomes to which the raw sequencing data were aligned.

For data associated with figures 1-3, the general workflow is alignment_example >> rmdup_folder >> sort_index >> coverage >> normalize_bedgraphs >> example_bedgraph_plot


For data associated with figures 4-5, a similar workflow takes place. However, instead "AA_genome_stats" was used to calculate the read coverage assigned to the paradoxus genome. 
