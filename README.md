
The files here are supporting information for Bondra and Rine 2023, "Context dependent function of the transcriptional regulator Rap1 in gene silencing and activation in Saccharomyces cerevisiae"

The source codes in this folder are for the analysis and plotting used in the manuscript above and NCBI GEO Series GSE227763
Included in this folder are the custom genomes to which the raw sequencing data were aligned. 

For data associated with figures 1-3, the general workflow is alignment_example (aligned to "fig1-3_genome.fa") >> rmdup_folder >> sort_index >> coverage >> normalize_bedgraphs >> example_bedgraph_plot


For data associated with figures 4-5, a similar workflow takes place through alignment (to hybrid_genome_forAA_fig4-5.fa). Then "AA_genome_stats" was used to calculate the read coverage assigned to the paradoxus genome. 

To define peak regions, "define_peaks_MACS.R" was called using a merged JRY15237_00_IP and JRY15238_00_IP as the IP, and a no-tag control JRY12923_IP as the Input. 

Peaks were assigned to ORFs using "peak_assignment" and subsetted into the 377 loci used in the paper using "subset peaks". 

Bedgraphs for Figure 4B were plotted using "HML_bedgraphs.R" and "MAT_bedgraphs.R". 
Fitting the coverage data from SI Appendix Dataset S5, with normalized peak coverage over time to the non-linear regression as outlined in Materials and Methods, was done using "final_fits.R" 
Non-linear regressions for Figure 4C,D were plotted using "decay_plots.R". 
Figure 5B was plotted using "distance_vs_restime.R"
"AA_plots.R" contains code for data analysis and plotting used in Figure 5, Figure S5, and Figure S6. 
