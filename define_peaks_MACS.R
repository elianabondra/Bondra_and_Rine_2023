samtools merge merged_IP.bam 237_00_IP-sorted.bam 238_00_IP-sorted.bam 

macs3 callpeak -t merged_IP.bam -c 12923_IP-sorted.bam -n merged_IP_over_notag -f BAMPE -g 2.4e7 -q 0.01 --keep-dup=auto -B --call-summits

bedtools slop -i merged_IP_over_notag_summits.bed -g yeast-sizes.txt -b 150 > peaks_over_notag.bed

awk -F\| '{if($1~/[a-zA-Z]+_cer/){print}}' peaks_over_notag.bed > cer_peaks_notag.bed  
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' cer_peaks_notag.bed > cer_notag_featurecounts.saf
