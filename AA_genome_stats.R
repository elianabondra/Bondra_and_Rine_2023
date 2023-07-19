#This code takes sorted bam files and counts the number of reads per chromosome. These data were then categorized and used for normalization purposes for Anchor Away experiments
samples=/home/ebondra/EXP_046/samfiles/samples.txt

for i in $(cat ${samples})
do
    samtools idxstats /home/ebondra/EXP_046/bamfiles/$i-sorted.bam | cut -f 1,3 > /home/ebondra/EXP_046/genome_stats/$i.txt
done
