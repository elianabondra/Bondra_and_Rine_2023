#This code converts .sam to .bam files - Files listed here are an example.
bamdir=/home/ebondra/EXP_070/bam
for i in 380_IP 380_IN 77_IP 77_IN E1-IP E-IN E2-IP S1-IP S2-IP S-IN
do

    samtools view -b -S -o $bamdir/${i}.bam /home/ebondra/EXP_070/sam/${i}.sam
    samtools sort $bamdir/${i}.bam -o $bamdir/${i}-sorted.bam
    samtools index $bamdir/${i}-sorted.bam

done
