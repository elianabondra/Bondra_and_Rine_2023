#Assumes that all sam files are in the same folder and all end with ".sam"     \
                                                                                
#Change directory to that provided as argument                                 
cd "$1"

#Iterates through all sam files in folder                                     
for file in *.sam
do
        output=${file%.sam}_fix.sam
        samtools fixmate -@13 -m "$file" "$output"
done

for file in *_fix.sam
do
        output=${file%.sam}_sorted.sam
        samtools sort -@13 -o "$output" "$file"
done

for file in *sorted.sam
do
        output=${file%.sam}_rmdup.sam
        samtools markdup -r -@13 "$file" "$output"
done
