#This program is used to generate a .bedgraph file of genomic coordinates with the number of times a base is covered in the 4th column
#The input .sam file and output .bedgraph file are provided as command-line arguments
#The program assumes that the only mapping locations are yeast chromosomes, labeled "I", "II", ..., "XIV", "MT"
#Will only work with paired-end reads
#Example call: python3 /home/ebondra/code/coverage.py /home/ebondra/EXP_070/sam/a.sam /home/ebondra/EXP_070/bedgraphs/a.bedgraph 0 500

import sys

file_in=open(sys.argv[1], 'r')
file_out=open(sys.argv[2], 'w')
min_length=int(sys.argv[3])
max_length=int(sys.argv[4])

lines = file_in.readlines()

#Generate a dictionary with each chromosome as a key and a subdictionary as the values
output_dict = {"I":{},"II":{},"III":{},"IV":{},"V":{},"VI":{},"VII":{},"VIII":{},"IX":{},
               "X":{},"XI":{},"XII":{},"XIII":{},"XIV":{},"XV":{},"XVI":{},"MT":{}}


#reads in the .sam file one line at a time                                                                                                                             
for line in lines:
        #Only reads lines after top material (start with @) and whose length is within the range specified above by min_length and max_length                          
        if not(line.split()[0].startswith("@")) and abs(int(line.split()[8]))>min_length and abs(int(line.split()[8]))<max_length:
                chrom = line.split()[2]
                start = int(line.split()[3])
                mate_start = int(line.split()[7])
                length = abs(int(line.split()[8]))
                end = min(start,mate_start)+length
                #Assumes paired-end reads--finds lower genomic coordinate and sets as starting point                                                                   
                #for adding 1 to the pileup and adds 1 at each point between the lower genomic coordinate                                                              
                #and the lower genomic coordinate + the length of the insert. Only adds 0.5 at each line                                                               
                # of the sam file since it will read both mates individually                                                                                           
                for i in range(min(start,mate_start),end):
                        if i in output_dict[chrom]:
                                output_dict[chrom][i] = output_dict[chrom][i]+0.5
                        else:
                                output_dict[chrom][i]=0.5

#The only feature in the output file besides the data is a specification that the type of file is "bedGraph"                                                           
file_out.write("type=bedGraph\n")
#Reads each item in the sorted main dictionary (of chromosomes) and then the sorted subdictionary (of coordinates)                                                     
for x in sorted(output_dict):
    for y in sorted(output_dict[x]):
        #creates a tab-delimited line with four columns: chromosome, location, location+1 (i.e. 0-width), and number of times that location is covered                 
        lineout = x + "\t" + str(y)+ "\t" + str(y+1) + "\t" + str(output_dict[x][y]) + "\n"
        file_out.write(lineout)
file_in.close()
file_out.close()
