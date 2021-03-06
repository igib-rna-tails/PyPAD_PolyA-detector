import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import csv
import re
from Bio import SeqIO
from Bio.Seq import MutableSeq
import filecmp
import subprocess
import os
import subprocess 
import argparse

pattern_PC6 = re.compile("PC")
pattern_6AU_forward = re.compile("[AT]{6}$")
pattern_6AU_reverse = re.compile("^[AT]{6}")
pattern_1AU_forward = re.compile("[AT]{1}$")
pattern_1AU_reverse = re.compile("^[AT]{1}")
pattern_1A_forward = re.compile("[A]{1}$")
pattern_1A_reverse = re.compile("^[A]{1}")

pattern_min6AU_forward =  re.compile("[AT]{6,}$")
pattern_min6AU_reverse =  re.compile("^[AT]{6,}")





## infile path == unmapped reads fastq
def selectminNnucleotides_trimby1nt(strandness = 'forward', selectminNnucteotides = True,
                                    mintail = pattern_min6AU_forward, pattern_loop = pattern_1AU_forward,
                                    infile_path = "unmapped_reads.fastq",
                           outfile_path = 'output.fastq'):
    '''
    Function can select reads having at least 6AU at 3' tail of RNA
    
    '''
    def take_reads_corresponding_to_headerlist(header_list_path = 'lista_headerow.txt', infile_path =  'unmapped_reads.fastq',outfile_path='unmapped_reads2.fastq',pattern=pattern_6AU_forward):
        with open(header_list_path) as f:
            pojedyncze_headery = f.readlines()
            pojedyncze_headery = [x.strip() for x in pojedyncze_headery]
            content2 = [x.split('_')[0] for x in pojedyncze_headery]
            rodzaje_headerow = set(content2)

        with open(outfile_path, 'w') as w:
            infile = SeqIO.parse(infile_path, "fastq")
            for read in infile:
                fastq = read.format("fastq")
                header = fastq.partition('\n')[0]
                sequence = fastq.splitlines()[1]
                strand = fastq.splitlines()[2]
                quality = fastq.splitlines()[3]
                ogon = re.search(pattern, str(sequence))
                if header in rodzaje_headerow:
                    header_new = header
                else:
                    header_new = 'unmappedreads'
                new_read =  header_new + '\n' + sequence + '\n' + strand + '\n'+ quality + '\n'
                w.write(new_read)

        os.system("grep -A 3 -E '^@' unmapped_reads2.fastq | sed '/--/d' > unmapped_reads.fastq")
    def move_pattern_to_header(infile_path, outfile_path, pattern):
        '''
        Function removes the pattern from the end of second line of fastq read, and move it to the header of read
        in the form:  @header_pattern.
        The function check length of the read and remove appropriate number of characters from fourth line of read
        (the line with the quality of read). 
        Reads without pattern have header with "_withoutanytail" at the end of the end of header.
        
        Input:
            infile_path: str
            outfile_path: str
            pattern: str
            
        Output: new fastq file (under output_path) with pattern (tail) mover from the sequence to header, 
                and truncated qualit line. 
        '''
        with open(outfile_path, 'w') as w:
            infile = SeqIO.parse(infile_path, "fastq")
            for read in infile:
                fastq = read.format("fastq")
                header = fastq.partition('\n')[0]
                sequence = fastq.splitlines()[1]
                strand = fastq.splitlines()[2]
                quality = fastq.splitlines()[3]
                pattern_str = re.compile(pattern)
                ogon = re.search(pattern_str, str(sequence))
                if ogon == None:
                    ogon = "withoutanytail"
                else:
                    ogon = ogon.group(0)
                    ogon_length = len(ogon)
                    sequence = str(sequence[0:(0-ogon_length)])
                    seq_length = len(sequence)
                    quality = str(quality[0:(0-ogon_length)])
                header=str(read.id) + "_" + ogon 
                new_read = '@'+header+'\n'+sequence+'\n'+strand+'\n'+ quality+'\n'
                w.write(new_read)
                
                
    pattern_1AT_forward =  re.compile("(A|T)$")
    pattern_1AT_reverse =  re.compile("^(A|T)")

    def rm_1nt_at_time_map_frombeginning(infile_path, 
                                         outfile_path = "unmapped_reads.fastq", 
                                         pattern=pattern_1AT_forward):
        '''
        Function cut one nucleotide A or T from the 3' end of RNA (end of the sequence) and move it to the header.
        It cut also one character from quality line in fastq
        
        
        '''
        move_pattern_to_header(infile_path, outfile_path = "unmapped_reads1.fastq", 
                                         pattern=pattern_1AT_forward)    
      #  os.system("grep -A 3 -E '_(A|T)+$' unmapped_reads1.fastq | sed '/--/d' | sed -n '1~4s/^@/>/p;2~4p' sample_2.fasta")
        os.system("grep -A 3 -E '_(A|T)+$' unmapped_reads1.fastq | sed '/--/d' > sample_2.fastq")
        #mapowanie
        os.system("sed -n '1~4s/^@/>/p;2~4p' sample_2.fastq > sample_2.fasta")
        os.system("rm sample_2.fastq")


        os.system('hisat2 --score-min L,0,-0.4 -f -x reference/reference sample_2.fasta --rna-strandness FR -S sample_2.sam --summary-file sample_2.txt --new-summary -p 16')
        os.system("rm sample_2.fasta")
        os.system("rm sample_2.txt")

        os.system("samtools sort sample_2.sam > sample_2.bam")
        os.system("rm sample_2.sam")

        #wyci??gam niezmapowane i konwertuj?? do fastq
        os.system("samtools view -b -f 4  sample_2.bam > unmapped.bam")
        os.system("bedtools bamtofastq -i unmapped.bam -fq unmapped.fastq")
        
        os.system("rm unmapped.bam")


        #tylko zmapowane - zostawiam we wsp??lnym pliku
        os.system("samtools view -b -F 4 sample_2.bam > mapped2.bam")
        os.system("rm sample_2.bam")
        os.system("bedtools bamtofastq -i mapped2.bam -fq mapped2.fastq")
        os.system("rm mapped2.bam")
        if filecmp.cmp('mapped2.fastq', 'all_mapped_po1nt.fastq') ==False:
            os.system("cat all_mapped_po1nt.fastq mapped2.fastq >> all_mapped_po1nt_before.fastq")
            os.system("cp all_mapped_po1nt_before.fastq all_mapped_po1nt.fastq")
            os.system("rm all_mapped_po1nt_before.fastq")
      
        else:
            print('No reads mapped in this iteration')
        #!rm unmapped_sample.fastq
        os.system("rm mapped2.fastq")
        
        if filecmp.cmp('unmapped.fastq', 'empty.fastq') == False:
            rm_1nt_at_time_map_frombeginning(infile_path = "unmapped.fastq", outfile_path = "unmapped_reads.fastq",
                              pattern=pattern_1AT_forward)
        
        else:
            
            print('Happy end')
            
    def rm_1nt_at_time_map_frombeginning_withoutpreselection(infile_path, 
                                         outfile_path = "unmapped_reads.fastq", 
                                         pattern_loop=pattern_1AT_forward):
        '''
        Function cut one nucleotide A or T from the 3' end of RNA (end of the sequence) and move it to the header.
        It cut also one character from quality line in fastq
        
        
        '''
        move_pattern_to_header(infile_path, outfile_path = "unmapped_reads1.fastq", 
                                         pattern=pattern_1AT_forward)    
        os.system("grep -A 3 -E '_(A|T)+$' unmapped_reads1.fastq | sed '/--/d' > sample_2.fastq")
        #mapowanie
        #!sed -n '1~4s/^@/>/p;2~4p' sample_2.fastq > sample_2.fasta
        
        os.system("hisat2 -q -x reference/reference sample_2.fastq --rna-strandness FR -S sample_2.sam --summary-file sample_2.txt --new-summary -p 16")
        os.system("rm sample_2.fastq")
        os.system("rm sample_2.txt")
        os.system("samtools sort sample_2.sam > sample_2.bam")
        os.system("rm sample_2.sam")
        #wyci??gam niezmapowane i konwertuj?? do fastq
        os.system("samtools view -b -f 4  sample_2.bam > unmapped.bam")
        os.system("bedtools bamtofastq -i unmapped.bam -fq unmapped.fastq")
        
        os.system("rm unmapped.bam")


        #tylko zmapowane - zostawiam we wsp??lnym pliku
        os.system("samtools view -b -F 4 sample_2.bam > mapped2.bam")
        os.system("bedtools bamtofastq -i mapped2.bam -fq mapped2.fastq")
        if filecmp.cmp('mapped2.fastq', 'all_mapped_po1nt.fastq') ==False:
            os.system("cat all_mapped_po1nt.fastq mapped2.fastq >> all_mapped_po1nt_before.fastq")
            os.system("cp all_mapped_po1nt_before.fastq all_mapped_po1nt.fastq")
            os.system("rm all_mapped_po1nt_before.fastq")
        else:
            print('No reads mapped in this iteration')
        os.system("rm mapped2.bam")
        #!rm unmapped_sample.fastq
        os.system("rm sample_2.bam")
        os.system("rm mapped2.fastq")
        
        if filecmp.cmp('unmapped.fastq', 'empty.fastq') == False:
            rm_1nt_at_time_map_frombeginning(infile_path = "unmapped.fastq", outfile_path = "unmapped_reads.fastq",
                              pattern=pattern_1AT_forward)
        
        else:
            
            print('Happy end')

    def organize_tail_after_loop(infile_path, outfile_path):

        with open(outfile_path, 'w') as w:
            infile = SeqIO.parse(infile_path, "fastq")
            for read in infile:
                fastq = read.format("fastq")
                header = fastq.partition('\n')[0]
                sequence = fastq.splitlines()[1]
                strand = fastq.splitlines()[2]
                quality = fastq.splitlines()[3]
                pattern_string = str(header)[1:4]
                pattern_str = re.compile(pattern_string)#"PC6"
                ogon = re.search(pattern_str, str(sequence))

                headerek = header.split("_")
                headerek2 = headerek[1:]
                nazwa_readu = headerek[0]
                headerek2.reverse()
                header_str = ''.join(headerek2)
                new_read = nazwa_readu + "_" + header_str + '\n' + sequence + '\n'+ strand + '\n' + quality +'\n'
                w.write(new_read)
    
    
    os.system("touch all_mapped_po1nt.fastq")
    os.system("touch empty.fastq")
    if __name__ == '__main__':
        parser = argparse.ArgumentParser(description='Welcome to PyPAD - a tool that detects polyadenylation in the available RNA-seq sequencing data . Written by lipinska@biol.uw.edu.pl and maintained at https://github.com/igib-rna-tails/PyPAD_PolyA-detector.')
        parser.add_argument('--strandness', type=str, choices=['forward', 'reverse'], default="forward", help="Forward reads (R1, reading from 5' to 3') or reverse reads (R2, reading from 3' to 5'")
        parser.add_argument('--selectminNnucteotides',  default=True, help='Option in PyPAD to pre-select reads having eg 6 nt in the tail before the proceduce of triming one by one nucleotide from the tail, and realign fastq file .')
        parser.add_argument('--mintail', choices=[pattern_min6AU_forward, pattern_min6AU_reverse], default=pattern_min6AU_forward, help='Option of pattern to preselect by --selectminNnucteotides.')
        parser.add_argument('--pattern_loop',choices=[pattern_1AU_forward, pattern_1AU_reverse, pattern_1A_forward, pattern_1A_reverse ],  default=pattern_1AU_forward, help="Pattern of nucleotides cut one by one from 3' tail")
        parser.add_argument('--infile_path', type=str, default="unmapped_reads.fastq", help='Path to the fastq file with unmapped reads to analyse. Please prepared data before the analysis with PyPAD')
        parser.add_argument('--outfile_path', type=str, default="output.fastq", help='Path to the output fastq file after PyPAD analysis. The fastq file contains reads with precise polyA tail sequence in the header.')
        args = parser.parse_args()
    if strandness == 'forward':
        if selectminNnucteotides == True:
            move_pattern_to_header(infile_path=infile_path, outfile_path='pattern_in_header.fastq', 
                                   pattern = mintail)
            os.system("grep -A 3 -E '_A|T{6,}$' pattern_in_header.fastq | sed '/--/d' > firststep2.fastq")
            os.system("sed -n '1~4s/^@/>/p;2~4p' firststep2.fastq > firststep2.fasta")
            os.system("hisat2  -f -x reference/reference firststep2.fasta --rna-strandness FR -S firststep2.sam --summary-file firststep2.txt --new-summary -p 16")

            os.system("samtools sort firststep2.sam > firststep2.bam")
            os.system("samtools view -b -F 4 firststep2.bam > firststep3.bam")
            os.system("bedtools bamtofastq -i firststep3.bam -fq mapped.fastq")
            os.system("grep -E '_A|T{6,}$' mapped.fastq > lista_headerow.txt")
            
            os.system("rm mapped.fastq")
            take_reads_corresponding_to_headerlist(header_list_path = 'lista_headerow.txt',
                                           infile_path= infile_path, #unmapped reads 
                                           outfile_path='unmapped_reads2.fastq',
                                                  pattern=  pattern_6AU_forward)
            rm_1nt_at_time_map_frombeginning(infile_path='unmapped_reads.fastq',
                                             outfile_path='all_mapped_po1nt.fastq', 
                                             pattern=pattern_1AT_forward)
            os.system("cp all_mapped_po1nt.fastq output_headertoorganise.fastq")
            os.system("rm all_mapped_po1nt.fastq")
            organize_tail_after_loop(infile_path='output_headertoorganise.fastq',
                                    outfile_path='output_tails_in_header.fastq')
                                                   
            os.system("rm firststep*")
            

        elif selectminNnucteotides == False:
            rm_1nt_at_time_map_frombeginning_withoutpreselection(infile_path='unmapped_reads.fastq',
                                             outfile_path='all_mapped_po1nt.fastq', 
                                             pattern_loop=pattern_1AT_forward)
            os.system("cp all_mapped_po1nt.fastq output_headertoorganise.fastq")
            os.system("rm all_mapped_po1nt.fastq")
            organize_tail_after_loop(infile_path='output_headertoorganise.fastq',
                                    outfile_path='output_tails_in_header.fastq')
            
        else:
            print("Incorrect value. Please select True of False" )
            
            
            
            
    elif strandness == 'reverse':
        if selectminNnucteotides == True:
                
                
            move_pattern_to_header(infile_path=infile_path, outfile_path='pattern_in_header.fastq', 
                                   pattern = pattern_min6AU_reverse)
            os.system("grep -A 3 -E '_A|T{6,}$' pattern_in_header.fastq | sed '/--/d' > firststep2.fastq")
            os.system("sed -n '1~4s/^@/>/p;2~4p' firststep2.fastq > firststep2.fasta")
            os.system("hisat2  -f -x reference/reference firststep2.fasta --rna-strandness R -S firststep2.sam --summary-file firststep2.txt --new-summary -p 16")

            os.system("samtools sort firststep2.sam > firststep2.bam")
            os.system("samtools view -b -F 4 firststep2.bam > firststep3.bam")
            os.system("bedtools bamtofastq -i firststep3.bam -fq mapped.fastq")
            os.system("grep -E '_A|T{6,}$' mapped.fastq > lista_headerow.txt")
            os.system("rm firststep*")
            os.system("rm mapped.fastq")
            take_reads_corresponding_to_headerlist(header_list_path = 'lista_headerow.txt',
                                           infile_path = infile_path, #unmapped reads 
                                           outfile_path ='unmapped_reads.fastq',
                                                  pattern =  pattern_6AU_reverse)
            rm_1nt_at_time_map_frombeginning(infile_path ='unmapped_reads.fastq',
                                             outfile_path ='all_mapped_po1nt.fastq', 
                                             pattern = pattern_1AU_reverse)
            os.system("cp all_mapped_po1nt.fastq output_headertoorganise.fastq")
            os.system("rm all_mapped_po1nt.fastq")
            organize_tail_after_loop(infile_path ='output_headertoorganise.fastq',
                                    outfile_path ='output.fastq')
                                                   

            

        elif selectminNnucteotides == False:
            rm_1nt_at_time_map_frombeginning_withoutpreselection(infile_path=infile_path,
                                             outfile_path='all_mapped_po1nt.fastq', 
                                             pattern_loop=pattern_1AU_reverse)
            os.system("cp all_mapped_po1nt.fastq output_headertoorganise.fastq")
            os.system("rm all_mapped_po1nt.fastq")
            organize_tail_after_loop(infile_path='output_headertoorganise.fastq',
                                    outfile_path='output.fastq')  
        else: 
            print("Incorrect value. Please select True of False" )   
    else:
        print("Incorrect strandness. Please select 'forward' (for Read 1) or 'reverse' (for Read 2)" )
    os.system("rm unmapped_reads1.fastq")
    os.system("rm unmapped_reads2.fastq")
    os.system("rm unmapped.fastq")
    os.system("rm pattern_in_header.fastq")
    os.system("rm output_headertoorganise.fastq")
    os.system("rm lista_headerow.txt")
    


##strandness = 'forward', selectminNnucteotides = True,   mintail = pattern_min6AU_forward, pattern_loop = pattern_1AU_forward,
                                   ## infile_path = "unmapped_reads.fastq",outfile_path = 'output.fastq'
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Welcome to PyPAD - a tool that detects polyadenylation in the available RNA-seq sequencing data . Written by lipinska@biol.uw.edu.pl and maintained at https://github.com/igib-rna-tails/PyPAD_PolyA-detector.')
    parser.add_argument('--strandness', type=str, choices=['forward', 'reverse'], default="forward", help="Forward reads (R1, reading from 5' to 3') or reverse reads (R2, reading from 3' to 5'")
    parser.add_argument('--selectminNnucteotides',  default=True, help='Option in PyPAD to pre-select reads having eg 6 nt in the tail before the proceduce of triming one by one nucleotide from the tail, and realign fastq file .')
    
    
    parser.add_argument('--mintail', choices=[pattern_min6AU_forward, pattern_min6AU_reverse], default=pattern_min6AU_forward, help='Option of pattern to preselect by --selectminNnucteotides.')
    parser.add_argument('--pattern_loop',choices=[pattern_1AU_forward, pattern_1AU_reverse, pattern_1A_forward, pattern_1A_reverse ],  default=pattern_1AU_forward, help="Pattern of nucleotides cut one by one from 3' tail")
    
    parser.add_argument('--infile_path', type=str, default="unmapped_reads.fastq", help='Path to the fastq file with unmapped reads to analyse. Please prepared data before the analysis with PyPAD')
    parser.add_argument('--outfile_path', type=str, default="output.fastq", help='Path to the output fastq file after PyPAD analysis. The fastq file contains reads with precise polyA tail sequence in the header.')
    
    
    args = parser.parse_args()
   
    selectminNnucleotides_trimby1nt(args.strandness, args.selectminNnucteotides, args.mintail, args.pattern_loop, args.infile_path, args.outfile_path)

