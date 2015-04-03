## Please use at your own discretion. test it on a small set. i have tested it thoroughly but you do it as well.
## this script creates pseudogenomes from a reference by using a tab delimited snp matrix
## input snp matrix format should be tab separated and in the following format
##chrom    postion    reference    sample1    sample2    sample3
##  no other extra colum should be in the file
## if you already have a VCF snp file for your samples please use the GATK reference maker instead
## java -Xmx2g -jar GenomeAnalysisTK.jar \
##   -R ref.fasta \
##   -T FastaAlternateReferenceMaker \
##   -o output.fasta \
##   -L input.intervals \
##   --variant input.vcf \
##   [--snpmask mask.vcf]

import sys,re,collections 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pprint import pprint
from argparse import ArgumentParser
from collections import defaultdict


def rreplace(s, old, new, occurrence): #taken from stackoverflow
    li = s.rsplit(old, occurrence)
    return new.join(li)

def find(s, ch): #taken from stackoverflow
    return [i for i, ltr in enumerate(s) if ltr == ch]
 
def arg_parse():
    parser=ArgumentParser(description="replaces sites in reference fasta with SNPs in SNP file")
    parser.add_argument("-i","--input_SNP_file",type=str,help="tab separated , chrom name should exactly match fasta header , no other extra colum should be in the file, the thrid colummn should be reference nucleotides")
    parser.add_argument("-r","--input_reference_file",type=str,help="input reference file")
    return parser

def main():
    parser=arg_parse()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    arg=parser.parse_args()
    #print "in"    
    
    snp_dict = defaultdict(dict)
  
    
    reference_dict = SeqIO.index(arg.input_reference_file, "fasta") # read in your reference that is to be changed to the per sample pseudogenome
   
        
    with open(arg.input_SNP_file,"r") as infile: # open snp matrix file
        header_str=""
        header_arr=[]
        line_arr=[]
        header_str=infile.readline() # read in the header line of snp matrix file
        #print header_str 
        header_arr=header_str.rstrip().split("\t") # split the header line on tabs
        #pprint(header_arr)
        header_arr.pop(0) # remove 'chrom' from header array
        header_arr.pop(0) # remove 'position' from header array
        #pprint(header_arr)
        
       # for itm in header_arr:
       #     snp_dict[itm.strip()]=""
        for line in infile: # loop and read the rest of the rows of snp matrix file
        # print(line)  
         if line.strip(): 
          line_arr=line.rstrip().split('\t')  # split the line on tabs
          #pprint(line_arr)
          count=1 
           # we will make use of preserved order of a python list here.
           #each column we encounter in line_arr is ordered like header_arr. line_arr[0] is the chrom name , line_arr[1] is the position
          for sample in header_arr: 
             count+=1
             if sample:
                # we are trying to create a multi dimensional array  # snp_dict[sample][chrom]=[position,SNP]
                if line_arr[0] not in snp_dict[sample]: # if chromosome name is not in dictionary add it and associate a list with it
                  snp_dict[sample][line_arr[0]]=[]
                snp_dict[sample][line_arr[0]].append([int(line_arr[1]),line_arr[count]])
                 
        # we have a snp dictionary now. time to replace the reference postions with the snps and create pseudo genomes         
        for sample in header_arr:
         outfl=sample+"_pseudo.fasta" # for each sample in the header array create a sample_pseudo.fasta file
         with open(outfl,"w") as ofl:
             chrom_rec=snp_dict[sample] # now we loop through the ssnp dictionary we created earlier. this step gets the chromosome dict for each sample
             
             for chrom_nm in chrom_rec:
                # print chrom_nm
                 seq_record=reference_dict[chrom_nm] # extract the chromosome record from the reference file
                 my_id=""
                 my_seq=""
                 my_seq_lst=[]
                 my_id=seq_record.id
                 my_seq=seq_record.seq
            
                 for pos_snp_lst in chrom_rec[chrom_nm]: # extract a list of [ postion and snp] for each chrom
                  #pprint(pos_snp_lst)
                  my_seq = my_seq[:pos_snp_lst[0]-1] + pos_snp_lst[1] + my_seq[pos_snp_lst[0]:] # replace the reference postion with SNP
                
                 
                 #write out the pesudogenome
                 ofl.write(">%s ref_%s\n%s\n" % (
                          sample,
                          seq_record.id,
                          my_seq
                          ))
    
                
if __name__=="__main__":
    main() 
