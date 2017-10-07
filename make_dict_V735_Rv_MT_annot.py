#!/usr/bin/python3

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio import Seq
import subprocess
import fileinput

def main():
    map_file = open("V735_to_Rv_MT_annot_mapping.csv",'w') # File to store the mapping
    #return_code = subprocess.call("blast_bin/blastp -out blast_results.xml -outfmt 5 -query CDC1551_proteins.fa -db old_CDC1551_blast_db/old_CDC1551 -evalue .001 -max_hsps 1 -max_target_seqs 50", shell=True)
    V735_records = SeqIO.parse("CDC1551_genes.fa","fasta")
    result_handle = open("blast_results.xml")
    blast_records = NCBIXML.parse(result_handle)
    V735_to_MT_dict = {} # Initialize dictionary to hold mapping
    V735_to_Rv_dict = {} # Initialize dictionary to hold mapping
    V735_to_annotation_dict = {} # Initialize dictionary to hold mapping
    i = -1
    for blast_record in blast_records:
        V735_record = next(V735_records)
        i = i+1
        #print(i)
        V735_num, Rv_num, annotation = V735_record.description.split(sep = '|', maxsplit=3)[0:3]
        annotation = annotation[1:].replace(',',';')
        V735_to_Rv_dict[V735_num] = Rv_num
        V735_to_annotation_dict[V735_num] = annotation
        if blast_record.alignments:
            evalue = ""
            MT_num = ""
            for alignment in blast_record.alignments:
                evalue = evalue + str(alignment.hsps[0].expect) + ' '
                MT_num = MT_num + alignment.title.split(sep = '|', maxsplit=2)[0] + ' '
            map_file.write(V735_num + ',' + Rv_num + ',' + MT_num + ',' + str(evalue) + ',' + annotation + '\n')
            V735_to_MT_dict[V735_num] = MT_num
        else:
            map_file.write(V735_num + ',' + Rv_num + ',' + 'None' + ',' + 'NA' + ',' + annotation + '\n')

if __name__ == "__main__":
    main()
