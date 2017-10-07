#!/usr/bin/python3

import sys, argparse, os.path
import pdb
from Bio import SeqIO

def main():
# Process Arguments

    parser = argparse.ArgumentParser(description='Produces *.rnt and *.ptt files from a Genbank *.gbff file.')
    parser.add_argument('input_file', type=str, help="Input file from Genbank (*.gbff)")
    parser.add_argument('-r','--rnt', type=str, help="*.rnt filename for output", default='')
    parser.add_argument('-p','--ptt', type=str, help="*.ptt filename for output", default='')

    args = parser.parse_args()

# Create .rnt and .ptt file from .gbff file. 
# This code is a slightly modified version of its original found at: https://www.biostars.org/p/116540/ (submitted by user rafa.rios.50).
    annotation_file = args.input_file 
    rnt_file = args.rnt if len(args.rnt) > 0 else os.path.split(os.path.splitext(args.input_file)[0])[1] + '.rnt' 
    ptt_file = args.ptt if len(args.ptt) > 0 else os.path.split(os.path.splitext(args.input_file)[0])[1] + '.ptt' 
    print(annotation_file, rnt_file, ptt_file, sep='\n')
    r = SeqIO.parse(annotation_file, "gb")
    strand = {1:'+', -1:'-'}
    fout = open(rnt_file, "w")
    for record in r:
       record.features = [f for f in record.features if f.type == "rRNA" or f.type == "tRNA"]
       fout.write("{0} - 0..{1}\n".format(record.description, len(record)))
       fout.write("{0} RNAs\n".format(len(record.features)))
       fout.write("Location\tStrand\tLength\tPID\tGene\tSynonym Code\tCOG\tProduct\n")
       fout.write("{0}\n".format("\t".join(['0..0','+','0','-','EmptyRNA','EmptyRNA',"-",'EmptyRNA necessary to avoid EDGE-pro bug'])))
       for f in record.features:
           gene_location = str(f.location.start)+".."+str(f.location.end)
           gene_name = f.qualifiers["locus_tag"][0]
           synonym_code = f.qualifiers["locus_tag"][0]
           fout.write("{0}\n".format("\t".join([gene_location,strand[f.strand],str(abs(f.location.start-f.location.end)),'-',gene_name,synonym_code,"-",f.qualifiers["product"][0]])))
    fout.close()
    
    r = SeqIO.parse(annotation_file, "gb")
    fout = open(ptt_file, "w")
    for record in r:
       record.features = [f for f in record.features if f.type == "CDS"]
       fout.write("{0} - 0..{1}\n".format(record.description, len(record)))
       fout.write("{0} proteins\n".format(len(record.features)))
       fout.write("Location\tStrand\tLength\tPID\tGene\tSynonym Code\tCOG\tProduct\n")
       fout.write("{0}\n".format("\t".join(['0..0','+','0','-','EmptyGene','EmptyGene',"-",'EmptyGene necessary to avoid EDGE-pro bug'])))
       for f in record.features:
           gene_location = str(f.location.start)+".."+str(f.location.end)
           #Remove non-numeric characters from gene_location
           gene_location = ''.join([c for c in gene_location if c in '1234567890.']) 
           gene_name = f.qualifiers["locus_tag"][0]
           synonym_code = f.qualifiers["locus_tag"][0]
           fout.write("{0}\n".format("\t".join([gene_location,strand[f.strand],str(abs(f.location.start-f.location.end)),'-',gene_name,synonym_code,"-",f.qualifiers["product"][0]])))

    fout.close()

if __name__ == "__main__":
    main()
