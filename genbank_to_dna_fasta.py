#!/usr/bin/python3

import sys, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

def main():
# Read in genbank file containing list of genes and output a fasta with the locus tag and Rv# (if it exists) for each gene.

# Process Arguments (identify input and output locations)
    parser = argparse.ArgumentParser(description='Produces a fasta file containing DNA gene sequences from a Genbank file and the corresponding genome (fasta format).')
    parser.add_argument('-b','--genbank_file', type=str, help="Input file from Genbank (*.gbk, *.gb)")
    parser.add_argument('-o','--out_file', type=str, help="Output file (fasta format)")
    args = parser.parse_args()

    out_file = open(args.out_file, 'w')
    r = SeqIO.parse(args.genbank_file, "gb")
    for record in r:
        record.features = [f for f in record.features if f.type == "CDS"] 
        for f in record.features:
            gene_name = f.qualifiers["locus_tag"][0] + '|'
            gene_seq = f.extract(record.seq)
            prot_desc = f.qualifiers["product"][0]
            try:
                Rv_name = f.qualifiers["db_xref"][0]
                Rv_name = Rv_name[12:] + '|'
            except KeyError:
                Rv_name = 'NA' + '|'
            out_record = SeqRecord(gene_seq,id=gene_name + Rv_name, name="name", description=prot_desc)
            SeqIO.write(out_record, out_file, "fasta")
    out_file.close()

if __name__ == "__main__":
    main()
