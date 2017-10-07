SHELL := /bin/bash

# Variables storing path and file names
EDGE_PRO_FOLDER := /home/will/src_and_bin/Other_Software/EDGE_pro_v1.3.1
BLAST_BIN := ../ncbi-blast-2.4.0+/bin
EDGE_PRO_OUTPUT_FOLDER := EDGE-pro_output
READ_DATA_FOLDER := RNAseq_rawReads

SAMPLE_NAMES := DR1 DR2 DR17 DR5 DR6 DR25 DR9 DR10 DR21 DR13 DR14 DR3 DR4 DR18 DR26 DR7 DR8 DR11 DR12 DR22 DR15 DR16
SAMPLE_NAMES_WPATH := $(addprefix $(EDGE_PRO_OUTPUT_FOLDER)/, $(join $(SAMPLE_NAMES), $(addprefix /, $(addsuffix .rpkm, $(SAMPLE_NAMES)))))

DOWNLOADS_PATH := Downloads
DOWNLOAD_FILES := Downloads/CDC1551_genomic.gbff Downloads/CDC1551_genome.fa
INTERMEDIATE_FILES := CDC1551.rnt CDC1551.ptt deseqFile
OUTPUT_FILES := Analysis_Output/7H9vshyp_low-read-rm.csv Analysis_Output/rnaseq_confusion_matrices_tf_hyp.xlsx Analysis_Output/lit_confusion_matrices_tf_hyp.xlsx
TUBERCULIST_FILES_WPATH := $(addprefix $(DOWNLOADS_PATH)/tuberculist_functional_category_20160722_, cell-wall-and-cell-processes.tsv conserved-hypothetical.tsv information-pathways.tsv insertion-seqs-and-phages.tsv intermediary-metabolism-and-respiration.tsv lipid-metabolism.tsv PE-PPE-families.tsv regulatory-proteins.tsv unknown.tsv virulence-detoxification-and-adaptation.tsv)

# Dummy targets
all: ${OUTPUT_FILES}

clean:
	#rm -f ${DOWNLOAD_FILES} ${INTERMEDIATE_FILES}
	@echo Note: alignment files in $(EDGE_PRO_OUTPUT_FOLDER) were not removed with clean

.PHONY: all clean

# Make directories
$(shell mkdir -p 2002_CDC1551_blast_db)

./Downloads/CDC1551_genomic.gbff: 
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Mycobacterium_tuberculosis/latest_assembly_versions/GCA_000669715.1_Myco_tube_cdc1551_V1/GCA_000669715.1_Myco_tube_cdc1551_V1_genomic.gbff.gz
	gunzip GCA_000669715.1_Myco_tube_cdc1551_V1_genomic.gbff.gz
	mv GCA_000669715.1_Myco_tube_cdc1551_V1_genomic.gbff Downloads/CDC1551_genomic.gbff
	rm -f GCA_000669715.1_Myco_tube_cdc1551_V1_genomic.gbff.gz

./Downloads/CDC1551_genome.fa:
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Mycobacterium_tuberculosis/latest_assembly_versions/GCA_000669715.1_Myco_tube_cdc1551_V1/GCA_000669715.1_Myco_tube_cdc1551_V1_genomic.fna.gz
	gunzip GCA_000669715.1_Myco_tube_cdc1551_V1_genomic.fna.gz
	mv GCA_000669715.1_Myco_tube_cdc1551_V1_genomic.fna Downloads/CDC1551_genome.fa
	rm -f GCA_000669715.1_Myco_tube_cdc1551_V1_genomic.fna.gz

# Download CDC1551 sequenced in 2002 genome and annotations (genbank).
./Downloads/2002_CDC1551_annot.gbk:
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Bacteria/Mycobacterium_tuberculosis_CDC1551_uid223/AE000516.gbk
	mv AE000516.gbk ./Downloads/2002_CDC1551_annot.gbk

# Make fasta files from genbank files for 2002_CDC1551 and CDC1551_genomic.gbff:
CDC1551_genes.fa: Downloads/CDC1551_genomic.gbff
	./genbank_to_dna_fasta.py -b Downloads/CDC1551_genomic.gbff -o CDC1551_genes.fa

2002_CDC1551_blast_db/2002_CDC1551_genelist.fa: Downloads/2002_CDC1551_annot.gbk
	./genbank_to_dna_fasta.py -b Downloads/2002_CDC1551_annot.gbk -o 2002_CDC1551_blast_db/2002_CDC1551_genelist.fa

# Making blast db for 2002_CDC1551 genome: 
2002_CDC1551_blast_db/2002_CDC1551.nsq: 2002_CDC1551_blast_db/2002_CDC1551_genelist.fa
	$(BLAST_BIN)/makeblastdb -in 2002_CDC1551_blast_db/2002_CDC1551_genelist.fa -parse_seqids -dbtype nucl -out 2002_CDC1551_blast_db/2002_CDC1551

#Find all homologs of 2014_CDC1551 proteins in 2002_CDC1551 proteome:
blast_results.xml: 2002_CDC1551_blast_db/2002_CDC1551.nsq CDC1551_genes.fa
	$(BLAST_BIN)/blastn -out blast_results.xml -task dc-megablast -outfmt 5 -query CDC1551_genes.fa -db 2002_CDC1551_blast_db/2002_CDC1551 -evalue .001 -max_hsps 1 -max_target_seqs 50

#Dictionary for mapping between 2014_CDC1551 genes and Rv and MT.
V735_to_Rv_MT_annot_mapping.csv: blast_results.xml CDC1551_genes.fa
	./make_dict_V735_Rv_MT_annot.py

# Process CDC1551_genomic.gbff into *.rnt and *.ptt for input into EDGE-pro
CDC1551.rnt CDC1551.ptt: Downloads/CDC1551_genomic.gbff
	./create_rnt_ptt.py Downloads/CDC1551_genomic.gbff -r CDC1551.rnt -p CDC1551.ptt

# Run EDGE-pro on multiple samples (specified in SAMPLE_NAMES)
$(SAMPLE_NAMES_WPATH): %: CDC1551.rnt CDC1551.ptt Downloads/CDC1551_genome.fa $(READ_DATA_FOLDER)
	$(eval D := $(notdir $(@:.rpkm=))) #Removes the path and suffix
	mkdir -p ${EDGE_PRO_OUTPUT_FOLDER}/$D
	$(eval READ1 = $(shell find $(READ_DATA_FOLDER)/*$D_R1*))
	$(eval READ2 = $(shell find $(READ_DATA_FOLDER)/*$D_R2*))
	${EDGE_PRO_FOLDER}/edge.pl -u ${READ1} -v ${READ2} -t 2 -o ${EDGE_PRO_OUTPUT_FOLDER}/$D/$D -g Downloads/CDC1551_genome.fa -r CDC1551.rnt -p CDC1551.ptt -s ${EDGE_PRO_FOLDER} > ${EDGE_PRO_OUTPUT_FOLDER}/$D/EDGE-pro_mapping_output
	cp ${EDGE_PRO_OUTPUT_FOLDER}/$D/$D.rpkm_0 ${EDGE_PRO_OUTPUT_FOLDER}/$D/$D.rpkm #These lines combine the two .rpkm files into one.
	tail -n +3 ${EDGE_PRO_OUTPUT_FOLDER}/$D/$D.rpkm_1 >> ${EDGE_PRO_OUTPUT_FOLDER}/$D/$D.rpkm

# Combine all EDGE-pro files into a single input file for DE-seq
#deseqFile: $(SAMPLE_NAMES_WPATH)
#	${EDGE_PRO_FOLDER}/additionalScripts/edgeToDeseq.perl $(SAMPLE_NAMES_WPATH) 

# Convert downloaded Tuberculist files into mapping between Rv# and functional category. 
# Note: Due to the setup of Tuberculist these files must unfortunately be downloaded by hand and named correctly.
all_functional_categories.csv: $(TUBERCULIST_FILES_WPATH)
	./make_all_functional_categories.py

# Use DE-seq to look for significantly up/downregulated genes
Analysis_Output/7H9vshyp_low-read-rm.csv: deseqFile V735_to_Rv_MT_annot_mapping.csv all_functional_categories.csv run_DEseq_analysis.R DEseq_analysis_functions.R
	mkdir -p Analysis_Output
	./run_DEseq_analysis.R

Downloads/tfoe.searchable_130115.xlsx:
	wget http://networks.systemsbiology.net/mtb/sites/default/files/file_attach/tfoe.searchable_130115.xlsx
	mv tfoe.searchable_130115.xlsx Downloads/tfoe.searchable_130115.xlsx

#Do Enrichment Analysis using RNA-seq data
Analysis_Output/rnaseq_confusion_matrices_tf_hyp.xlsx: Analysis_Output/7H9vshyp_low-read-rm.csv Downloads/tfoe.searchable_130115.xlsx
	./rnaseq_enrichment_analysis.py

#Do enrichment analysis using previously published microarray data
Downloads/1-s2.0-S147297920400023X-mmc1.xls:
	############# This link worked as of 2017-10-06 but the publisher will probably move it so if this fails you will need to track down 1-s2.0-S147297920400023X-mmc1.xls from the supplement of Voskuil et al 2004 and place it into the Downloads/ folder.
	wget https://ars.els-cdn.com/content/image/1-s2.0-S147297920400023X-mmc1.xls
	mv 1-s2.0-S147297920400023X-mmc1.xls Downloads/1-s2.0-S147297920400023X-mmc1.xls

Downloads/MMI_2779_sm_sup.xlsx:
	################ You will need to track down MMI_2779_sm_sup.xlsx from the supplement of Betts et al 2002 and place it into the Downloads/ folder'

Analysis_Output/lit_confusion_matrices_tf_hyp.xlsx: Analysis_Output/7H9vshyp_low-read-rm.csv Downloads/tfoe.searchable_130115.xlsx Downloads/1-s2.0-S147297920400023X-mmc1.xls Downloads/MMI_2779_sm_sup.xlsx
	./lit_microarray_enrichment_analysis.py

