# regulon-enrichment
As provided, this code replicates the analysis only from the output of EDGE-pro through the regulon enrichment test. For convenience, the EDGE-pro output file is provided with this code as deseqFile. You can also download the individual output EDGE-pro files from GEO: GSE104599 (*.rpkm). 

To replicate the regulon enrichment test of the RNA-seq and published literature datasets you will need to take the following steps. 

1. Modify the path names for the various applications used in this project. These can be changed by modifying the top few lines of the Makefile.

2. Track down the supplemental data for the prior publications (Voskuil et al 2004, and Betts et al 2002) and place them into the Downloads/ folder. See Makefile for expected filenames. In the Makefile I have provided a url that works on my machine to the Voskuil et al 2004 dataset but you will have to find the Betts et al 2002 supplement on your own (publisher did not use persistent links). 

3. Run ```make```. See makefile_dependency_graph.png and makefile_full_dependency_graph.png for a useful graph of the dependecy structure of the Makefile. This takes about 1 hour to run on my laptop (16 GB RAM, Intel i7-7700HQ).

4. Optional: For testing purposes you may wish to change the number of Monte-Carlo samples used for calculation with the regulon enrichment test. To do this edit num_MC_samp at the top of rnaseq_enrichment_analysis.py and lit_microarray_enrichment_analysis.py.

If, in addition, you want make to replicate the analysis done by EDGE-pro you will need to take the following steps:

1. Uncomment the following two lines from the Makefile: 
```bash
#deseqFile: $(SAMPLE_NAMES_WPATH)
#	${EDGE_PRO_FOLDER}/additionalScripts/edgeToDeseq.perl $(SAMPLE_NAMES_WPATH) 
```

2. Ensure you have EDGE-pro v1.3.1 installed. You will furthermore need to update EDGE_PRO_FOLDER at the top of the Makefile to be the location of these scripts.

3. Download the raw read files from GEO (GSE104599). Update READ_DATA_FOLDER to be the location for all the *.fastq (uncompressed) files.

4. Optional: If you wish to keep the large alignment files (sam format) generated by EDGE-pro out of the default folder you will need to update EDGE_PRO_OUTPUT_FOLDER. You could also make EDGE_PRO_OUTPUT_FOLDER a softlink to an external location.

If you have issues or feedback, please feel free to contact me:

Will Matern - maternwill@gmail.com
