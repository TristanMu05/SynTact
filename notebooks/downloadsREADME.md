You need to download the following datasets and place them in the data/raw/ folder before running the pipeline:

1. gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz
Source: GENCODE v48
Description: Comprehensive human gene annotation (reference chromosomes, scaffolds, haplotypes)
Download page:
🔗 https://www.gencodegenes.org/human/release_48.html

Direct download link:
🔗 https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.chr_patch_hapl_scaff.annotation.gtf.gz

2. hpa_subcellular_location.tsv
Source: Human Protein Atlas (HPA)
Description: Subcellular location annotations (used to detect membrane-bound proteins)
Download page:
🔗 https://www.proteinatlas.org/humanproteome/subcellular/data#locations

Instructions:

Download the .zip file from the site

Extract it

Copy only the file: hpa_subcellular_location.tsv into data/raw/

3. tcga_ov_expression.tsv.gz
Source: UCSC Xena / TCGA OV
Description: Gene expression RNA-seq (TPM) for ovarian cancer from TCGA
Dataset: STAR - TPM (n=429) GDC Hub
Docs: https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Ovarian%20Cancer%20(OV)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

Direct dataset page:
🔗 https://xenabrowser.net/datapages/?dataset=TCGA-OV.star_tpm.tsv&host=https%3A%2F%2Fgdc.xenahubs.net

Instructions:

Download the .tsv.gz file from the page

Rename it (if needed) to tcga_ov_expression.tsv.gz

Move it to data/raw/

4. GTEx_gene_tpm.gct.gz
Source: GTEx v10
Description: Gene expression (TPM) across normal human tissues
Download page:
🔗 https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression

Dataset name:
GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz

Gene Expression Transcripts Per Million (TPM) from RNASeQC v2.4.2

Direct download link:
🔗 https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_tpm.gct.gz

Instructions:

Rename the file to: GTEx_gene_tpm_v10.gct.gz

Move it to data/raw/