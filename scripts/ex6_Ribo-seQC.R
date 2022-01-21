#Install and load packages
install.packages("devtools", repos='http://cran.us.r-project.org')
install.packages("rmarkdown", repos='http://cran.us.r-project.org')
install.packages("knitr", repos='http://cran.us.r-project.org')
library("rmarkdown")
library("knitr")
library("devtools")
install_github(repo = "lcalviell/Ribo-seQC")
library("RiboseQC")

#Set working directory
setwd("<insert here your working directory>")

#Prepare genome file (done once)
prepare_annotation_files(annotation_directory = "./",
                         twobit_file = "Rattus_norvegicus.Rnor_6.0.dna.toplevel.2bit",
                         gtf_file = "Rattus_norvegicus.Rnor_6.0.104.gtf",
                         scientific_name = "Rattus.norvegicus",
                         annotation_name = "Rnor_6",
                         export_bed_tables_TxDb = F,
                         forge_BSgenome = T,
                         create_TxDb = T)

#Get mapped genomes without undesired RNAs as bam files
genome_bam <- c("Biever_Neuropil_Poly_1_clpd_tr.fastq_no_r_sno_sn_t_RNA_Rnor_6_sorted.bam",
                "Biever_Neuropil_Poly_2_clpd_tr.fastq_no_r_sno_sn_t_RNA_Rnor_6_sorted.bam",
                "Biever_Neuropil_Poly_3_clpd_tr.fastq_no_r_sno_sn_t_RNA_Rnor_6_sorted.bam",
                "Biever_Somata_Poly_1_clpd_tr.fastq_no_r_sno_sn_t_RNA_Rnor_6_sorted.bam",
                "Biever_Somata_Poly_2_clpd_tr.fastq_no_r_sno_sn_t_RNA_Rnor_6_sorted.bam",
                "Biever_Somata_Poly_3_clpd_tr.fastq_no_r_sno_sn_t_RNA_Rnor_6_sorted.bam")

#Generate quality control plots
RiboseQC_analysis(annotation_file ="Rattus_norvegicus.Rnor_6.0.104.gtf_Rannot",
                  bam_files = genome_bam,
                  fast_mode = T,
                  report_file = "Rnor_Biever_QC.html",
                  sample_names = c("Neuropil_Poly_1", "Neuropil_Poly_2", "Neuropil_Poly_3",
                                   "Somata_Poly_1", "Somata_Poly_2", "Somata_Poly_3"),
                  dest_names = c("Neuropil_Poly_1", "Neuropil_Poly_2", "Neuropil_Poly_3",
                                 "Somata_Poly_1", "Somata_Poly_2", "Somata_Poly_3"),
                  write_tmp_files = F)