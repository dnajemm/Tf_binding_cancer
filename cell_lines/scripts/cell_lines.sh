################################
# Load snakemake module
################################
module load snakemake/9.4.0


#######################################################################
# COAD cell lines RNA-seq and WGBS data download and analysis pipeline
#######################################################################
#########################################################################################
# HCT166 cell lines RNA-seq download and analysis pipeline
#########################################################################################

# Download the RNA-seq data for HCT116 cell lines from SRA using fasterq-dump done
for sample in RNA_HCT116_1 RNA_HCT116_2 ; do
    snakemake --cores 4 --latency-wait 180 --use-envmodules \
        -s /data/scripts/pipelines/1.0/getdata_pipeline.smk \
        -C sample=$sample geo=y quality=y \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD
done

# launch RNAseq pipeline for HCT116 cell lines done
for sample in RNA_HCT116_1 RNA_HCT116_2 ; do
    snakemake --cores 4 --latency-wait 150 --use-envmodules \
        -s /data/collotr/pipelines/RNAseq_pipeline.smk \
        -C sample=$sample track=y trim=y map=y multiqc=y all_genes=n counts=y \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD
done

# moved the rnaseq to /data/papers/yao_natCommun_2014 and do a symbolic link of the files from there to here 
ln -s /data/papers/yao_natCommun_2014/data/* /data/najemd/TF_binding_cancer/cell_lines/COAD/data
ln -s /data/papers/yao_natCommun_2014/reads/* /data/najemd/TF_binding_cancer/cell_lines/COAD/reads
ln -s /data/papers/yao_natCommun_2014/tracks/* /data/najemd/TF_binding_cancer/cell_lines/COAD/tracks
ln -s /data/papers/yao_natCommun_2014/expression/* /data/najemd/TF_binding_cancer/cell_lines/COAD/expression


#########################################################################################
# NCM460 cell line RNA-seq and RRBS data download and analysis pipeline
#########################################################################################

# Download the RNA-seq data and WGBS data for NCM460 and HCT116 cell lines from SRA using fasterq-dump done
for sample in RNA_NCM460_1 RNA_NCM460_2 RNA_NCM460_3 RRBS_NCM460_1 RRBS_HCT116_1; do
    snakemake --cores 4 --latency-wait 180 --use-envmodules \
        -s /data/scripts/pipelines/1.0/getdata_pipeline.smk \
        -C sample=$sample geo=y quality=y \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD
done

# launch RNAseq pipeline for NCM460 cell lines done
for sample in RNA_NCM460_1 RNA_NCM460_2 RNA_NCM460_3; do
    snakemake --cores 4 --latency-wait 150 --use-envmodules \
        -s /data/collotr/pipelines/RNAseq_pipeline.smk \
        -C sample=$sample track=y trim=y map=y multiqc=y all_genes=n counts=y \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD
done

# moved the rnaseq to /data/papers/shanS_natCommun_2026/ and do a symbolic link of the files from there to here 
ln -s /data/papers/shanS_natCommun_2026/data/* /data/najemd/TF_binding_cancer/cell_lines/COAD/data
ln -s /data/papers/shanS_natCommun_2026/reads/* /data/najemd/TF_binding_cancer/cell_lines/COAD/reads
ln -s /data/papers/shanS_natCommun_2026/tracks/* /data/najemd/TF_binding_cancer/cell_lines/COAD/tracks
ln -s /data/papers/shanS_natCommun_2026/expression/* /data/najemd/TF_binding_cancer/cell_lines/COAD/expression

# launch trimgalore pipeline for RRBS data for NCM460 and HCT116 cell lines done
for sample in RRBS_NCM460_1 RRBS_HCT116_1;do
    snakemake --cores 4 --latency-wait 180 --use-envmodules \
    -s /data/scripts/pipelines/1.0/trimgalore_pipeline.smk \
    -C sample=$sample multiqc=y \
    prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD
done

# launch bismark pipeline for RRBS data for NCM460 and HCT116 cell lines done
for sample in RRBS_NCM460_1 RRBS_HCT116_1; do
    snakemake --cores 4 --latency-wait 150 --use-envmodules \
    -s /data/scripts/pipelines/1.0/bismark_pipeline_2.smk \
    -C sample=$sample subset=n trimmed=y \
    prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD
done

# looking at the M bias in the RRBS data for NCM460 and HCT116 cell lines we have to clip the R2 reads to 5 bp in order to remove the bias in the methylation calls and coverage files for both cell lines.
# relaunch the bismark pipeline for RRBS data for NCM460 and HCT116 cell lines with the R2 reads clipped to 5 bp done
mkdir -p /data/najemd/TF_binding_cancer/cell_lines/COAD/logs

for sample in RRBS_NCM460_1 RRBS_HCT116_1; do
    while (( $(jobs -r | wc -l) >= 2 )); do
        wait -n
    done

    snakemake \
        --cores 4 \
        --latency-wait 180 \
        --use-envmodules \
        --directory /data/tmp/najemd/${sample} \
        -s /data/scripts/pipelines/1.0/trimgalore_pipeline.smk \
        -C sample="$sample" multiqc=y \
            args="-q 20 --stringency 2 --rrbs --non_directional --clip_R2 5" \
            prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD \
        > /data/najemd/TF_binding_cancer/cell_lines/COAD/logs/${sample}_trim.log 2>&1 &
done
wait
echo "Trimming done"

# ran the bismark pipeline done 
for sample in RRBS_NCM460_1 RRBS_HCT116_1; do
    while (( $(jobs -r | wc -l) >= 2 )); do
        wait -n
    done

    snakemake \
        --cores 4 \
        --retries 2 \
        --latency-wait 180 \
        --use-envmodules \
        --directory /data/tmp/najemd/${sample} \
        -s /data/collotr/pipelines/bismark_pipeline.smk \
        -C sample="$sample" subset=n trim=y \
            prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD \
        > /data/najemd/TF_binding_cancer/cell_lines/COAD/logs/${sample}_bismark.log 2>&1 &
done
wait
echo "All done"

# moved the RRBS to /data/papers/yu_BMCGenomics_2023/ and do a symbolic link of the files from there to here 
ln -s /data/papers/yu_BMCGenomics_2023/data/* /data/najemd/TF_binding_cancer/cell_lines/COAD/data
ln -s /data/papers/yu_BMCGenomics_2023/reads/* /data/najemd/TF_binding_cancer/cell_lines/COAD/reads
ln -s /data/papers/yu_BMCGenomics_2023/tracks/* /data/najemd/TF_binding_cancer/cell_lines/COAD/tracks
ln -s /data/papers/yu_BMCGenomics_2023/expression/* /data/najemd/TF_binding_cancer/cell_lines/COAD/expression
ln -s /data/papers/yu_BMCGenomics_2023/methylation/* /data/najemd/TF_binding_cancer/cell_lines/COAD/methylation

# Analysis 
# using the deseq2 script to analyze the RNA-seq data for HCT116 and NCM460 cell lines 
# step 1 create design file :
cat > /data/najemd/TF_binding_cancer/cell_lines/COAD/study_design.txt << 'EOF'
sample	condition	replicate
RNA_HCT116_1	HCT116	1
RNA_HCT116_2	HCT116	2
RNA_NCM460_1	NCM460	1
RNA_NCM460_2	NCM460	2
RNA_NCM460_3	NCM460	3
EOF

# launch the deseq2 script 
Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Deseq2.R \
    -n HCT116_vs_NCM460 \
    -d /data/najemd/TF_binding_cancer/cell_lines/COAD/study_design.txt \
    -c /data/najemd/TF_binding_cancer/cell_lines/COAD/expression \
    -o /data/najemd/TF_binding_cancer/cell_lines/COAD/results/expression/deseq2 \
    -t COAD

Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Expression_boxplots_TCGA_celllines.R COAD HCT116_tumor NCM460_healthy

#################################################################################################################################################################################

##########################################################################
# PRAD   cell lines RNA-seq and WGBS data download and analysis pipeline 
##########################################################################
mkdir -p /data/najemd/TF_binding_cancer/cell_lines/PRAD/

# Download the WGBS data for LNCaP and PrEC cell lines from SRA using fasterq-dump done 
for sample in WGBS_LNCaP_1 WGBS_LNCaP_2 WGBS_PreC_1 WGBS_PreC_2; do
    snakemake --cores 4 --latency-wait 180 --use-envmodules \
        -s /data/scripts/pipelines/1.0/getdata_pipeline.smk \
        -C sample=$sample geo=y quality=y \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/PRAD
done


# launch bismark pipeline for WGBS data for LNCaP and PreC cell lines HAVE NTO RELAUNCH IT with clipR2 three prime 1 and from 5' at 2  launched in WGBS_PRAD screen
mkdir -p /data/najemd/TF_binding_cancer/cell_lines/PRAD/logs

for sample in WGBS_LNCaP_1 WGBS_LNCaP_2 WGBS_PreC_1 WGBS_PreC_2; do
    while (( $(jobs -r | wc -l) >= 2 )); do
        wait -n
    done

    snakemake \
        --cores 4 \
        --retries 2 \
        --latency-wait 180 \
        --use-envmodules \
        --directory /data/tmp/najemd/${sample} \
        -s /data/collotr/pipelines/bismark_pipeline.smk \
        -C sample="$sample" subset=n trim=y trim_args="--clip_R2 2 --three_prime_clip_R2 1" \
            prepath=/data/najemd/TF_binding_cancer/cell_lines/PRAD \
        > /data/najemd/TF_binding_cancer/cell_lines/PRAD/logs/${sample}_bismark.log 2>&1 &
done
wait
echo "All done"


# Download the RNAseq data for LNCaP and PrEC cell lines from SRA using fasterq-dump done 
for sample in RNA_LNCaP_1 RNA_LNCaP_2 RNA_LNCaP_3 RNA_PreC_1 RNA_PreC_2 RNA_PreC_3; do
    snakemake --cores 4 --latency-wait 180 --use-envmodules \
        -s /data/scripts/pipelines/1.0/getdata_pipeline.smk \
        -C sample=$sample geo=y quality=y \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/PRAD
done

# launch RNAseq pipeline for PreC and LNCaP cell lines done
for sample in RNA_PreC_1 RNA_PreC_2 RNA_PreC_3 RNA_LNCaP_1 RNA_LNCaP_2 RNA_LNCaP_3 ; do
    snakemake --cores 4 --latency-wait 150 --use-envmodules \
        -s /data/collotr/pipelines/RNAseq_pipeline.smk \
        -C sample=$sample track=y trim=y map=y multiqc=y all_genes=n counts=y \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/PRAD
done

# moved the rnaseq to /data/papers/FrenchHJ_unpublished_2018/ and do a symbolic link of the files from there to here 
ln -s /data/papers/FrenchHJ_unpublished_2018/data/* /data/najemd/TF_binding_cancer/cell_lines/PRAD/data
ln -s /data/papers/FrenchHJ_unpublished_2018/reads/* /data/najemd/TF_binding_cancer/cell_lines/PRAD/reads
ln -s /data/papers/FrenchHJ_unpublished_2018/tracks/* /data/najemd/TF_binding_cancer/cell_lines/PRAD/tracks
ln -s /data/papers/FrenchHJ_unpublished_2018/expression/* /data/najemd/TF_binding_cancer/cell_lines/PRAD/expression


# Analysis 
# using the deseq2 script to analyze the RNA-seq data for PreC and LNCaP cell lines 
# step 1 create design file :
printf "sample\tcondition\treplicate\n\
RNA_PreC_1\tPreC\t1\n\
RNA_PreC_2\tPreC\t2\n\
RNA_PreC_3\tPreC\t3\n\
RNA_LNCaP_1\tLNCaP\t1\n\
RNA_LNCaP_2\tLNCaP\t2\n\
RNA_LNCaP_3\tLNCaP\t3\n" \
> /data/najemd/TF_binding_cancer/cell_lines/PRAD/study_design.txt

grep "PRAD" /data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs.tsv | cut -f3 | sort -u | wc -l
# 26

# launch the deseq2 script done
Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Deseq2.R \
    -n PreC_vs_LNCaP \
    -d /data/najemd/TF_binding_cancer/cell_lines/PRAD/study_design.txt \
    -c /data/najemd/TF_binding_cancer/cell_lines/PRAD/expression \
    -o /data/najemd/TF_binding_cancer/cell_lines/PRAD/results/expression/deseq2 \
    -t PRAD

Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Expression_boxplots_TCGA_celllines.R PRAD LNCaP_tumor PreC_healthy


#######################################################################
# BRCA cell lines RNA-seq and WGBS data download and analysis pipeline
#######################################################################
#########################################################################################
# HCC1954 and HMEC cell lines RNA-seq and WGBS data get from papers and analysis pipeline
#########################################################################################
  
mkdir -p /data/najemd/TF_binding_cancer/cell_lines/BRCA/data    

# get data and RNAseq done
run_sample() {
    local sample=$1
    local prepath=/data/najemd/TF_binding_cancer/cell_lines/BRCA
    local logdir=$prepath/logs

    echo "[$(date)] Starting getdata: $sample"
    snakemake --cores 4 --latency-wait 180 --use-envmodules \
        -s /data/scripts/pipelines/1.0/getdata_pipeline.smk \
        -C sample=$sample geo=y quality=y \
        prepath=$prepath \
        > $logdir/${sample}_getdata.log 2>&1

    if [[ $? -eq 0 ]]; then
        echo "[$(date)] getdata done for $sample — launching RNAseq"
        snakemake --cores 4 --latency-wait 150 --use-envmodules \
            -s /data/collotr/pipelines/RNAseq_pipeline.smk \
            -C sample=$sample track=y trim=y map=y multiqc=y all_genes=n counts=y \
            prepath=$prepath \
            > $logdir/${sample}_rnaseq.log 2>&1
        echo "[$(date)] RNAseq done for $sample"
    else
        echo "[$(date)] ERROR: getdata failed for $sample — skipping RNAseq" >&2
    fi
}

mkdir -p /data/najemd/TF_binding_cancer/cell_lines/BRCA/logs

for sample in RNA_HMEC_1 RNA_HMEC_2 RNA_HCC1954_1 RNA_HCC1954_2; do
    run_sample $sample &

    if (( $(jobs -r | wc -l) >= 2 )); then
        wait
    fi
done
wait
echo "All samples done"

# moved the rnaseq to /data/papers/FrenchHJ_unpublished_2018/ and do a symbolic link of the files from there to here 
ln -s /data/papers/hon_genomeRes_2012/data/RNA* /data/najemd/TF_binding_cancer/cell_lines/BRCA/data
ln -s /data/papers/hon_genomeRes_2012/reads/RNA* /data/najemd/TF_binding_cancer/cell_lines/BRCA/reads
ln -s /data/papers/hon_genomeRes_2012/tracks/RNA* /data/najemd/TF_binding_cancer/cell_lines/BRCA/tracks
ln -s /data/papers/hon_genomeRes_2012/expression/RNA* /data/najemd/TF_binding_cancer/cell_lines/BRCA/expression

# Analysis 
# using the deseq2 script to analyze the RNA-seq data for PreC and LNCaP cell lines 
# step 1 create design file :
printf "sample\tcondition\treplicate\n\
RNA_HMEC_1\tHMEC\t1\n\
RNA_HMEC_2\tHMEC\t2\n\
RNA_HCC1954_1\tHCC1954\t1\n\
RNA_HCC1954_2\tHCC1954\t2\n" \
> /data/najemd/TF_binding_cancer/cell_lines/BRCA/study_design.txt

grep "BRCA" /data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs.tsv | cut -f3 | sort -u | wc -l
# 30

# launch the deseq2 script done
Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Deseq2.R \
    -n HMEC_vs_HCC1954 \
    -d /data/najemd/TF_binding_cancer/cell_lines/BRCA/study_design.txt \
    -c /data/najemd/TF_binding_cancer/cell_lines/BRCA/expression \
    -o /data/najemd/TF_binding_cancer/cell_lines/BRCA/results/expression/deseq2 \
    -t BRCA


Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Expression_boxplots_TCGA_celllines.R BRCA HMEC_healthy HCC1954_tumor


# get data and launch bismark pipeline  done 
# Did it first and found it should be clipped at 2bp for r1 so redid it 
run_sample() {
    local sample=$1
    local prepath=/data/najemd/TF_binding_cancer/cell_lines/BRCA
    local logdir=$prepath/logs

    echo "[$(date)] Starting getdata: $sample"
    snakemake --cores 4 --latency-wait 180 --use-envmodules \
        -s /data/scripts/pipelines/1.0/getdata_pipeline.smk \
        -C sample=$sample geo=y quality=y \
            prepath=$prepath \
        > $logdir/${sample}_getdata.log 2>&1

    if [[ $? -eq 0 ]]; then
        echo "[$(date)] getdata done for $sample — launching Bismark"
        snakemake \
            --cores 4 \
            --retries 2 \
            --latency-wait 180 \
            --use-envmodules \
            --directory /data/tmp/najemd/${sample} \
            -s /data/collotr/pipelines/bismark_pipeline.smk \
            -C sample=$sample subset=n trim=y \
                trim_args="-q 20  --clip_R1 2" \
                prepath=$prepath \
            > $logdir/${sample}_bismark.log 2>&1
        echo "[$(date)] Bismark done for $sample"
    else
        echo "[$(date)] ERROR: getdata failed for $sample — skipping Bismark" >&2
    fi
}

mkdir -p /data/najemd/TF_binding_cancer/cell_lines/BRCA/logs

for sample in WGBS_HCC1954_1 WGBS_HCC1954_2 WGBS_HMEC_1 WGBS_HMEC_2; do
    run_sample $sample &

    if (( $(jobs -r | wc -l) >= 2 )); then
        wait
    fi
done
wait
echo "All BRCA WGBS samples done"

# plot the boxplots for the methylation data for HCC1954 and HMEC cell lines and TCGA BRCA samples
Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Methylation_boxplots_TCGA_celllines.R BRCA HCC1954_tumor HMEC_healthy

# plot the boxplots for the methylation data and expression data combined for HCC1954 and HMEC cell lines and TCGA BRCA samples
Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/combined_boxplot_TCGA_celllines.R BRCA HCC1954_tumor HMEC_healthy

#######################################################################
# LUAD cell lines RNA-seq and WGBS data download and analysis pipeline
#######################################################################