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

# dowload the HCT116 WGBS 
# get data and launch bismark pipeline  done
# found that R2 should be clipped at 3 

sample=WGBS_HCT116_1
mkdir -p /data/najemd/TF_binding_cancer/cell_lines/COAD/logs

snakemake --cores 4 --latency-wait 180 --use-envmodules \
    -s /data/scripts/pipelines/1.0/getdata_pipeline.smk \
    -C sample=${sample} geo=y quality=y \
    prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD \
    > /data/najemd/TF_binding_cancer/cell_lines/COAD/logs/${sample}_getdata.log 2>&1

sample=WGBS_HCT116_1 # done
snakemake --cores 4 --latency-wait 180 --use-envmodules \
    -s /data/collotr/pipelines/bismark_pipeline.smk \
    -C sample="$sample" subset=n trim=y \
        trim_args="--clip_R2 3" \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/COAD \
    > /data/najemd/TF_binding_cancer/cell_lines/COAD/logs/${sample}_bismark.log 2>&1 

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

Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Expression_boxplots_TCGA_celllines.R COAD NCM460_healthy HCT116_tumor
Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Methylation_boxplots_TCGA_celllines.R COAD NCM460_healthy HCT116_tumor
Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/combined_boxplot_TCGA_celllines.R COAD NCM460_healthy HCT116_tumor

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


# launch bismark pipeline for WGBS data for LNCaP and PreC cell lines HAVE NTO RELAUNCH IT with clipR2 three prime 1 and from 5' at 2   done
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


# plot the boxplots for the methylation data for HCC1954 and HMEC cell lines and TCGA PRAD samples
Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/Methylation_boxplots_TCGA_celllines.R PRAD PreC_healthy LNCaP_tumor

# plot the boxplots for the methylation data and expression data combined for HCC1954 and HMEC cell lines and TCGA PRAD samples
Rscript /data/najemd/TF_binding_cancer/cell_lines/scripts/combined_boxplot_TCGA_celllines.R PRAD PreC_healthy LNCaP_tumor


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


# get data and launch bismark pipeline  done done FOR THE _3 
# Did it first and found it should be clipped at 2bp for r1 so redid it 
export _JAVA_OPTIONS="-Djava.io.tmpdir=/data/tmp/najemd/"
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

for sample in  WGBS_HCC1954_2 WGBS_HMEC_1 WGBS_HMEC_2 WGBS_HCC1954_3 ; do
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

#######################
# ChIP HCC1954
#######################
# BANP
# Download the ChIP data for HCC1954 from SRA using fasterq-dump done 
export _JAVA_OPTIONS="-Djava.io.tmpdir=/data/tmp/najemd/"
for sample in BANP_HCC1954_1 BANP_HCC1954_2; do
    snakemake --cores 4 --latency-wait 180 --use-envmodules \
        -s /data/scripts/pipelines/1.0/getdata_pipeline.smk \
        -C sample=$sample geo=y quality=y \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/BRCA \
        2>&1 | tee /data/najemd/TF_binding_cancer/cell_lines/BRCA/logs/${sample}_getdata.log
done

# launch the chipseq pipeline done
export _JAVA_OPTIONS="-Djava.io.tmpdir=/data/tmp/najemd/"
PROJECT="/data/najemd/TF_binding_cancer/cell_lines/BRCA"
MAX_PARALLEL=2

cd "$PROJECT"

samples=$(awk '!/^#/ && NF > 0 {print $1}' "$PROJECT/SAMPLEFILE" | grep -E '^(BANP)_')

running=0

for sample in $samples; do

    echo "===================================="
    echo "Processing: $sample"
    echo "===================================="

    snakemake -s /data/bancquaa/pipelines/chipseq_pipeline.smk \
        -C sample="$sample" \
           prepath="$PROJECT" \
           trim=y \
           map=y \
           track=y \
        --cores 2 &

    running=$((running + 1))

    # Wait when 2 samples are running
    if [ "$running" -ge "$MAX_PARALLEL" ]; then
        wait
        running=0
    fi

done

# Wait for remaining samples
wait

echo "All samples processed successfully."

# NRF1
# Download the ChIP data for HCC1954 from SRA using fasterq-dump done 
export _JAVA_OPTIONS="-Djava.io.tmpdir=/data/tmp/najemd/"
for sample in NRF1_HCC1954_1 NRF1_HCC1954_2 inputNRF1_HCC1954_1; do
    snakemake --cores 4 --latency-wait 180 --use-envmodules \
        -s /data/scripts/pipelines/1.0/getdata_pipeline.smk \
        -C sample=$sample geo=y quality=y \
        prepath=/data/najemd/TF_binding_cancer/cell_lines/BRCA \
        2>&1 | tee /data/najemd/TF_binding_cancer/cell_lines/BRCA/logs/${sample}_getdata.log
done

# launch the chipseq pipeline done
export _JAVA_OPTIONS="-Djava.io.tmpdir=/data/tmp/najemd/"
PROJECT="/data/najemd/TF_binding_cancer/cell_lines/BRCA"
MAX_PARALLEL=3

cd "$PROJECT"

samples=$(awk '!/^#/ && NF > 0 {print $1}' "$PROJECT/SAMPLEFILE" | grep -E '^(NRF1|inputNRF1)_')

running=0

for sample in $samples; do

    echo "===================================="
    echo "Processing: $sample"
    echo "===================================="

    snakemake -s /data/bancquaa/pipelines/chipseq_pipeline.smk \
        -C sample="$sample" \
           prepath="$PROJECT" \
           trim=y \
           map=y \
           track=y \
        --cores 2 &

    running=$((running + 1))

    # Wait when 3 samples are running
    if [ "$running" -ge "$MAX_PARALLEL" ]; then
        wait
        running=0
    fi

done

# Wait for remaining samples
wait

echo "All samples processed successfully."

# launch the macs2 pipeline done
export PATH="/data/bardeta/scripts:$PATH"
PROJECT="/data/najemd/TF_binding_cancer/cell_lines/BRCA"
MACS_PIPELINE="/data/bancquaa/pipelines/macs_pipeline.smk"
MAX_PARALLEL=4

cd "$PROJECT"

samples=$(awk '!/^#/ && NF > 0 {print $1}' "$PROJECT/SAMPLEFILE" | grep -E '^(BANP|NRF1)_')
running=0

for sample in $samples; do
    echo "===================================="
    echo "MACS2 peak calling: $sample"
    echo "===================================="

    snakemake -s "$MACS_PIPELINE" \
        -C sample="$sample" \
           prepath="$PROJECT" \
           type=narrow \
        --cores 2 &

    running=$((running + 1))

    if [ "$running" -ge "$MAX_PARALLEL" ]; then
        wait
        running=0
    fi
done

wait
echo "All samples peak-called successfully."

###############################
# Create a Design File for each TF
###############################

awk -F'\t' '
  /^#/ { next }
  {
    name = $1
    rep = name; sub(/.*_/, "", rep)
    study = name; sub(/_.*/, "", study)
    cond = name; sub(/^[^_]*_/, "", cond); sub(/_[^_]*$/, "", cond)
    print name "\t" cond "\t" rep >> "/data/najemd/TF_binding_cancer/cell_lines/BRCA/" study "_sample_design.txt"
  }
' <(grep -E '^(NRF1|BANP)_' /data/najemd/TF_binding_cancer/cell_lines/BRCA/SAMPLEFILE)

for f in /data/najemd/TF_binding_cancer/cell_lines/BRCA/*_sample_design.txt; do
    sed -i '1s/^/sample\tcondition\treplicate\n/' "$f"
done

###############################
# Peak scatter done
###############################

for design in $(grep -E '^(BANP|NRF1)_' /data/najemd/TF_binding_cancer/cell_lines/BRCA/SAMPLEFILE | cut -f1 | cut -d'_' -f1 | sort -u); do
    /data/bardeta/scripts/peak_scatter_all_pipeline.sh \
        -n "$design" \
        -d /data/najemd/TF_binding_cancer/cell_lines/BRCA/${design}_sample_design.txt \
        -m macs
done

###############################
# Peak Count 
###############################

for design in $(grep -E '^(BANP|NRF1)_' /data/najemd/TF_binding_cancer/cell_lines/BRCA/SAMPLEFILE | cut -f1 | cut -d'_' -f1 | sort -u); do
    /data/najemd/Delphine_Cancer/peak_count_pipeline.sh \
        -n "$design" \
        -d /data/najemd/TF_binding_cancer/cell_lines/BRCA/${design}_sample_design.txt \
        -m macs
done

###########################
# Look if y significant motifs of BRCA fall within the peaks in HCC1954 for BANP and NRF1 done IM HERE SHOULD ADD THE INFO TO THE METH FILE 
###########################
# merge the peaks for the replicates for each TF
cat /data/najemd/TF_binding_cancer/cell_lines/BRCA/peaks/BANP_HCC1954_1_peaks_macs.bed \
    /data/najemd/TF_binding_cancer/cell_lines/BRCA/peaks/BANP_HCC1954_2_peaks_macs.bed \
| sort -k1,1 -k2,2n \
| bedtools merge > /data/najemd/TF_binding_cancer/cell_lines/BRCA/peaks/BANP_merged_peaks.bed
cat /data/najemd/TF_binding_cancer/cell_lines/BRCA/peaks/NRF1_HCC1954_1_peaks_macs.bed \
    /data/najemd/TF_binding_cancer/cell_lines/BRCA/peaks/NRF1_HCC1954_2_peaks_macs.bed \
| sort -k1,1 -k2,2n \
| bedtools merge > /data/najemd/TF_binding_cancer/cell_lines/BRCA/peaks/NRF1_merged_peaks.bed

# overlap my consistent motifsi n BRCA with 
for tf in BANP NRF1; do
    awk -F'\t' 'NR>1 && $4=="BRCA" && $24=="'$tf'"' \
        /data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_cpgi.tsv \
    | awk -F'\t' '{split($2, a, ":"); print a[1]"\t"a[2]"\t"a[3]"\t"$2"\t"$1"\t"a[4]}' \
    | bedtools intersect \
        -a stdin \
        -b /data/najemd/TF_binding_cancer/cell_lines/BRCA/peaks/${tf}_merged_peaks.bed \
        -wa -u \
    > /data/najemd/TF_binding_cancer/results/integration/motifs_${tf}_HCC1954_in_peaks.bed

    total=$(awk -F'\t' 'NR>1 && $4=="BRCA" && $24=="'$tf'"' \
        /data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_cpgi.tsv | wc -l)
    overlap=$(wc -l < /data/najemd/TF_binding_cancer/results/integration/motifs_${tf}_HCC1954_in_peaks.bed)
    echo "$tf: $overlap / $total BRCA motif instances overlap with peaks"
done

# BANP: 0 / 6 BRCA motif instances overlap with peaks
# NRF1: 3 / 30 BRCA motif instances overlap with peaks

################################################################################################################################################
# Adding the HCC1954 and HMEC cell lines information to the meth file if the gene was concordant and if it falls in HCC1954 peaks of the TF and the orientation of the genes if upstream or downstream 
################################################################################################################################################
Rscript -e '
library(data.table)

main <- fread("/data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_cpgi.tsv")
conc <- fread("/data/najemd/TF_binding_cancer/cell_lines/BRCA/results/BRCA_concordant_genes_HCC1954_HMEC.tsv")
tss <- fread("/data/genome/annotations/hg38_tss.bed",
             col.names = c("tss_chr","tss_start","tss_end","gene_name","ensembl_gene","ensembl_tx","strand"))

# Keep one strand per gene
tss_strand <- unique(tss[, .(gene_name, strand)])

# Join strand
main <- merge(main, tss_strand, by.x = "gene", by.y = "gene_name", all.x = TRUE)

# Concordance
main[, concordant_HCC1954_HMEC := fifelse(
  cancer == "BRCA",
  paste(probe, gene) %in% paste(conc$probe, conc$gene),
  NA
)]

# Orientation
main[, tss_orientation := fcase(
  strand == "+" & dist_to_tss > 0, "downstream",
  strand == "+" & dist_to_tss < 0, "upstream",
  strand == "-" & dist_to_tss > 0, "upstream",
  strand == "-" & dist_to_tss < 0, "downstream",
  dist_to_tss == 0, "at_TSS",
  default = NA_character_
)]

main[, strand := NULL]

fwrite(main,
       "/data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_cpgi.tsv",
       sep = "\t")
cat("Upstream:", sum(main$tss_orientation == "upstream", na.rm=TRUE), "\n")
cat("Downstream:", sum(main$tss_orientation == "downstream", na.rm=TRUE), "\n")
'

# adding the HCC1954 peaks info 
Rscript -e '
library(data.table)

main <- fread("/data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_cpgi.tsv")

nrf1 <- fread("/data/najemd/TF_binding_cancer/results/integration/motifs_NRF1_HCC1954_in_peaks.bed",
              header = FALSE)
banp <- fread("/data/najemd/TF_binding_cancer/results/integration/motifs_BANP_HCC1954_in_peaks.bed",
              header = FALSE)

main[, in_NRF1_HCC1954_peaks := fifelse(
  cancer == "BRCA" & TF == "NRF1",
  motif_instance_id %in% nrf1$V4,
  NA
)]

main[, in_BANP_HCC1954_peaks := fifelse(
  cancer == "BRCA" & TF == "BANP",
  motif_instance_id %in% banp$V4,
  NA
)]

fwrite(main,
       "/data/najemd/TF_binding_cancer/results/integration/meth_expr_consistent_pairs_cancermine_cpgi.tsv",
       sep = "\t")
cat("NRF1 BRCA in peaks:", sum(main$in_NRF1_HCC1954_peaks, na.rm=TRUE), "\n")
cat("BANP BRCA in peaks:", sum(main$in_BANP_HCC1954_peaks, na.rm=TRUE), "\n")
'



###########################################
# MCF-7 and MCF10A cell line launched in MCF screen for MCF10A withoutsubset and for MCF10A  in screen MCF10A
###########################################
# found that MCF-7 should have R2 clipped at 7 while RRBS doesnt need clipping
export _JAVA_OPTIONS="-Djava.io.tmpdir=/data/tmp/najemd/"
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
            -s /data/scripts/pipelines/1.0/bismark_pipeline.smk \
            -C sample=$sample subset=n trim=y delete=y \
                prepath=$prepath \
            > $logdir/${sample}_bismark.log 2>&1
        echo "[$(date)] Bismark done for $sample"
    else
        echo "[$(date)] ERROR: getdata failed for $sample — skipping Bismark" >&2
    fi
}

mkdir -p /data/najemd/TF_binding_cancer/cell_lines/BRCA/logs

for sample in RRBS_MCF10A ; do
    run_sample $sample &

    if (( $(jobs -r | wc -l) >= 2 )); then
        wait
    fi
done
wait
echo "All BRCA WGBS samples done"







#######################################################################
# LUAD cell lines RNA-seq and WGBS data download and analysis pipeline NOT LAUNCHED YET
#######################################################################
# get data and launch bismark pipeline   
run_sample() {
    local sample=$1
    local prepath=/data/najemd/TF_binding_cancer/cell_lines/LUAD
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
            -C sample=$sample subset=y trim=y \
                prepath=$prepath \
            > $logdir/${sample}_bismark.log 2>&1
        echo "[$(date)] Bismark done for $sample"
    else
        echo "[$(date)] ERROR: getdata failed for $sample — skipping Bismark" >&2
    fi
}

mkdir -p /data/najemd/TF_binding_cancer/cell_lines/LUAD/logs

for sample in WGBS_A549_1 WGBS_A549_2; do
    run_sample $sample &

    if (( $(jobs -r | wc -l) >= 2 )); then
        wait
    fi
done
wait
echo "All LUAD WGBS samples done"