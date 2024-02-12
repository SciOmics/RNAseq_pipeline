
caffeinate &

ROOT_FOLDER=/Volumes/8TB_SSD/Moce_RNA_ATAC/RNAseq
FASTQ_FOLDER=/Volumes/8TB_SSD/Moce_RNA_ATAC/RNAseq/fastq_files
REFERENCE_GENOME=/Volumes/8TB_SSD/Moce_RNA_ATAC/RNAseq/Ens_grch38_ht2/genome
SAMPLE_NAMES=(
    A172_control \
    A172_moce \
    BxPC3_control \
    BxPC3_moce \
    HCC1806_control \
    HCC1806_moce \
    HT144_control \
    HT144_moce \
    J82_control \
    J82_moce \
    MB2141_control \
    MB2141_moce \
    NCIH23_control \
    NCIH23_moce \
    PANC1_control \
    PANC1_moce \
    SN12C_control \
    SN12C_moce 
)
THREADS=52

#assess quality with fastqc
mkdir $ROOT_FOLDER/fastqc_results
fastqc -o ${ROOT_FOLDER}/fastqc_results -t ${THREADS} ${FASTQ_FOLDER}/*.fastq*


#remove adapter sequences and reads with <20 PHRED score
mkdir $ROOT_FOLDER/trimmed_files
mkdir $ROOT_FOLDER/trimmed_files/reports
for sample in $SAMPLE_NAMES
	do
	    trim_galore \
        --cores THREADS \
        --paired ${FASTQ_FOLDER}/${sample}_R1.fastq* \
        ${FASTQ_FOLDER}/${sample}_R2.fastq* \
        -o ${ROOT_FOLDER}/trimmed_files
	done
mv ${ROOT_FOLDER}/trimmed_files/*.txt ${ROOT_FOLDER}/trimmed_files/reports


#map reads
mkdir $ROOT_FOLDER/bam_files
for sample in $SAMPLE_NAMES
	do
        hisat2 \
            -p $THREADS \
            -x /Volumes/8TB_SSD/ref_genomes/Ens_grch38_ht2/genome \
            -1 /Volumes/8TB_SSD/Moce_RNA_ATAC/RNAseq/trimmed_files/${sample}_R1_val_1.fq.gz \
            -2 /Volumes/8TB_SSD/Moce_RNA_ATAC/RNAseq/trimmed_files/${sample}_R2_val_2.fq.gz \
            -S ${ROOT_FOLDER}/bam_files/${sample}.sam 
        samtools view -b \
            -@ $THREADS \
            ${ROOT_FOLDER}/bam_files/${sample}.sam > ${ROOT_FOLDER}/bam_files/${sample}.bam 
        rm ${ROOT_FOLDER}/bam_files/${sample}.sam
 	done


#sort and index reads
for sample in $SAMPLE_NAMES
    do
        samtools sort \
            -@ $THREADS \
            ${ROOT_FOLDER}/bam_files/${sample}.bam \
            -o ${ROOT_FOLDER}/bam_files/${sample}_sorted.bam
        samtools index \
            -@ $THREADS \
            ${ROOT_FOLDER}/bam_files/${sample}_sorted.bam
        rm ${ROOT_FOLDER}/bam_files/${sample}.bam
    done


#generate count table
featureCounts \
    -T $THREADS \
    -p \
    -s 0 \
    -a /Volumes/8TB_SSD/Moce_RNA_ATAC/RNAseq/Ens_grch38_ht2/Homo_sapiens.GRCh38.104.gtf \
    -o ${ROOT_FOLDER}/counts.out \
    ${ROOT_FOLDER}/bam_files/*_sorted.bam


#MultiQC Report
multiqc ${ROOT_FOLDER}

