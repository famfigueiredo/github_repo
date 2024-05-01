# Re-running Quantseq data from 2018/2020
## Alex's 2018 runs


# 15/04/2024


# Pool1_TKD180900288_HNJTWCCXY_L7_1.fq.gz = 03_lane3_read1.fq.gz
### Lane 3
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane3_alex \
-b /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane3_alex/demux-summary/lane3.lostreads.fq.gz \
-s /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane3_alex/demux-summary/lane3.summary.txt \
/data/ffi007/01_quantseq/03_data/05_useful-files/lane3-barcodes.txt \
/data/ffi007/01_quantseq/03_data/01_source-files/03_lane3-read1.fq.gz


fastqc -o /data/ffi007/01_quantseq/03_data/02_demuxed-files/  -t 72 -f fastq *.fq.gz
multiqc .


#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane3.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane3_alex/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane3_alex/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane3-bbduklog.txt

done

fastqc -o /data/ffi007/01_quantseq/03_data/03_clean-files/lane3_alex/fastqc -t 72 -f fastq *.fq.gz


# 16/04/2024


#!/bin/bash

ulimit -n 10000
ulimit -n

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane3.txt |
while read line
do
STAR --runThreadN 70 --genomeDir /data/ffi007/01_quantseq/02_star_genome --readFilesCommand zcat \
 --readFilesIn /data/ffi007/01_quantseq/03_data/03_clean-files/lane3_alex/${line}.clean.fq.gz \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/ffi007/01_quantseq/03_data/07_STAR-out/heart_data_ensembl/lane3_alex/${line}_ |& tee -a starlog-lane3-alex.txt \

done


rename 's/_Aligned.sortedByCoord.out/_sorted/' *.bam

#!/bin/bash

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane3.txt | \
parallel -j 57 "htseq-count -m intersection-nonempty -i gene_id -t exon -s yes -f bam -r pos \
/data/ffi007/01_quantseq/03_data/07_STAR-out/heart_data_ensembl/lane3_alex/{}_sorted.bam \
/data/ffi007/01_quantseq/03_data/05_useful-files/genomes/Salmo_salar.Ssal_v3.1.111.gtf \
>/data/ffi007/01_quantseq/03_data/08_ht-seq-count/lane3_alex/{}_readcount.txt" |& tee -a htseqlog-lane3-alex.txt



# 17/04/2024

## Demultiplexing

### Pool2_TKD180900289_HNJWCCCXY_L7_1.fq.gz = 04_lane4_read1.fq.gz ###

demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane4_alex \
-b /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane4_alex/demux-summary/lane4.lostreads.fq.gz \
-s /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane4_alex/demux-summary/lane4.summary.txt \
/data/ffi007/01_quantseq/03_data/05_useful-files/lane4-barcodes.txt \
/data/ffi007/01_quantseq/03_data/01_source-files/04_lane4-read1.fq.gz


fastqc -o /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane4_alex/  -t 72 -f fastq *.fq.gz
multiqc .


# 18/04/2024


## Trimming

#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane4.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane4_alex/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane4_alex/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane4-bbduklog.txt

done

fastqc -o /data/ffi007/01_quantseq/03_data/03_clean-files/lane4_alex/fastqc -t 72 -f fastq *.fq.gz

# 19/04/2024


## STAR alignment

#!/bin/bash

ulimit -n 10000
ulimit -n

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane4.txt |
while read line
do
STAR --runThreadN 70 --genomeDir /data/ffi007/01_quantseq/02_star_genome --readFilesCommand zcat \
 --readFilesIn /data/ffi007/01_quantseq/03_data/03_clean-files/lane4_alex/${line}.clean.fq.gz \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/ffi007/01_quantseq/03_data/07_STAR-out/heart_data_ensembl/lane4_alex/${line}_ |& tee -a starlog-lane4-alex.txt \

done


rename 's/_Aligned.sortedByCoord.out/_sorted/' *.bam


## HTseq-count

#!/bin/bash

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane4.txt | \
parallel -j 57 "htseq-count -m intersection-nonempty -i gene_id -t exon -s yes -f bam -r pos \
/data/ffi007/01_quantseq/03_data/07_STAR-out/heart_data_ensembl/lane3_alex/{}_sorted.bam \
/data/ffi007/01_quantseq/03_data/05_useful-files/genomes/Salmo_salar.Ssal_v3.1.111.gtf \
>/data/ffi007/01_quantseq/03_data/08_ht-seq-count/lane4_alex/{}_readcount.txt" |& tee -a htseqlog-lane4-alex.txt



# 30/04/2024


### Pool1_TKD180900288_HNJWCCCXY_L8_1.fq.gz = 05_lane5-read1.fq.gz

## Demultiplexing

demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5_alex \
-b /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5_alex/demux-summary/lane5.lostreads.fq.gz \
-s /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5_alex/demux-summary/lane5.summary.txt \
/data/ffi007/01_quantseq/03_data/05_useful-files/lane5-barcodes.txt \
/data/ffi007/01_quantseq/03_data/01_source-files/05_lane5-read1.fq.gz


fastqc -o /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5_alex/  -t 72 -f fastq *.fq.gz
multiqc .


# 01/05/2024

## Trimming

#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane5.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5_alex/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5_alex/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane5-bbduklog.txt

done

## Repeating the trimming process since I got an error on several samples 
# pigz: skipping: /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5_alex/194_ivld_4wpc_s1_L5.fq.gz: corrupted -- incomplete deflate data
# pigz: abort: internal threads error


#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane5.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5_alex/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/03_clean-files/lane5_alex/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane5-bbduklog_v2.txt

done






# 19/04/2024


## STAR alignment

#!/bin/bash

ulimit -n 10000
ulimit -n

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane4.txt |
while read line
do
STAR --runThreadN 70 --genomeDir /data/ffi007/01_quantseq/02_star_genome --readFilesCommand zcat \
 --readFilesIn /data/ffi007/01_quantseq/03_data/03_clean-files/lane4_alex/${line}.clean.fq.gz \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/ffi007/01_quantseq/03_data/07_STAR-out/heart_data_ensembl/lane4_alex/${line}_ |& tee -a starlog-lane4-alex.txt \

done


rename 's/_Aligned.sortedByCoord.out/_sorted/' *.bam


## HTseq-count

#!/bin/bash

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane4.txt | \
parallel -j 57 "htseq-count -m intersection-nonempty -i gene_id -t exon -s yes -f bam -r pos \
/data/ffi007/01_quantseq/03_data/07_STAR-out/heart_data_ensembl/lane3_alex/{}_sorted.bam \
/data/ffi007/01_quantseq/03_data/05_useful-files/genomes/Salmo_salar.Ssal_v3.1.111.gtf \
>/data/ffi007/01_quantseq/03_data/08_ht-seq-count/lane4_alex/{}_readcount.txt" |& tee -a htseqlog-lane4-alex.txt


















