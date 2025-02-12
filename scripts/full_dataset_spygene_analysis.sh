# Re-running all sequencing data through Spygene, as a whole

## Moving data into the server. 15/10/2024
### Lane 1. Do not have demuxed files stored for lane 1, so have to move source
rsync -havP 01_lane1-read1.fq.gz ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/01_source-files

### Lane 2
rsync -havP *.fq.gz ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/02_demuxed-files/

### Lane 3
/Volumes/Backup/PhD/QuantSeq_02_unzipped-source-files/03_lane3-october2018/02_demuxed-files  # source
rsync -havP *.fq.gz ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/02_demuxed-files/

### Lane 4
/Volumes/Backup/PhD/QuantSeq_02_unzipped-source-files/04_lane4-october2018/02_demuxed-files  # source
rsync -havP *.fq.gz ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/02_demuxed-files/

### Lane 5
/Volumes/Backup/PhD/QuantSeq_02_unzipped-source-files/05_lane5-october2018/02_demuxed-files  # source
rsync -havP *.fq.gz ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5/

### Lane 6 - FAILED
/Volumes/Backup/PhD/QuantSeq_02_unzipped-source-files/06_lane6-october2018/02_demuxed-files  # source
rsync -havP *.fq.gz ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane6/


### Demultiplexing Lane 1
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane1/ \
-b /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane1/lane1.lostreads.fq.gz \
-s /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane1/demux-summary/lane1.summary.txt \
/data/ffi007/01_quantseq/03_data/05_useful-files/lane1-barcodes.txt \
/data/ffi007/01_quantseq/03_data/01_source-files/01_lane1-read1.fq.gz


# 16/10/2024

### Lane 6
/Volumes/Backup/PhD/QuantSeq_02_unzipped-source-files/06_lane6-october2018/02_demuxed-files  # source
rsync -havP *.fq.gz ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane6/


## Trimming

#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames_lanes1_to_4.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lanes_1_4/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane1_4-bbduklog.txt

done


# 17/10/2024
# input files are misnamed. eomes should be gata3.
# correcting demuxed files to include the correct names and re-running bbduk
mv 136_eomes_6wpc_s1_L2.fq.gz 136_gata3_6wpc_s1_L2.fq.gz
mv 137_eomes_6wpc_s2_L2.fq.gz 137_gata3_6wpc_s2_L2.fq.gz
mv 138_eomes_6wpc_s3_L2.fq.gz 138_gata3_6wpc_s3_L2.fq.gz
mv 139_eomes_6wpc_s4_L2.fq.gz 139_gata3_6wpc_s4_L2.fq.gz
mv 140_eomes_6wpc_s5_L2.fq.gz 140_gata3_6wpc_s5_L2.fq.gz


diff -y -W 70 --suppress-common-lines demuxed_list.txt new_clean_list.txt

146_eomes_4wpc_h6_L2.fq.gz        <
153_eomes_6wpc_h6_L2.fq.gz        <
160_eomes_6wpc_hk6_L2.fq.gz       <
167_eomes_6wpc_s6_L2.fq.gz        <
174_eomes_4wpc_hk6_L2.fq.gz       <
181_eomes_4wpc_s6_L2.fq.gz        <

# these files are missing in the clean list, and their are misnamed in the demultiplexed
# these are all gata3 files that were misnamed from the beginning in the Excel WORKSHEET
# CHANGED ALL THE WRONG EOMES TO GATA3 IN THE WORKSHEET ON 17/10/2024

## CORRECT FILENAMES
146_gata3_4wpc_h6_L2        
153_gata3_6wpc_h6_L2        
160_gata3_6wpc_hk6_L2       
167_gata3_6wpc_s6_L2        
174_gata3_4wpc_hk6_L2       
181_gata3_4wpc_s6_L2       


# ## Trimming missing GATA3 files

# #!/bin/bash
# cat /data/ffi007/01_quantseq/03_data/05_useful-files/missing_gata3.txt |
# while read line
# do
#         /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
#         in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/$line.fq.gz \
#         out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lanes_1_4/missing_GATA3/$line.clean.fq.gz \
#         ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
#         k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a missing_gata3-bbduklog.txt

# done

# I can't be sure that the demultiplexed files were done so using the correct barcodes,
# so I have to re-demultiplex lane 2

# 18/10/2024
### Lane 2. Moving source to Spygene so I can demultiplex it
rsync -havP 03_lane2-read1.fq.gz ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/01_source-files


## Demultiplexing lane 2
### Lane 2
demuxFQ -c -d -e -i -t 1 -l 9 \
-o /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane2 \
-b /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane2/lane2.lostreads.fq.gz \
-s /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane2/lane2.summary.txt \
/data/ffi007/01_quantseq/03_data/05_useful-files/lane2-barcodes.txt \
/data/ffi007/01_quantseq/03_data/01_source-files/03_lane2-read1.fq.gz


# 19/10/2024
## Trimming lane 2
#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane2.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane2/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lane2_v2/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane2-bbduklog.txt

done

# 20/10/2024

fastqc -o /data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lanes_1_4/fastqc -t 60 -f fastq *.fq.gz

## Trimming lane 5
#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane5.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane5-bbduklog.txt

done


## Trimming lane 6
#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-lane6.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane6/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane6/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane6-bbduklog.txt

done

# moved clean lanes 1-4 to 03_clean_files/
# need to fastqc lanes 5 and 6, and multiqc everything

fastqc -o /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane5/fastqc -t 72 -f fastq *.clean.fq.gz
fastqc -o /data/ffi007/01_quantseq/03_data/02_demuxed-files/lane6/fastqc -t 72 -f fastq *.clean.fq.gz

# fastqc'ed and multiqc'ed everything. moved multiqc_report.html to Mac.
# moved all clean files into 03_clean-files
# need to map to ENSEMBL Ssal v3.1 annotation


#!/bin/bash

ulimit -n 10000

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-full-dataset.txt |
while read line
do
STAR --runThreadN 70 --genomeDir /data/ffi007/01_quantseq/02_star_genome --readFilesCommand zcat \
 --readFilesIn /data/ffi007/01_quantseq/03_data/03_clean-files/${line}.clean.fq.gz \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/ffi007/01_quantseq/03_data/07_STAR-out/full_dataset/${line}_ |& tee -a starlog-full-dataset.txt \

done

# 22/10/2024

rename 's/_Aligned.sortedByCoord.out/_sorted/' *.bam

## HTseq-count

#!/bin/bash

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames-full-dataset.txt | \
parallel -j 57 "htseq-count -m intersection-nonempty -i gene_id -t exon -s yes -f bam -r pos \
/data/ffi007/01_quantseq/03_data/07_STAR-out/full_dataset/sorted/{}_sorted.bam \
/data/ffi007/01_quantseq/03_data/05_useful-files/genomes/Ssal_v3.1/ENSEMBL/Salmo_salar.Ssal_v3.1.113.gtf \
>/data/ffi007/01_quantseq/03_data/08_ht-seq-count/{}_readcount.txt" |& tee -a htseqlog_fulldataset.txt


# 06/11/2024
# Moving count data from Spygene to laptop
rsync -havP --info=progress2 "ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/08_ht-seq-count/*_readcount.txt" "/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/full_dataset/readcounts/"


# Forgot to deal with lane 7. Moving source to Spygene from Backup
rsync -havPn --info=progress2 "/Volumes/Backup/PhD/QuantSeq_library3_source-file_January2024/X204SC24010084-Z01-F001.tar" "ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/01_source-files/"

# Extracting source file
tar -xvf X204SC24010084-Z01-F001.tar

## Renaming dry-run
while IFS=$'\t' read -r old_name new_name; do
  # Check if the file with the old name exists
  if [[ -f "$old_name" ]]; then
    # Dry run: display the rename command instead of executing it
    echo "mv \"$old_name\" \"$new_name\""
  else
    echo "File $old_name not found"
  fi
done < renaming_lane7.txt


## Renaming demuxed files named by Novogene
while IFS=$'\t' read -r old_name new_name; do
  # Check if the file with the old name exists
  if [[ -f "$old_name" ]]; then
    # Rename the file to the new name
    mv "$old_name" "$new_name"
    echo "Renamed $old_name to $new_name"
  else
    echo "File $old_name not found"
  fi
done < renaming_lane7.txt

## Removing spaces and dashes
for file in *; do
  # Replace spaces with underscores
  mv "$file" "${file// /_}"
done

for file in *dna_vaccine*.fq.gz; do
  mv "$file" "${file//dna_vaccine/dnavaccine}"
done

for file in *iv-hd*.fq.gz; do
  mv "$file" "${file//iv-hd/ivhd}"
done

for file in *iv-ld*.fq.gz; do
  mv "$file" "${file//iv-ld/ivld}"
done

sed 's/\.fq\.gz$//' lane7_filenames.txt > output_file.txt


## Trimming

#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames_lane7.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane7/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lane7/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane7-bbduklog.txt

done


fastqc -o /data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lane7/fastqc -t 72 -f fastq *.clean.fq.gz


## Mapping
#!/bin/bash

ulimit -n 10000

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames_lane7.txt |
while read line
do
STAR --runThreadN 70 --genomeDir /data/ffi007/01_quantseq/02_star_genome --readFilesCommand zcat \
 --readFilesIn /data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lane7/${line}.clean.fq.gz \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/ffi007/01_quantseq/03_data/07_STAR-out/lane7/${line}_ |& tee -a starlog-lane7.txt \

done

# 07/11/2024

/data/ffi007/01_quantseq/03_data/07_STAR-out/lane7$ multiqc .

rename 's/_Aligned.sortedByCoord.out/_sorted/' *.bam

### Need to redo trimming and mapping because filename_lane7.txt was wrong
sed -i 's/dna vaccine/dnavaccine/g; s/iv-hd/ivhd/g; s/iv-ld/ivld/g' filenames_lane7.txt


## Trimming - lane7

#!/bin/bash
cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames_lane7.txt |
while read line
do
        /data/ffi007/01_quantseq/01_bbmap/bbduk.sh \
        in=/data/ffi007/01_quantseq/03_data/02_demuxed-files/lane7/$line.fq.gz \
        out=/data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lane7/$line.clean.fq.gz \
        ref=/data/ffi007/01_quantseq/03_data/05_useful-files/polyA.fa,/data/ffi007/01_quantseq/03_data/05_useful-files/adapters.fa \
        k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20 |& tee -a lane7-bbduklog.txt

done


fastqc -o /data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lane7/fastqc -t 72 -f fastq *.clean.fq.gz


## Mapping - lane7
#!/bin/bash

ulimit -n 10000

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames_lane7.txt |
while read line
do
STAR --runThreadN 70 --genomeDir /data/ffi007/01_quantseq/02_star_genome --readFilesCommand zcat \
 --readFilesIn /data/ffi007/01_quantseq/03_data/02_demuxed-files/clean_lane7/${line}.clean.fq.gz \
 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 \
 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD \
 --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/ffi007/01_quantseq/03_data/07_STAR-out/lane7/${line}_ |& tee -a starlog-lane7.txt \

done

rename 's/_Aligned.sortedByCoord.out/_sorted/' *.bam

## HTseq-count - lane7

#!/bin/bash

cat /data/ffi007/01_quantseq/03_data/05_useful-files/filenames_lane7.txt | \
parallel -j 57 "htseq-count -m intersection-nonempty -i gene_id -t exon -s yes -f bam -r pos \
/data/ffi007/01_quantseq/03_data/07_STAR-out/lane7/sorted/{}_sorted.bam \
/data/ffi007/01_quantseq/03_data/05_useful-files/genomes/Ssal_v3.1/ENSEMBL/Salmo_salar.Ssal_v3.1.113.gtf \
>/data/ffi007/01_quantseq/03_data/08_ht-seq-count/lane7/{}_readcount.txt" |& tee -a htseqlog_lane7.txt

/data/ffi007/01_quantseq/03_data/08_ht-seq-count/lane7$ multiqc .

ffi007@uit-mac-1312 ~ % scp 'ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/08_ht-seq-count/lane7/*_readcount.txt' \
/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/full_dataset/











