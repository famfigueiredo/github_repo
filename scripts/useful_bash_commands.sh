https://explainshell.com/

### Useful bash commands ###

# Cut/Paste current command behind cursor
Ctrl + U, Ctrl + Y

### SSH keys
ls ~/.ssh/*.pub # list all your public keys 
ssh-agent sh -c 'ssh-add; ssh-add -l' # RSA fingerprint
ssh-agent sh -c 'ssh-add; ssh-add -L' # public key
pbcopy < ~/.ssh/id_rsa.pub # copy RSA key to notepad

ssh ffi007@spygene.uit.no # login to Spygene

### File transfer from remote to local
rsync -havP $src $dst
# rsync arguments ... -h: human-readable format; -P: combines --progress and --partial; -n: dry run; -a: archive mode; -v: verbose; -z: compression; -e: specify the remote shell to use (i.e. ssh); --info: progress bar
rsync -havP . ffi007@spygene.uit.no:/data/ffi007/01_quantseq/03_data/january2020/02_unzipped-source-files/02_lane2-december2019/02_demuxed-files

### Kill all screen sessions
killall screen # DO NOT SUDO! #

### Check HD space
du -lh  # directory size
df -lh  # system-wide space availability
du -h -d1 # level-specific directory size

# Another way to get a report about the disk usage of the first-level subdirectories is to use the --max-depth option:
du -h --max-depth=1 /data

### Refresh path
hash -r

### Check versions of software in your PATH
which -a software

### In the case of needing to use an older version of Java #
sudo apt-get update
sudo apt-get openjdk-8-jre # installing older version of Java #
sudo update-alternatives --config java # adding older version of Java to root path #

### Creating soft links to data #
ln -s /PATH/TO/ORIGINAL/FOLDER/*.fq.gz . # This creates a soft link of all the .fq.gz files in original directory, in current directory #

### Changing directory ownership #
sudo chown -R ffi007:ffi007 . # Changing ownership of all current directories (-R) #
sudo chown -R ffi007:science folder/

### Permission numbers. Order goes user/group/users #
    0 = ---
    1 = --x
    2 = -w-
    3 = -wx
    4 = r-
    5 = r-x
    6 = rw-
    7 = rwx


### Remember to change permissions for bash script execution #
chmod 755 bashcript.sh
chmod +x bashcript.sh  # to make it executable

### And for the filelist.txt #
chmod 755 filelist.txt

### Find a file #
find /data -name uniprot_modelfish_backtranseq.fasta # find the fasta file in the /data directory and all sub-directories #

### Print line X of GFF file #
sed -n Xp GFF_FILE.gff

### Add text before first line in text file #
sed -i '1 i\gene_name transcript_id FPKM_CHN' *.txt

### Add .fq.gz to each line in filenames.txt #
sed -i 's/$/.fq.gz/' filenames.txt

### Selecting lines of interest in a GFF file #
grep transcript_id inputGFFfile > outputGFFfile # all lines containing 'transcript_id' #
grep -A20 -B20 target file.txt  # 20 lines before and after the target

### Counting lines in a file #
wc -l someGFFfile.gff

### Show the first line (column headers) in an organized manner #
cat *  |  head -n 1 | tr “\t” “\n” | nl

### Compare two files #
diff -y -W 70 --suppress-common-lines alpha1 alpha2  # comparing files alpha 1 and 2, while suppressing common lines

### Selecting lines from a txt file #
awk '$2 > 0' FILE.txt > SOMEOUTPUTFILE.txt # Selecting lines where column 2 has a value > 0 #

### Print selected columns from a txt file #
awk '{print $x,$y}' inputfile > outputfile # For x and y columns #

### Reading a few lines from a .gz file #
zcat file.fq.gz | head

### Read count in fq.gz file #
zcat file.fq.gz | echo $((`wc -l`/4))
# or
zcat file.fq.gz | echo $(( $(wc -l)/4 ))

### Screen management #
screen -r -d sessioname # reopen attached sessions #
screen -S sessioname # create session with custom name #
screen -X -S sessioname kill # kill screen session #

### Rename files large number of files #
rename 's/abc/xyz/' *.bam # Rename .bam files where 'abc' shows up, change it to 'xyz' #
rename 's/old_pattern/new_pattern/' * | echo  # dry run rename

### Find a specific file/program #
find / -iname multiqc # might have to sudo #

### Indexing bam file using samtools #
samtools index bamfile.bam

### View bam file #
samtools view -h bamfile.bam

### Delete all files in a directory except important.sh
rm !(important.sh)


### To write the output of a command to a file, there are basically 10 commonly used ways.
# Overview:

#           || visible in terminal ||   visible in file   || existing
#   Syntax  ||  StdOut  |  StdErr  ||  StdOut  |  StdErr  ||   file
# ==========++==========+==========++==========+==========++===========
#     >     ||    no    |   yes    ||   yes    |    no    || overwrite
#     >>    ||    no    |   yes    ||   yes    |    no    ||  append
#           ||          |          ||          |          ||
#    2>     ||   yes    |    no    ||    no    |   yes    || overwrite
#    2>>    ||   yes    |    no    ||    no    |   yes    ||  append
#           ||          |          ||          |          ||
#    &>     ||    no    |    no    ||   yes    |   yes    || overwrite
#    &>>    ||    no    |    no    ||   yes    |   yes    ||  append
#           ||          |          ||          |          ||
#  | tee    ||   yes    |   yes    ||   yes    |    no    || overwrite
#  | tee -a ||   yes    |   yes    ||   yes    |    no    ||  append
#           ||          |          ||          |          ||
#  n.e. (*) ||   yes    |   yes    ||    no    |   yes    || overwrite
#  n.e. (*) ||   yes    |   yes    ||    no    |   yes    ||  append
#           ||          |          ||          |          ||
# |& tee    ||   yes    |   yes    ||   yes    |   yes    || overwrite
# |& tee -a ||   yes    |   yes    ||   yes    |   yes    ||  append

# List:

    command > output.txt

    # The standard output stream will be redirected to the file only, it will not be visible in the terminal. If the file already exists, it gets overwritten.

    command >> output.txt

    # The standard output stream will be redirected to the file only, it will not be visible in the terminal. If the file already exists, the new data will get appended to the end of the file.

    command 2> output.txt

    # The standard error stream will be redirected to the file only, it will not be visible in the terminal. If the file already exists, it gets overwritten.

    command 2>> output.txt

    # The standard error stream will be redirected to the file only, it will not be visible in the terminal. If the file already exists, the new data will get appended to the end of the file.

    command &> output.txt

    # Both the standard output and standard error stream will be redirected to the file only, nothing will be visible in the terminal. If the file already exists, it gets overwritten.

    command &>> output.txt

    # Both the standard output and standard error stream will be redirected to the file only, nothing will be visible in the terminal. If the file already exists, the new data will get appended to the end of the file..

    command | tee output.txt

    # The standard output stream will be copied to the file, it will still be visible in the terminal. If the file already exists, it gets overwritten.

    command | tee -a output.txt

    # The standard output stream will be copied to the file, it will still be visible in the terminal. If the file already exists, the new data will get appended to the end of the file.

    (*)

    # Bash has no shorthand syntax that allows piping only StdErr to a second command, which would be needed here in combination with tee again to complete the table. If you really need something like that, please look at "How to pipe stderr, and not stdout?" on Stack Overflow for some ways how this can be done e.g. by swapping streams or using process substitution.

    command |& tee output.txt

    # Both the standard output and standard error streams will be copied to the file while still being visible in the terminal. If the file already exists, it gets overwritten.

    command |& tee -a output.txt

    # Both the standard output and standard error streams will be copied to the file while still being visible in the terminal. If the file already exists, the new data will get appended to the end of the file.


### Removing files
    # To remove a folder with all its contents (including all interior folders):

    rm -rf /path/to/directory

    # To remove all the contents of the folder (including all interior folders) but not the folder itself:

    rm -rf /path/to/directory/*

   # or, if you want to make sure that hidden files/directories are also removed:

    rm -rf /path/to/directory/{*,.*}

    # To remove all the "files" from inside a folder(not removing interior folders):

    rm -f /path/to/directory/{*,.*}

    ## Delete all file except file1 ##
    rm  !(file1)
     
    ## Delete all file except file1 and file2 ##
    rm  !(file1|file2) 
     
    ## Delete all file except all zip files ##
    rm  !(*.zip)
     
    ## Delete all file except all zip and iso files using '*' wildcard ##
    rm  !(*.zip|*.iso)
     
    ## You set full path too ##
    rm /Users/ffi007/!(*.zip|*.iso|*.mp3)
     
    ## Pass options to the rm ##
    rm [options]  !(*.zip|*.iso)
    rm -v  !(*.zip|*.iso) #verbose option
    rm -f  !(*.zip|*.iso) #force option
    rm -v -i  !(*.php)  #confirm and verbose option



for x in *L6_Aligned.sortedByCoord.out.bam
do
    mv "$x" "${x%_Aligned.sortedByCoord.out.bam}"
done

# To download FASTA files from NCBI. Go to page, view source, copy 'ncbi_uidlist', and paste to query string.
# Use 

wget -O gata3.fq \
'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=180038037&db=nuccore&report=fasta&retmode=text&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=100000000'


# File transfer to FTP
curl -k --ssl-reqd -u uit.temp:'PASSWORD' -T test.txt ftp://neptuno.sparos.pt/  # this worked

for file in /data/ffi007/01_quantseq/03_data/02_demuxed-files/*.fq.gz; do
  curl -k --ssl-reqd -u uit.temp:'PASSWORD' -T "$file" ftp://neptuno.sparos.pt/
done

curl -k --ssl-reqd -v -u uit.temp:'PASSWORD' -T test_sample.fq.gz ftp://neptuno.sparos.pt/