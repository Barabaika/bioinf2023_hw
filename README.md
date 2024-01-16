# Bioinformatics 2023 homeworks

**before start:**
 - create dir `./data` 
 - unzip `data_archive.zip` contents to data
 - for hw3, need to load data (detales are in HW 3 section)

**requirements**
 - biopython (`pip install biopython`)
 - tqdm (`pip install tqdm`)


## HW 1:
**answer:** in data/out_hw1.txt \
**reproduse:** `python hw_1_contigs.py data/Carsonella_ruddii_reads_paired_reads_left.fastq data/Carsonella_ruddii_reads_paired_reads_right.fastq -o hw1_test_out.txt -k 51` \
**notes:** Very slow: :(. about 1500 minutes to run. In output there are meny contigs == reads and only about 3 - longer sequences

## HW 2:
**answer:** YNMFPQVYEP, for N = 10000 \
**score:** 37/90 \
**reproduse:** `python hw_2_massspec.py data/Spectrum_task2.txt -n 20 --cyclic` \
**notes:** tried N from 10 to 10000, with N grows we able to find proteins with better score.To improve- if there are multiple sequences with equal score i output only one score, need to output all
- N=    10 - score = 23
- N=  1000 - score = 36
- N= 10000 - score = 37

## HW 3:
load data and store it to ./data as chrX_human.fa and chrX_mouse.fa: \
`wget -O data/chrX_human.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrX.fa.gz ` \
`wget -O data/chrX_mouse.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/chrX.fa.gz` \
`gzip -d chrX_human.fa.gz`
`gzip -d chrX_mouse.fa.gz`
**code to run full data:** `python hw_3_2break_dist.py -k 50 --max_distance 3 --min_syntency_block_size 10` \
**code to run test sequence:** `python hw_3_2break_dist.py -k 3 --max_distance 2 --min_syntency_block_size 1 --test` \
**notes:** can not process full sequences (150 mln letters). Run out of memory

