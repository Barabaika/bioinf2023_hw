# Bioinformatics 2023 homeworks

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
**code to run full data:** `python hw_3_2break_dist.py -k 50 --max_distance 3 --min_syntency_block_size 10` \
**code to run test sequence:** `python hw_3_2break_dist.py -k 50 --max_distance 3 --min_syntency_block_size 10 --test` \
data loaded from: \
human - https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/ \
mouse - https://hgdownload.soe.ucsc.edu/goldenPath/mm39/chromosomes/ \
**notes:** 

