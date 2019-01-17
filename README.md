# finding_loci_based_off_primers
Uses BLAST and R scripts to find genomic position of loci based off their primer sequences and expected size

### Things you need
1. A genome in fasta format
2. A tab/space-delimited list of your primers with locus name in the 1st column, F primer (5' to 3') in the 2nd column, R primer (5' to 3') in the 3rd column, and expected size including primers (bp) in the 4th. e.g.
```
2F91	GCATTTCTGGGCTGTAACAT	 AAAGGACAATGTAATTGGTG 		99
4E8	TGTGGGAAACCAGAGGAAA	 CAGGGGAAAAATAGAGAGGG 		99
4H2	GTTTTGACCTACAGAGAGAG	 GATAGGACAATCAACAGGCT 			214
ADCbm	GATGTGAGTAACCAGCCACT	 ATAACACAGGAGCGGTGA 		177
Ase64	CCACCTTTCATACTGGGGAG	 TTCAGCCAGTCAGTGTAGCC 			399
ASu15	AATAGATTCAGGTGCTTTTTCC	 GGTTTTTGAGAAAATTATACTTTCAG 			124
CcaTgu21	GGCAGACATGATTGCATCC	 TCTCAGTGGTCATTGGAAAGTG 			205
DkiB119	CATACAACTTCATGACTACCATAGCAC	 TCCATAGTGACATAGAACGAGCTG 			230
Escu6	CATAGTGATGCCCTGCTAGG	 GCAAGTGCTCCTTAATATTTGG 			110
GCSW42	GGCTTCTCTGGTTGCATGTC	 ACAGTAATCCCCAGCCATCA 			215
Ind28	CCCAGGAAGTATCCCAGAA	 CCTCCAATGCTTTAGTGACC 			166
Mjg1	CCCGGGAAAGGCTTCGTCTTC	 GGAGATTTTATATCGGTGGC 			172
Pau01	TTAGAAGTGAAAGGCTTG	 GAGGAATAAAAACAATGC 			110
Pau04	AATAAAGCAGATACTGAG	 ACAGGTAAACCAGAGCAG 			131
Pau06	CTGAGGTTCAAAGTTTCC	 ACCAGCCATCCTTATGC 			102
Pau07	CTTTCCTTGACTGAAGTG	 TTAGCTTCATTTCCAGTC 			106
Pau09	ATGATGTAGTCAGAGTCG	 TATTTTGCAACCTTCTTG 			114
Pau23	CAGGGCATTTACAGATTTTCC	 CTCATCTTGCACACTGCTG 			188
Pau26	AAATTACTACAGTGTTACGGTGAAAA	 GGGACCACCAAGAACTTCAA 			170
Pau66	TGGGCCAGTTTATACCCTCT	 ATGAAAGGGTTCCATGATGC 		106
Pau67	CCAGGAAAGGTGCTCAGAGT	 TGTCTGTGTTGGCCTGATCT 			181
Pocc6	TCACCCTCAAAAACACACACA	 ACTTCTCTCTGAAAAGGGGAGC 			184
Ppi2	CACAGACCATTCGAAGCAGA	 GCTCCGATGGTGAATGAAGT 			324
TG02	TGTGTGTTGACAGTATTCTCTTGC	 TTTAAACCTAATAAACGTCACACAGTC 			260
TG11	ACAAACTAAGTACATCTATATCTGAAG	 TAAATACAGGCAACATTGG 			224
Tgu01	TGCGGTCTGTATGGAAATAGTC	 CTTGCAATACTCTCTGCCTCA 			192
Tgu03	TCTCTCTGCTAGGGATAAACAGTG	 TGCTCCCTCCCTCCAGTAAC 			158
TguEST09	AACCCAACCAACAAAATTGG	 CCAACTATCAGTTTTACAAGGCATAC 			173
```
3. BLAST, samtools and R (with the libraries dplyr and readr) installed and in your path
4. An idea of the tolerance around the size of locus you want to detect (e.g. for our microsatellite example, we'd expect an individual to be within 50 bp of the average locus size)
5. A params.txt file with the name of your genome in the first line, the name of your tab-delimited list of loci (see #2) in the second line, and the 'slush' in length of loci you are interested in (see #4) e.g.
```
blob_genome_15Jan2019.fa
primers.txt
50
```

### How do you run it?
Make sure the BLAST, samtools and R (with dplyr and readr previously installed) are in your $PATH. Make sure finding_loci_based_off_primers.sh and summarising_blast.R are in the same working directory as params.txt, and then execute by:
```
bash finding_loci_based_off_primers.sh
```

### What do you get?
_**best_blast_matches.txt**_: A tab-delimited output file with the best blast matches for each locus and their genomic location. The "best" scaffolds for each locus are evaluated by their matchrank (between 0 and 4), from adding up the following points:
1. nomatches: if a locus matches to the scaffold 2 or more times, the scaffold gets 1 point, otherwise 0 (we expect 2 matches for 'good' locations - corresponding to each of the two primers matching to the scaffold)
2. max_evalue: this column summarises the max_evalue across the locus matches for a given scaffold. If this is the smallest e-value across scaffolds for that locus then that scaffold wins a point.
3. min_bitscore: this columns summarised the min_bitscore for the matches for a given scaffold. If this is the largest bit score seen across scaffolds for that locus then that scaffold wins 1 point, otherwise 0.
4. Length within the acceptable bounds (from the tab-delimited file, plus/minus slush) gets 1 point, otherwise 0.
```
locus_name scaffold_name nomatches max_evalue min_bitscore startpos endpos length matchrank
2F91 scaffold274_size858921 2 0.028 40.1 818015 818110 96 4
4E8 scaffold280_size1803993_1_925356 1 0.44 36.2 801410 801427 18 2
4H2 scaffold16_size8633371 1 6.7 32.2 5067374 5067389 16 2
4H2 scaffold18_size4552767 1 6.7 32.2 1347858 1347873 16 2
4H2 scaffold427_size884515_1_848045 1 6.7 32.2 575045 575060 16 2
```
For each locus, the scaffolds with the highest equal matchrank are written out. We'll next pull out the fasta for these matches (with some slush) so you can evaluate whether you think these are good matches or not.

_**best_blast_matches.fa**_: a fasta file containing the the genomic sequences (with padding on each size equal to expected locus length + slush). Sequences are named: locus_name.scaffold_name.startpos.endpos


Credit to IDT (https://sg.idtdna.com/pages/education/decoded/article/tips-for-using-blast-to-locate-pcr-primers) for the idea behind this pipeline
