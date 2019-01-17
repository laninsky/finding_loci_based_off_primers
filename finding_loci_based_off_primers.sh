genomefile=`tail -n+1 params.txt | head -n1`;
primerfile=`tail -n+2 params.txt | head -n1`;
tolerance=`tail -n+3 params.txt | head -n1`;

makeblastdb -in $genomefile -dbtype nucl

noloci=`wc -l $primerfile | awk '{print $1}'`;

bylocus () {
  locusline=`tail -n+$i $primerfile | head -n1`;
  firstseq=`echo $locusline | awk '{print $2}'`
  secondseq=`echo $locusline | awk '{print $3}'`
  echo ${firstseq}"NNNNNNNNNN"${secondseq} > $i.queryseq

  blastn -task blastn-short -db $genomefile -query $i.queryseq -outfmt '6 qseqid sseqid sstart send bitscore evalue' > $i.blast
}

# for-loop allowing bylocus to be parallelized
for i in `seq 1 $noloci`;
  do bylocus "$i" & done
wait

Rscript summarising_blast.R

rm *queryseq

noloci=`wc -l best_blast_matches.txt | awk '{print $1}'`;

samtools faidx $genomefile

for i in `seq 2 $noloci`;
  # extracting scaffold and locus names
  do locusline=`tail -n+$i best_blast_matches.txt | head -n1`;
  scaffoldname=`echo $locusline | awk '{print $2}'`;
  locusname=`echo $locusline | awk '{print $1}'`;
  
  # extracting genomic coordinates, adding on padding, and making finding
  # the maximum length of the scaffold
  tempstartpos=`echo $locusline | awk '{print $6}'`;
  tempendpos=`echo $locusline | awk '{print $7}'`;
  maxscaffoldlength=`grep -A1 $scaffoldname $genomefile | tail -n 1 | wc | awk '{print $3}'`;
  locuslength=`grep $locusname $primerfile | awk '{print $4}'`;  
  tempstartpos=$((tempstartpos-locuslength-tolerance));
  tempendpos=$((tempendpos+locuslength+tolerance));

  # if the 'padded' start position is less than 1, resetting to 1
  if (( 1 > $tempstartpos )); then tempstartpos=1; fi;
  # if the 'padded' end position is larger than maxscaffoldlength, resetting to maxscaffoldlength  
  if (( $tempendpos > $maxscaffoldlength )); then tempendpos=$maxscaffoldlength; fi;
  
  # printing sequence name to file
  echo ">"${locusname}.${scaffoldname}.${tempstartpos}.${tempendpos} >> best_blast_matches.fa;
  
  # printing genomic slice to file
  samtools faidx $genomefile $scaffoldname:$tempstartpos-$tempendpos | sed 1,1d >> best_blast_matches.fa;
done  
