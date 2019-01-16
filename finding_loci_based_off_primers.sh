genomefile=`tail -n+1 params.txt | head -n1`;
primerfile=`tail -n+2 params.txt | head -n1`;
tolerance=`tail -n+3 params.txt | head -n1`;

makeblastdb -in $genomefile -dbtype nucl

noloci=`wc -l $primerfile | awk '{print $1}'`;

bylocus () {
locusline=`tail -n+i $primerfile | head -n1`;
firstseq=`echo $locusline | awk '{print $2}'`
secondseq=`echo $locusline | awk '{print $3}'`
echo ${firstseq}"NNNNNNNNNN"${secondseq} > $i.queryseq

blastn -task blastn-short -db $genomefile -query $i.queryseq -outfmt '6 qseqid sseqid sstart send bitscore' > $i.blast
}

# for-loop allowing bylocus to be parallelized
for i in `seq 1 $noloci`;
do bylocus "$i" & done
wait

Rscript summarising_blast.R
