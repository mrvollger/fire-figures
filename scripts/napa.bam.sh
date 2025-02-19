 samtools merge \
 --write-index \
 -o data/gm12878.napa.bam \
 -L Tables/NAPA.100k.bed \
 -b scripts/napa.gm12878.bam.fofn \
 -f

