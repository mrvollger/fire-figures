
mkdir -p ../data/ctcf-motifs/
for x in MA0139.1 MA1102.2 MA1929.1 MA1930.1; do
	wget -O ../data/ctcf-motifs/$x.tsv.gz http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022//hg38/$x.tsv.gz
done

zcat ../data/ctcf-motifs/*tsv.gz | bedtools sort | bedtools merge | bgzip > ../data/ctcf-motifs.bed.gz

