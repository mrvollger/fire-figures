ls results/*/fire/*bam | parallel -n 1 -k 'printf "{.}\t"; samtools view -c -@ 8 {}' > Tables/fiber-counts.txt
ls results/*/fiber-calls/FIRE.bed.gz | parallel -n 1 -k 'printf "{.}\t"; zcat {} | wc -l' > Tables/fire-counts.txt
