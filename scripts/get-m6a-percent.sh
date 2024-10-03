cat gm12878.bam.fofn | parallel -n 1 \
  $'printf \'{/.}\t\'; \
    samtools view -s 0.1 {} -@ 16 -u | ft extract --all - -s | hck -F fiber -F total_m6a_bp -F total_AT_bp \
    | sed -r \'s#/[0-9]+/ccs##g\' \
    | head -n 100000 \
    | sort \
    | datamash -H -g 1 sum 2 sum 3 count 1 \
    | grep -v Group'
