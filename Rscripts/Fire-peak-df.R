source("Rscripts/utils.R")
#library(regioneR)
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)

if(F){

    # VCF stuff
    vcf="data/GM12878.vcf.gz"
    vcf_df = fread(cmd=glue("bcftools view -f 'PASS' -i 'GT=\"HET\"' {vcf} | grep -v ^##"))  %>% 
        mutate(start = POS-1) %>%
        rename(
            chrom="#CHROM",
            end=POS
        )

    pbsv = "data/GM12878.pbsv.vcf.gz"
    KVsep <- fixed(";")  #key-value separator
    pbsv_vcf = fread(cmd=glue("bcftools view -f 'PASS' -i 'GT=\"HET\"' {pbsv} | grep -v ^##")) %>%
    mutate(start = POS-1) %>%
        rename(
            chrom="#CHROM",
        ) %>%
        select(-ALT) %>%
        mutate(KVpairs = str_split(INFO, KVsep)) %>%
        unnest(KVpairs) %>%
        separate(KVpairs, into = c("key", "value"), sep="\\=") %>%
        spread(key, value) %>%
        mutate(
            end = as.numeric(END)
        )
    pbsv_ends = pbsv_vcf %>%
        pivot_longer(
            cols=c("start", "end")
        ) %>%
        select(chrom,name,value) %>%
        mutate(
            start = value, 
            end = start + 1,
        ) %>%
        filter(!is.na(end))
        #%>%unnest()
    ### END VCF stuff

    FAI = fread("~/assemblies/hg38.analysisSet.chrom.sizes") %>%
        filter(!grepl("_", V1))
    colnames(FAI) = c("chrom", "end")
    FAI$chrom = factor(FAI$chrom, levels=FAI$chrom)
    FAI$name = FAI$chrom
    FAI$size = FAI$end
    FAI$start= 0

    FAI_NO_SEX = FAI %>% filter(!grepl("chr[XYM]", chrom)) %>% filter(chrom!="chrEBV")
    FAI_NO_SEX

    # annotations 
    unreliable_df = my_read_bed("results/GM12878/coverage/unreliable-coverage-regions.bed.gz") %>% bed_merge()
    blacklist_df = my_read_bed("data/encode_blacklist_ENCFF356LFX.bed.gz") %>% bed_merge()
    imprinted=my_read_bed("data/lcl_dmr_coordinates_Akbari.bed.gz")
    ALTs = my_read_bed("data/GRCh38-alt-locations.bed.gz")
    sds=my_read_bed("data/SDs.merged.hg38.bed.gz")
    SDs=sds
    SDs_all = fread("data/SDs.bed.gz")
    colnames(SDs_all)[1:4] = c("bin", "chrom", "start", "end")
    SD_size = sum(SDs$end - SDs$start)
    G_size = sum(FAI$end - FAI$start)
    encode=my_read_bed("data/ENCODE3_consensus_DHS_ENCFF503GCK.tsv.gz")
    tss=my_read_bed("data/gencode.v42.annotation_TSS.gff3.gz")
    colnames(tss)[1:6] = c("chrom", "start", "end", "name", "score", "strand")
    cage_df = my_read_bed("data/GM12878_cage_peaks_merged.bed.gz")
 
    cage_tss_with_direction_df = fread("data/TSS-CAGE-FIRE-intersect.bed.gz")
    colnames(cage_tss_with_direction_df)[1:6] = c("chrom", "start", "end", "name", "score", "strand")
    
    # data 
    dnase_peaks=my_read_bed("data/ENCFF762CRQ_DNase_peaks.bed.gz")
    dnase=my_read_bed("data/bedgraph_annotations/ENCFF960FMM_dnase_signal.bed.gz")
    colnames(dnase)[4] = "dnase_sig"
    
    #atac_peaks = my_read_bed("data/10X_GM12878_peaks_max_cov.bed.gz")
    atac_peaks = my_read_bed("data/scATAC_GM12878_peaks_MACS2.bed.gz")
    atac = my_read_bed("data/10X_GM12878_aggr_scATAC.bg.gz")
    colnames(atac)[4] = "atac_sig"

    # 5mc 
    #chr1    10468   10469   92.4    hap1    43      40      3       93.0
    cpg_df = my_read_bed("data/GM12878-5mC.hap.bed.gz");
    colnames(cpg_df) = c("chrom", "start", "end", 
        "raw_hap1_5mC_percent", "hap1_5mC_percent",
        "raw_hap2_5mC_percent", "hap2_5mC_percent"
    )
    

    # CTCF 
    c1=my_read_bed("data/CTCF_peak_ENCFF356LIU.bed.gz")
    c2=my_read_bed("data/CTCF_peak_ENCFF960ZGP.bed.gz")
    ctcf_peaks_df = bind_rows(c1,c2) %>% bed_merge(
        max_ctcf = max(V7)
    )
    ctcf_peaks_df

    ctcf_motifs = my_read_bed("data/ctcf-motifs.bed.gz")


    in_file="results/GM12878/FDR-peaks/FDR-FIRE-peaks.bed.gz"
    df=fread(in_file)
    df$peak_cov = df$coverage
    df$acc_percent = df$fire_coverage/df$peak_cov
    colnames(df)[1:5]=c("chrom", "start","end","ostart","oend")
    df$ID = paste0(df$chrom, "_", df$start, "_", df$end)

    vcf="data/GM12878.vcf.gz"
    gm12878_vcf_df = fread(cmd=glue("bcftools view -f 'PASS' -i 'GT=\"HET\"' {vcf} | grep -v ^##"))  %>% 
    mutate(start = POS-1) %>%
    rename(
        chrom="#CHROM",
        end=POS
    )

    df = df %>%
        filter(acc_percent >= MIN_FRAC_ACC) %>%
        filter(pass_coverage) %>%
        arrange(-acc_percent) %>%
        mutate(
            n=seq(n()),
            min_percent_acc=min(acc_percent)
        ) %>%
        arrange(-n) %>%
        mutate(
            Mbp=cumsum(end-start)/1e6,
        ) %>%
        bed_map(encode,
            encode_count=length(core_end),
            encode_anno = paste0(unique(component), collapse=";")
        ) %>%
        bed_map(sds, sd_count=length(end)) %>%
        bed_map(tss,
            TSS=length(name),
            is_TSS = n()>0,
            TSS_strand = case_when(
                length(unique(strand)) == 1 ~ unique(strand)[1],
                TRUE ~ ".",
            ),
        ) %>%
        bed_map(dnase_peaks, 
            is_dnase_peak = n()>0,
        ) %>%
        bed_map(ctcf_peaks_df, 
            is_ctcf_peak = n()>0,
            ctcf_max = max(max_ctcf)[[1]],
        ) %>%
        bed_map(
            ctcf_motifs,
            has_ctcf_motif = n()>0,
        ) %>%
        bed_map(
            imprinted,
            imprinted=n()>0
        ) %>%
        bed_map(
            unreliable_df,
            is_unreliable = n()>0,
        ) %>%
        bed_map(
            blacklist_df,
            is_blacklist = n()>0,
        ) %>%
        bed_map(
            ALTs,
            is_alt = n()>0,
        ) %>%
        bed_map(
            SDs,
            is_SD = n()>0,
        ) %>%
        bed_map(
            cage_df,
            is_cage_peak = n()>0,
        ) %>%
        bed_map(atac_peaks, 
            is_atac_peak = n()>0,
        ) %>%
        bed_map(
            SDs_all,
            SD_max_frac_match = max(fracMatch),
        ) %>% 
        bed_map(
            gm12878_vcf_df,
            is_het = n()>0,
            n_snvs = n(),
        ) %>%
        bed_map(
            cage_tss_with_direction_df,
            is_cage_tss = n()>0,
            cage_tss_strand = case_when(
                length(unique(strand)) == 1 ~ unique(strand)[1],
                TRUE ~ ".",
            ),
        ) %>%
        bed_map(cpg_df,
            max_hap_5mC_diff = max(abs(hap1_5mC_percent - hap2_5mC_percent))[1],
            max_hap1_5mC_percent = max(hap1_5mC_percent)[1],
            max_hap2_5mC_percent = max(hap2_5mC_percent)[1],
            ave_hap_5mC_diff = mean(abs(hap1_5mC_percent - hap2_5mC_percent))[1],
            ave_hap1_5mC_percent = mean(hap1_5mC_percent)[1],
            ave_hap2_5mC_percent = mean(hap2_5mC_percent)[1],
        ) %>%
        bed_map(atac, atac_max = max(atac_sig)) %>%
        bed_map(dnase, dnase_max = max(dnase_sig)) %>%
        bed_map(vcf_df,
            n_SNVs = n(),
        ) %>% 
        bed_map(pbsv_ends,
            n_SVs = n(),
        ) %>%
        replace_na(
            list(
                encode_count = 0,
                sd_count=0,
                n_SNVs=0,
                n_SVs=0,
                TSS=0,
                is_TSS=F,
                is_atac_peak=F,
                is_ctcf_peak=F,
                is_dnase_peak=F,
                has_ctcf_motif=F,
                imprinted=F,
                is_alt=F,
                is_SD=F,
                is_blacklist=F,
                is_cage_peak=F,
                is_unreliable=F,
                is_cage_tss=F
            )
        ) %>%
        arrange(-acc_percent, -n) %>%
        mutate(
            has_snv = n_SNVs > 0,
            has_sv = n_SVs > 0,
            has_variant = has_snv | has_sv,
            group = case_when(
                floor(acc_percent * 20) / 20 >= 0.9 ~ 0.9,
                TRUE ~ floor(acc_percent * 20) / 20,
            )
        ) %>%
        data.table()
    unique(df$group)
    table(df$TSS_strand)
    table(df$cage_tss_strand)

    df = df %>%
        mutate(
            re_strand = case_when(
                is.na(TSS_strand) & !is.na(cage_tss_strand) ~ cage_tss_strand,
                !is.na(TSS_strand) & is.na(cage_tss_strand) ~ TSS_strand,
                (cage_tss_strand != TSS_strand) ~ cage_tss_strand,
                TRUE ~ cage_tss_strand,
            ), 
            #re_stand = cage_tss_strand,
            downstream_ssd = case_when(
                re_strand == "-" ~ FIRE_start_ssd,
                re_strand == "+" ~ FIRE_end_ssd,
                TRUE ~ NA,
            ),
            upstream_ssd = case_when(
                re_strand == "+" ~ FIRE_start_ssd,
                re_strand == "-" ~ FIRE_end_ssd,
                TRUE ~ NA,
            ),
        )
    #cond = df$is_cage_tss & (df$cage_tss_strand != ".") & (df$is_ctcf_peak == F) & df$group > 0.5
    
    #Build a GRanges from your matrix
    ranges <- toGRanges(df[,c("chrom", "start", "end")])

    #Get the sequences and compute the GC content
    freqs = alphabetFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, ranges))
    df$GC_frac = (freqs[,'C'] + freqs[,'G'])/rowSums(freqs)

    # GC correction
    df$psize=df$end-df$start
    fit = lm(acc_percent ~ `GC_frac` * log10(dnase_max), data=df[is_dnase_peak==T,])
    summary(fit)
    df$dnase_max_gc_corrected = predict(fit, newdata=df)

    fit = lm(acc_percent ~ log10(atac_max), data=df[is_atac_peak==T,])
    summary(fit)
    fit = lm(acc_percent ~ `GC_frac` * log10(atac_max), data=df[is_atac_peak==T,])
    summary(fit)
    df$atac_max_gc_corrected = predict(fit, newdata=df) 

    # add hap specific peaks 
    hap_peaks = my_read_bed("results/GM12878/hap1-vs-hap2/FIRE.hap.differences.bed") %>%
        filter(fire_coverage/coverage >= MIN_FRAC_ACC) %>%
        mutate(
            p_adjust = p.adjust(p_value, method="BH"),
        ) %>%
        dplyr::select(
            chrom, start, end, 
            hap1_frac_acc, hap1_acc, hap2_acc, hap2_frac_acc, hap1_nacc, hap2_nacc,
            p_value, p_adjust, diff
        )

    fire_df = merge(df, hap_peaks, by=c("chrom", "start", "end"), all.x=T)

    system("mkdir -p Rdata")
    fire_df %>%
        fwrite("Rdata/fire-peaks-with-annotations.tsv.gz", sep="\t", quote=F)        
    con <- pipe("pigz -p8 > Rdata/df.fire-peaks.gz", "wb")
    save(
        fire_df,
        unreliable_df,
        blacklist_df,
        encode,
        dnase_peaks,
        atac_peaks,
        ctcf_peaks_df,
        ctcf_motifs,
        sds,
        pbsv_vcf,
        pbsv_ends,
        vcf_df,
        tss,
        imprinted,
        SDs,
        ALTs,
        FAI,
        FAI_NO_SEX,
        file = con
    ); close(con)
} else{
    load("Rdata/df.fire-peaks.gz")
}

