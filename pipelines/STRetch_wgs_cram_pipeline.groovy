// Bpipe pipeline to detect pathogenic STR expansions from whole genome sequencing data
// Takes a mapped cram as input and extracts relevant reads for remapping to STR decoys

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

hg38_bed = "~/git/STRetch-paper/reference-data/hg38.simpleRepeat_period1-6.bed"
hg38_genome = "~/git/STRetch-paper/reference-data/hg38.genome"
hg38_ref = "/humgen/atgu1/fs03/shared_resources/genomes/Homo_sapiens_assembly38.fasta"

@produce("hg38.simpleRepeat_period1-6.slop.bed")
str_targets_hg38 = {

    doc "Create bed file of region likely to contain STR reads and their mates"

    SLOP=800

        exec """
            $bedtools slop -b $SLOP -i $hg38_bed -g $hg38_genome | $bedtools merge > $output.bed
        """
}

extract_reads_region_cram = {

    doc "Extract reads from bam region + unaligned"

    def fastaname = get_fname(REF)

    produce(branch.sample + '_L001_R1.fastq.gz', branch.sample + '_L001_R2.fastq.gz') {
        exec """
            cat <( $samtools view -hu -T $hg38_ref -L $input.bed $input.cram ) 
                <( $samtools view -u -f 4 -T $hg38_ref $input.cram ) | 
            $samtools collate -Ou -n 128 - $output.prefix | 
            $bedtools bamtofastq -i - -fq >(gzip -c > $output1.gz) -fq2 >(gzip -c > $output2.gz)
        """
    }
}

run {
    str_targets_hg38 +
    '%.cram' * [
        set_sample_info +
        extract_reads_region_cram +
        align_bwa + index_bam +
        median_cov_target +
        STR_coverage +
        STR_locus_counts 
    ] +
    estimate_size
}
