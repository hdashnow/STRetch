// Bpipe pipeline to detect pathogenic STR expansions from whole genome sequencing data
// Takes a mapped bam as input and extracts relevant reads for remapping to STR decoys

// Load system configuration and other settings
load 'pipeline_config.groovy'

// Load Bpipe pipeline stages
load 'pipeline_stages.groovy'

@filter("STRregions")
extract_region_unmapped = {
    exec """
            set -o pipefail

            cat <( $samtools view -hu -L $input.bed $input.bam )
                <( $samtools view -u -f 4 $input.bam ) > $output.bam
    """
}

@filter("rnamesorted")
rname_sort = {
    exec """
        $samtools collate -O -n 128 $input.bam $output.prefix > $output.bam
    """
}

extract_fastq = {
    def R1temp = branch.sample + '_L001_R1.fastq'
    def R2temp = branch.sample + 'L001_R2.fastq'

    produce(branch.sample + '_L001_R1.fastq.gz', branch.sample + '_L001_R2.fastq.gz') {
        exec """
            set -o pipefail

            $bedtools bamtofastq -i $input.bam -fq $R1temp -fq2 $R2temp &&
            gzip -c $R1temp > $output1.gz &&
            gzip -c $R2temp > $output2.gz
        """
    }
}

gzip = {
    from(branch.sample + '_L001_R1.fastq', branch.sample + '_L001_R1.fastq') produce(branch.sample + '_L001_R1.fastq.gz', branch.sample + '_L001_R2.fastq.gz') {
        exec """
            
        """
    }
}

run {
    str_targets +
    '%.bam' * [
        set_sample_info +

        extract_region_unmapped +
        rname_sort +
        extract_fastq +

        align_bwa + index_bam +
        median_cov_target +
        STR_coverage +
        STR_locus_counts 
    ] +
    estimate_size
}
