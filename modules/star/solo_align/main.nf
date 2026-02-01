process STARSOLO {
    tag "$meta.id"
    label 'process_high'

    container 'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4'

    input:
    tuple val(meta), path(reads)
    path(index)
    path(whitelist)

    output:
    tuple val(meta), path('*d.out.bam')          , emit: bam
    tuple val(meta),  path('*.Solo.out')         , emit: counts
    tuple val(meta),  path('*Log.final.out')     , emit: log_final
    tuple val(meta),  path('*Log.out')           , emit: log_out
    tuple val(meta),  path('*Log.progress.out')  , emit: log_progress
    tuple val(meta),  path('*/Gene/Summary.csv') , emit: summary
    tuple val(meta), path('*.tab')               , optional:true, emit: tab
    path "versions.yml"                          , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def (forward, reverse) = reads.collate(2).transpose()
    def zcat = reads[0].getExtension() == "gz" ? "--readFilesCommand zcat": ""

    """

    if [[ $whitelist == *.gz ]]; then
        gzip -cdf $whitelist > whitelist.uncompressed.txt
    else
        cp $whitelist whitelist.uncompressed.txt
    fi

    STAR \\
        --genomeDir $index \\
        --readFilesIn ${reverse.join( "," )} ${forward.join( "," )} \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        --soloType CB_UMI_Simple \\
        --soloCBwhitelist whitelist.uncompressed.txt \\
        $zcat \\
        $args 

    if [ -d ${prefix}.Solo.out ]; then
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip {} \\;
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.Solo.out/Gene
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.Solo.out/Gene/Summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}