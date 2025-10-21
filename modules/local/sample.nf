process SAMPLE_FASTQ {
    tag "$meta.id"

    input:
      tuple val(meta), path(reads)

    output:
      tuple val(meta), path("*.sampled.fastq.gz")

    when:
      params.sample_reads && params.sample_reads.toInteger() > 0

    script:
      def nreads = params.sample_reads.toInteger()
      """
      zcat ${reads} | paste - - - - | shuf -n ${nreads} | tr '\\t' '\\n' | gzip > ${meta.id}.sampled.fastq.gz
      """
}
