process FASTQC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // --- START: patched preamble ---
    def args    = task.ext.args   ?: ''
    def suffix  = task.ext.suffix ?: null
    // If prefix is given, it wins. Else use meta.id or meta.id.suffix
    def prefix  = task.ext.prefix ?: (suffix ? "${meta.id}.${suffix}" : "${meta.id}")
    // --- END: patched preamble ---

    // Preserve compound extensions like .fastq.gz rather than only ".gz"
    def mkNewName = { String base, Path p ->
        def name = p.name
        def dot  = name.indexOf('.')
        def ext  = dot >= 0 ? name.substring(dot) : ''   // includes leading dot(s)
        return "${base}${ext}"
    }

    // Build list of (old -> new) names; for PE, add _1/_2 to base
    def old_new_pairs =
        (reads instanceof Path || reads.size() == 1)
            ? [[ reads, mkNewName(prefix, (reads instanceof Path ? reads : reads[0])) ]]
            : reads.withIndex().collect { p, i -> [ p, mkNewName("${prefix}_${i+1}", p) ] }

    def rename_to     = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ _old, newn -> newn }.join(' ')

    // FastQC memory budgeting as in nf-core
    def memory_in_mb  = task.memory ? task.memory.toUnit('MB') / task.cpus : null
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)

    """
    printf "%s %s\\n" ${rename_to} | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s "\$old_name" "\$new_name"
    done

    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        --memory ${fastqc_memory} \\
        ${renamed_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """

    stub:
    // --- START: patched stub to mirror suffix/prefix behaviour ---
    def suffix  = task.ext.suffix ?: null
    def prefix  = task.ext.prefix ?: (suffix ? "${meta.id}.${suffix}" : "${meta.id}")
    // --- END: patched stub ---
    """
    touch ${prefix}_fastqc.html
    touch ${prefix}_fastqc.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
