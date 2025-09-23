#!/usr/bin/env nextflow

process FORMAT_COUNT_MATRIX {
    tag "${meta.id ?: file(infile).name}"
    label 'process_single'

    input:
    tuple val(meta), path(infile)

    output:
    tuple val(meta), path("*.formatted.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def base = infile.getBaseName().replaceAll(/\.tsv$|\.txt$|\.csv$|\.gz$/,'')
    def catcmd = infile.name.endsWith('.gz') ? "gzip -cd ${infile}" : "cat ${infile}"

    """
    set -euo pipefail

    ${catcmd} | awk -F '\\t' -v OFS='\\t' '
        NR==1 {
            # Header: drop col1 (gene_id header) and col2 (gene_name header)
            # Print only sample IDs (cols 3..NF)
            out = (NF>=3) ? \$3 : ""
            for (i=4; i<=NF; i++) out = out OFS \$i
            print out
            next
        }
        {
            # Data: keep gene_id (col1), drop gene_name (col2), keep samples (3..NF) rounded
            out = \$1
            for (i=3; i<=NF; i++) out = out OFS int(\$i+0.5)
            print out
        }
    ' > "${base}.formatted.tsv"
    """
}
