nextflow.enable.dsl=2

process DELTATE {
  tag { meta.id ?: 'deltate' }
  label 'process_medium'
  container '107laserfish/deltate-bioconductor:latest'

  containerOptions "--platform=linux/amd64 -v ${task.workDir}:/work -w /work --entrypoint ''"


  input:
    tuple val(meta),  path(rna_counts)
    tuple val(meta2), path(ribo_counts)
    tuple val(meta3), path(sample_info)
    val batch_flag

  output:
    tuple val(meta), path("Results/*"), emit: results
    tuple val(meta), path("deltate_run.log"), emit: log, optional: true

  script:
  def riboArg   = "${ribo_counts.getName()}"
  def rnaArg    = "${rna_counts.getName()}"
  def sampleArg = "${sample_info.getName()}"

  """
  set -euo pipefail

  mkdir -p Results

  echo "== CWD =="; pwd
  echo "== LS work dir =="; ls -l
  echo "== Preview =="
  head -n2 "${riboArg}" || true
  head -n2 "${rnaArg}"  || true
  head -n2 "${sampleArg}" || true

  Rscript --vanilla /home/ruser/scripts/DTEG.R \\
    "/work/${riboArg}" \\
    "/work/${rnaArg}"  \\
    "/work/${sampleArg}" \\
    ${batch_flag} | tee /work/deltate_run.log
  """
}