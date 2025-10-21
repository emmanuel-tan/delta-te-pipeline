    nextflow.enable.dsl=2

    include { TRIMMOMATIC as TRIMMOMATIC_RNASEQ }               from './modules/nf-core/trimmomatic/main'                                                                                                                                                        
    include { BOWTIE2_BUILD }                                   from './modules/nf-core/bowtie2/build/main'                                                                                                                                                    
    include { BOWTIE2_ALIGN as BOWTIE2_REMOVE_CONTAMINANTS }    from './modules/nf-core/bowtie2/align/main'                                                                                                                                                    
    include { FASTQC as FASTQC_RIBO_RAW }                       from './modules/nf-core/fastqc/main'                                                                                                                                                    
    include { FASTQC as FASTQC_RIBO_TRIM }                      from './modules/nf-core/fastqc/main'                                                                                                                                                    
    include { FASTQC as FASTQC_RIBO_FILTER }                    from './modules/nf-core/fastqc/main'                                                                                                                                                    
    include { MULTIQC as MULTIQC_RIBO }                         from './modules/nf-core/multiqc/main'

    include { TRIMMOMATIC as TRIMMOMATIC_RIBOSEQ }              from './modules/nf-core/trimmomatic/main'                                                                                                                                                        
    include { FASTQC as FASTQC_RNA_RAW }                        from './modules/nf-core/fastqc/main'                                                                                                                                                    
    include { FASTQC as FASTQC_RNA_TRIM }                       from './modules/nf-core/fastqc/main'                                                                                                                                                    
    include { FASTQC as FASTQC_RNA_FILTER }                     from './modules/nf-core/fastqc/main'                                                                                                                                                    
    include { MULTIQC as MULTIQC_RNA }                          from './modules/nf-core/multiqc/main'                                                                                                                                                    

    include { SALMON_INDEX }                                    from './modules/nf-core/salmon/index/main'
    include { SALMON_QUANT as SALMON_QUANT_RIBO }               from './modules/nf-core/salmon/quant/main'
    include { SALMON_QUANT as SALMON_QUANT_RNA }                from './modules/nf-core/salmon/quant/main'

    include { MULTIQC }                                         from './modules/nf-core/multiqc/main'

    include { TXIMETA_TXIMPORT as TXIMMPORT_RIBO }              from './modules/nf-core/tximeta/tximport/main'
    include { TXIMETA_TXIMPORT as TXIMMPORT_RNA }               from './modules/nf-core/tximeta/tximport/main'

    include { CUSTOM_TX2GENE }                                  from './modules/nf-core/custom/tx2gene/main'
    include { CUSTOM_TX2GENE as CUSTOM_TX2GENE_RIBO }           from './modules/nf-core/custom/tx2gene/main'
    include { CUSTOM_TX2GENE as CUSTOM_TX2GENE_RNA }            from './modules/nf-core/custom/tx2gene/main'

    include { FORMAT_COUNT_MATRIX as FORMAT_RNA_COUNTS }        from './modules/local/formatting.nf'
    include { FORMAT_COUNT_MATRIX as FORMAT_RIBO_COUNTS }       from './modules/local/formatting.nf'

    include { SAMPLE_FASTQ as SAMPLE_FASTQ_RNA } from './modules/local/sample.nf'
    include { SAMPLE_FASTQ as SAMPLE_FASTQ_RIBO} from './modules/local/sample.nf'

    include { DELTATE } from './modules/local/deltate.nf'

    workflow {
        
        riboseq_ch = Channel
            .fromPath(params.riboseq_samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def meta = [ id: row.sample_id as String, single_end: true ]
                tuple(meta, file(row.fastq_path))
        }

        rrna_ref_ch = Channel
            .fromPath(params.abundantReference)
            .map { fasta ->
            def ref_id = fasta.baseName.replaceFirst(/\.fa(sta)?(\.gz)?$/, '')
            tuple( [ id: ref_id, single_end: true ], fasta )
        }
        
        save_unaligned_ch = Channel.value(true)
        sort_bam_ch       = Channel.value(false)

        // sampled_ribo_ch = SAMPLE_FASTQ_RIBO(riboseq_ch) 

        BOWTIE2_BUILD( rrna_ref_ch )
        // TRIMMOMATIC_RIBOSEQ(sampled_ribo_ch)
        TRIMMOMATIC_RIBOSEQ(riboseq_ch)

        ch_bowtie2_index_val = BOWTIE2_BUILD.out.index
            .groupTuple()
            .map { meta2, files -> tuple(meta2, files) }
            .first()                                
        ch_rrna_fasta_val = rrna_ref_ch.first()   

        BOWTIE2_REMOVE_CONTAMINANTS(
            TRIMMOMATIC_RIBOSEQ.out.trimmed_reads,
            ch_bowtie2_index_val,
            ch_rrna_fasta_val,
            save_unaligned_ch,
            sort_bam_ch
        )

        // FASTQC_RIBO_RAW     (   sampled_ribo_ch  )
        FASTQC_RIBO_RAW     (   riboseq_ch  )
        FASTQC_RIBO_TRIM    (   TRIMMOMATIC_RIBOSEQ.out.trimmed_reads   )
        FASTQC_RIBO_FILTER  (   BOWTIE2_REMOVE_CONTAMINANTS.out.fastq   )

        rnaseq_ch = Channel
            .fromPath(params.rnaseq_samplesheet)
            .splitCsv(header:true)
            .filter { row -> row?.fastq_path && row.fastq_path.trim() }
            .map { row ->
                def meta = [ id: row.sample_id as String, single_end: true ]
                tuple(meta, [ file(row.fastq_path) ])
        }

        // sampled_rna_ch = SAMPLE_FASTQ_RNA(rnaseq_ch)

        // rnaseq_ch = Channel
        //     .fromPath(params.rnaseq_samplesheet)
        //     .splitCsv(header:true)
        //     .filter { row -> row?.fastq_r1_path && row?.fastq_r2_path }
        //     .map { row ->
        //         def meta = [ id: row.sample_id as String, single_end: false ]
        //         tuple(meta, [ file(row.fastq_r1_path), file(row.fastq_r2_path) ])
        // }

        TRIMMOMATIC_RNASEQ( rnaseq_ch )

        // FASTQC_RNA_RAW( sampled_rna_ch )
        FASTQC_RNA_RAW( rnaseq_ch )
        FASTQC_RNA_TRIM( TRIMMOMATIC_RNASEQ.out.trimmed_reads )

        // tx_fa_ch            = Channel.fromPath(params.humanReferenceTx)            
        // genome_fa_ch        = Channel.fromPath(params.humanReferenceGenome)        
        // gtf_ch              = Channel.fromPath(params.humanReferenceAnnot).first()
        // alignment_mode_ch   = Channel.value(false)  
        // libtype_ch          = Channel.value( 'A' )

        // rpf_reads = BOWTIE2_REMOVE_CONTAMINANTS.out.fastq
        // rna_reads = TRIMMOMATIC_RNASEQ.out.trimmed_reads

        // if (params.salmon_index) {
        //     idxPath = Channel.fromPath(params.salmon_index)
        //     idx_ch = Channel.fromPath(params.salmon_index).map { it -> it }   
        //     idx_val = idx_ch.first()
        // }
        // else {
        //     SALMON_INDEX(genome_fa_ch, tx_fa_ch)
        //     idx_val = SALMON_INDEX.out.index.first()
        // }

        // gtf_val  = gtf_ch.first()
        // txfa_val = tx_fa_ch.first()

        // SALMON_QUANT_RIBO(
        //     rpf_reads,
        //     idx_val,
        //     gtf_val,
        //     txfa_val,
        //     alignment_mode_ch,
        //     libtype_ch
        // )

        // SALMON_QUANT_RNA(
        //     rna_reads,
        //     idx_val,
        //     gtf_val,
        //     txfa_val,
        //     alignment_mode_ch,
        //     libtype_ch
        // )

        // ribo_quants_list = SALMON_QUANT_RIBO.out.results       
        //     .map { meta, qdir -> qdir }                   
        //     .collect()                                    
        //     .map { lst -> tuple([ id: 'ribo' ], lst) }
        // rna_quants_list = SALMON_QUANT_RNA.out.results
        //     .map { meta, qdir -> qdir }
        //     .collect()
        //     .map { lst -> tuple([ id: 'rna' ], lst) }
        

        // id_choice_ch  = Channel.value(params.tx2gene_id ?: 'gene_id')
        // mode_ch       = Channel.value(params.tx_id_version ?: 'gene_name')
        // quant_type_ch = Channel.value('salmon')

        // gtf_for_tx2gene = gtf_ch.map { f -> tuple([ id:'annot' ], f) }

        // CUSTOM_TX2GENE_RIBO(
        //     gtf_for_tx2gene,
        //     ribo_quants_list,
        //     quant_type_ch,
        //     id_choice_ch,
        //     mode_ch
        // )

        // TXIMMPORT_RIBO(
        //     ribo_quants_list, 
        //     CUSTOM_TX2GENE_RIBO.out.tx2gene, 
        //     quant_type_ch
        // )
        
        // CUSTOM_TX2GENE_RNA(
        //     gtf_for_tx2gene,
        //     rna_quants_list,
        //     quant_type_ch,
        //     id_choice_ch,
        //     mode_ch
        // )

        // TXIMMPORT_RNA(
        //     rna_quants_list, 
        //     CUSTOM_TX2GENE_RNA.out.tx2gene, 
        //     quant_type_ch
        // )

        // rna_counts = TXIMMPORT_RNA.out.counts_gene
        // ribo_counts = TXIMMPORT_RIBO.out.counts_gene

        // FORMAT_RNA_COUNTS( rna_counts )
        // FORMAT_RIBO_COUNTS( ribo_counts )

        // rna_mat_ch  = FORMAT_RNA_COUNTS.out.tsv .map { meta, f -> tuple([id: 'deltate'], f) }
        // ribo_mat_ch = FORMAT_RIBO_COUNTS.out.tsv.map { meta, f -> tuple([id: 'deltate'], f) }

        // sample_info_ch = Channel
        //     .fromPath(params.deltate_samplesheet, checkIfExists: true)
        //     .map { f -> tuple([id: 'deltate'], f) }
        //     .first()
        // batch_flag_ch  = 0

        // DELTATE(
        //     rna_mat_ch,
        //     ribo_mat_ch,
        //     sample_info_ch,
        //     batch_flag_ch
        // )

        def onlyPath = { ch -> ch.map { meta, f -> f } }

        ch_multiqc_files = Channel
            .empty()
            // fastqc artifacts
            .mix( onlyPath(FASTQC_RIBO_RAW.out.zip) )
            .mix( onlyPath(FASTQC_RIBO_RAW.out.html) )
            .mix( onlyPath(FASTQC_RIBO_TRIM.out.zip) )
            .mix( onlyPath(FASTQC_RIBO_TRIM.out.html) )
            .mix( onlyPath(FASTQC_RIBO_FILTER.out.zip) )
            .mix( onlyPath(FASTQC_RIBO_FILTER.out.html) )
            // trimmomatic logs
            .mix( onlyPath(TRIMMOMATIC_RIBOSEQ.out.summary) )
            .mix( onlyPath(TRIMMOMATIC_RIBOSEQ.out.trim_log) )
            // bowtie2 logs
            .mix( onlyPath(BOWTIE2_REMOVE_CONTAMINANTS.out.log) )
            .mix( onlyPath(FASTQC_RNA_RAW.out.zip),     onlyPath(FASTQC_RNA_RAW.out.html) )
            .mix( onlyPath(FASTQC_RNA_TRIM.out.zip),    onlyPath(FASTQC_RNA_TRIM.out.html) )
            // .mix( onlyPath(TRIMMOMATIC_RNASEQ.out.summary), onlyPath(TRIMMOMATIC_RNASEQ.out.trim_log) )
            // .mix( onlyPath(SALMON_QUANT_RIBO.out.json_info) )
            // .mix( onlyPath(SALMON_QUANT_RIBO.out.lib_format_counts) )
            // .mix( onlyPath(SALMON_QUANT_RIBO.out.results) )
            // .mix( onlyPath(SALMON_QUANT_RNA.out.json_info) )
            // .mix( onlyPath(SALMON_QUANT_RNA.out.lib_format_counts) )   
            // .mix( onlyPath(SALMON_QUANT_RNA.out.results) )   
            .collect()
            .ifEmpty([])

        def none = Channel.value([])

        MULTIQC(
        ch_multiqc_files,
        none,
        none,
        none,
        none,
        none
        )
    }