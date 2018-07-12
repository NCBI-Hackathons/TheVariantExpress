#!/usr/bin/env nextflow

params.reads = "/projects/b1059/data/fastq/WI/rna/miska_rna_seq/*/*/*_L00{1,2}_R1_001.fastq.gz"
params.ref_gemone = ""
params.transcriptome = "$baseDir/test_data/c.elegans.cdna.ncrna.fa"
params.output = "results"
params.multiqc = "$baseDir/multiqc"
params.fragment_len = '250'
params.fragment_sd = '50'
params.bootstrap = '100'
params.experiment = "$baseDir/experiment_info.txt"
params.email = ""

File fq_file = new File(params.fqs)

log.info """\
              SEPSIS   P I P E L I N E  
         ===================================
         ref_gemone   : ${ params.ref_gemone }
         transcriptome: ${ params.transcriptome }
         fqs          : ${ params.fqs }
         output       : ${ params.output }
         fragment_len : ${ params.fragment_len } 
         fragment_sd  : ${ params.fragment_sd }
         bootstrap    : ${ params.bootstrap }
         experiment   : ${ params.experiment }
         email        : ${ params.email }

         """
         .stripIndent()


transcriptome_file = file(params.transcriptome)
multiqc_file = file(params.multiqc)
exp_file = file(params.experiment)
/*
 * Make sure files exist
 */

if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"

if( !exp_file.exists() ) exit 1, "Missing Experiment parameters file: ${exp_file}"



Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_1_ch; read_2_ch; read_3_ch; read_4_ch }

process qc_index {

    tag "$transcriptome_file.simpleName"

    input:
        file transcriptome from transcriptome_file

    output:
        file 'index' into index_ch

        """
        salmon index -t $transcriptome -i index
        """
}

process kal_index {

    input:
        file transcriptome_file

    output:
        file "transcriptome.index" into transcriptome_index

    script:
        //
        // Kallisto mapper index
        //
        """
        kallisto index -i transcriptome.index ${transcriptome_file}
        """
}

process kal_mapping {

	publishDir "${params.output}/quantification/kallisto", mode: "copy", overwrite: 'true'

    tag "SM"

    input:
        file index from transcriptome_index
        set SM, file(fq) from read_1_ch

    output:
        file "kallisto_${SM}" into kallisto_out_dirs

    script:
    //
    // Kallisto tools mapper
    //
    def single = fq instanceof Path
    if( !single ){
        """
        mkdir kallisto_${SM}
        kallisto quant --bootstrap ${params.bootstrap} -i ${index} -o kallisto_${SM} ${fq[0]} ${fq[1]}
        """
    }
    else {
        """
        mkdir kallisto_${SM}
        kallisto quant --single -l ${params.fragment_len} -s ${params.fragment_sd} --bootstrap ${params.bootstrap} -i ${index} -o kallisto_${SM} ${fq}
        """
    }
}

process salmon_quant {

    tag "${ SM }"

    input:
        file index from index_ch
        set SM, file(fq) from read_2_ch

    output:
        file(SM) into quant_ch

    script:
    def single = fq instanceof Path
    if( !single ){
        """
           salmon quant -l IU -i index -1 ${fq[0]} -2 ${fq[1]} -o ${SM}
        """
    }
    else {
        """
           salmon quant -i index -l U -r ${fq} -o ${SM}
        """
    }
}

process fastqc {

    tag "${ SM }"

    input:
        set SM, file(fq) from read_3_ch

    output:
        file("${SM}_log") into fastqc_ch

    script:
        """
        mkdir -p ${SM}_log
        fastqc -o ${SM}_log -f fastq -q ${fq}
        """
}

process multiqc {

	publishDir "${params.output}/report", mode: "copy", overwrite: 'true'

    input:
        file('*') from quant_ch.mix(fastqc_ch).collect()
        file(config) from multiqc_file

    output:
        file('multiqc_report.html')

    script:
        """
        cp $config/* .
        echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
        multiqc .
        """
}

process sleuth {

	publishDir "${params.output}/quantification", mode: "copy", overwrite: 'true'

    input:
        file 'kallisto/*' from kallisto_out_dirs.collect()   
        file exp_file

    output: 
        file 'sleuth_object.so'
        file 'sleuth_significant.csv' into exp_phe
        file 'gene_table.csv'

    script:
        //
        // Setup sleuth R dependancies and environment
        //
     
        """
        sleuth.R kallisto ${exp_file}
        """
}

workflow.onComplete {
    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
    println summary
    def outlog = new File("${params.output}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }
    // mail summary
    if (params.email) {
        ['mail', '-s', 'SEmRNA-seq-nf', params.email].execute() << summary
    }
}


