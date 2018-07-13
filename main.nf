#!/usr/bin/env nextflow

/*
========================================================================================
                         NCBI-Hackathons/TheVariantExpress
========================================================================================
 NCBI-Hackathons/TheVariantExpress Analysis Pipeline. Started 2018-06-21.
 https://github.com/NCBI-Hackathons/TheVariantExpress
 #### Authors
 Cong Chen < >
 Matthew Dapas < >
 Joseph Subida < >
 Octavious Talbot < >
 Chad Travis <   >
 Ye Wang <ye.wang@northwestern.edu>

----------------------------------------------------------------------------------------
*/

params.reads = "${baseDir}/test_data/fastq/*_{1,2}.fastq.gz"
params.ref_gemone = "${baseDir}/test_data/human_genome/subsample_refgenome.fa"
params.transcriptome = "${baseDir}/test_data/human_transcriptome/hg_transcriptome.fa"
params.output = "results"
params.multiqc = "${baseDir}/multiqc"
params.fragment_len = '250'
params.fragment_sd = '50'
params.bootstrap = '100'
params.experiment = "${baseDir}/info.txt"
params.email = ""

/*
Params for clinvar
 */
params.term = "sepsis"
params.variant_type = "NULL"
params.significance = 1
params.risk_factor = "NULL"

term = Channel.from(params.term)
variant_type = Channel.from(params.variant_type)
significance = Channel.from(params.significance)
risk_factor = Channel.from(params.risk_factor)


log.info """\
         ===================================
                TheVariantExpress
         ===================================
         ref_gemone   : ${ params.ref_gemone }
         transcriptome: ${ params.transcriptome }
         reads        : ${ params.reads }
         output       : ${ params.output }
         fragment_len : ${ params.fragment_len } 
         fragment_sd  : ${ params.fragment_sd }
         bootstrap    : ${ params.bootstrap }
         experiment   : ${ params.experiment }
         email        : ${ params.email }

         """
         .stripIndent()

genome_file = file(params.ref_gemone)
transcriptome_file = file(params.transcriptome)
multiqc_file = file(params.multiqc)
exp_file = file(params.experiment)
/*
 * Make sure files exist
 */

if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"

if( !exp_file.exists() ) exit 1, "Missing Experiment parameters file: ${exp_file}"

if( !genome_file.exists() ) exit 1, "MIssing reference genome file: ${genome_file}"



Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read_1_ch; read_2_ch; read_3_ch; read_4_ch; read_5_ch }


/*
=========================
   Pull clinvar data
=========================
 */

process get_clinvar_db {

    publishDir "${params.output}/clinvar", mode: "copy", overwrite: 'true'

    echo true

    input:
    val term
    val variant_type
    val significance
    val risk_factor 

    output:
    file "Clinvar_variables.tsv"

    script:
    """
    clinvar_db.R ${term} ${variant_type} ${significance} ${risk_factor}
    """
}

/*
====================
   Quantification
====================
 */

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

    tag "${SM}"

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

    echo true

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


/*
====================================
  Call variants from RNA-seq data
====================================
 */

process index_1 {

    publishDir "${params.output}/call/index_file_1", mode: "copy"

    tag "$genome_file.simpleName"

    input:
        file genome_file


    output:
        file 'index_file_1/*' into align_1_ch 

        """
        mkdir index_file_1

        STAR \\
        --runMode genomeGenerate \\
        --genomeDir index_file_1 \\
        --genomeFastaFiles ${genome_file}  \\
        --runThreadN 16
        """
}

process alignment_1 {

    publishDir "${params.output}/call/align_1", mode: "copy"

    tag "${SM}"

    input:
    set SM, file(fq) from read_4_ch
    val path from align_1_ch

    output:
    file "align_1/*" into index_2_ch


    """
    mkdir align_1
    STAR --genomeDir ${baseDir}/${params.output}call/index_file_1 --readFilesIn ${fq[0]} ${fq[1]} --runThreadN 8
    """
}

process index_2 {

    publishDir "${params.output}/call", mode: "copy"

    tag "$genome_file.simpleName"

    input:
        file genome_file
        val path from index_2_ch


    output:
        file 'index_file_2/*' into align_2_ch 

        """
        mkdir index_file_2

        STAR --runMode genomeGenerate --genomeDir index_file_2 --genomeFastaFiles ${genome_file} --sjdbFileChrStartEnd ${baseDir}/${params.output}call/align_1/SJ.out.tab --sjdbOverhang 75 --runThreadN 8
        """
}

process alignment_2 {

    publishDir "${params.output}/call/align_2", mode: "copy"

    tag "${SM}"

    input:
    set SM, file(fq) from read_5_ch
    val path from align_1_ch

    output:
    file "align_1/*" into next


    """
    mkdir align_2
    STAR --genomeDir ${baseDir}/${params.output}call/index_file_2 --readFilesIn ${fq[0]} ${fq[1]} --runThreadN 8
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


