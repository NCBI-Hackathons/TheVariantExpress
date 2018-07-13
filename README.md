![Logo](images/Logo.png)

# List the genetic variants associated with diseases,  all aboard the Variant Express!

Integrating genetic variation data with gene expression data is a powerful approach for discovering causal variants and uncovering relevant biological pathways. However, such integration typically requires multiple modes of data. The Variant Express is a new pipeline that performs association testing between gene expression estimates and genetic variants identified from RNA-Seq data alone.

The Variant Express simultaneously quantifies gene expression and identifies coding variants from RNA-Seq data and then tests for association between variants and differentially expressed genes. The Variant Express can query ClinVar to limit association testing to pathogenic variants, or can accept a user-supplied list of custom variants.

# How Does The Variant Express Do It??

Overview Diagram

![Pipeline](images/pipeline.png)

# How to use The Variant Express

### Dependencies

- R >= 3.3
- Python >= 3.6.4
- zlib >= 1.2.8

##### R Packages
- optparse
- MatrixEQTL
- httr
- jsonlite
- dplyr
- rentrez

See all [dependencies](https://github.com/NCBI-Hackathons/TheVariantExpress/blob/master/DEPENDENCIES)

### Installation

1. `git clone git@github.com:NCBI-Hackathons/TheVariantExpress.git && cd TheVariantExpress`
2. `git clone git@github.com:AndersenLab/SEmRNA-seq-nf.git`
3. `wget -qO- https://get.nextflow.io | bash`

### Running

    ./nextflow run main.nf -profile local -resume --reads "<directory_path>" --email "<your_email>" --output "output" --transcriptome "<transcriptome_file_location>"

##### Parameters

- `-profile`: The server environment; can be `local` or `quest`
- `--reads`: Path to directory of FASTQ files (Single and/or pair reads)
- `--email`: Email address to send notifications and results
- `--output`: Name of output directory
- `--transcriptome`: Transcriptome used for reads mapping


# Software Workflow Diagram

![UpdatedPipeline](https://docs.google.com/drawings/d/e/2PACX-1vRj84kE1cPLvzOBnkFm1tWz4ZjQWhGTybDpKjc9rBf2huzqlTTA3ViRTK6sJX6qW4ra-3TqnIGJPmKk/pub?w=960&h=720)


# Example Output

![img](images/Snp_11_vs_Gene_06.png)

Hackathon team: Cong Chen, Matthew Dapas, Joseph Subida, Octavious Talbot, Chad Travis, Ye Wang

