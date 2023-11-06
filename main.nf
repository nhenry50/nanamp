#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nanamp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nhenry50/nanamp
----------------------------------------------------------------------------------------
*/

process PRIMER_TRIMMING {

    label 'process_single'

    conda "bioconda::cutadapt=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
        'biocontainers/cutadapt:3.4--py39h38f01e4_1' }"

    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("trimmed.fastq.gz"), emit: fastqs
    path("trimming.log"), emit: log

    """
    PRIMER_F=${params.primer_f}
    PRIMER_R=${params.primer_r}

    cutadapt \
        --discard-untrimmed \
        -O \${#PRIMER_F} \
        -g \${PRIMER_F} \
        -e ${params.cutadapt_error} \
        --revcomp \
        --report=minimal \
        ${reads} 2> tmp.log | \
        cutadapt \
            --discard-untrimmed \
            -O \${#PRIMER_R} \
            -a \${PRIMER_R} \
            -e ${params.cutadapt_error} \
            --report=minimal \
            - 2>> tmp.log |
        gzip > trimmed.fastq.gz

    awk 'NR == 2 {printf "${meta.id} "\$2" "} NR == 4 {print \$7}' tmp.log > trimming.log

    """
}

process QUALITY_FILTERING {

    label 'process_single'

    conda "bioconda::vsearch=2.21.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("filtered.fastq.gz"), emit: fastqs
    path("filtering.log"), emit: log

    """
    vsearch --fastq_filter ${reads} \
        --fastq_maxee_rate ${params.maxee_rate} \
        --fastq_qmax 90 \
        --fastq_minlen ${params.min_length} \
        --fastq_maxlen ${params.max_length} \
        --fastqout filtered.fastq.gz \
        --fasta_width 0 2> tmp.log

    awk 'NR == 5 {print "${meta.id}",\$1}' tmp.log > filtering.log
    """
}

process PREP_FOR_CLUST {

    label 'process_single'

    conda "bioconda::vsearch=2.21.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"
    
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.id}_derep.fasta.gz")

    """
    vsearch --fastx_uniques ${reads} --sizeout --relabel_sha1 --fastaout - | \
        awk '/>/ {print \$0";sample=${meta.id}"}; !/>/{print}' | \
        gzip > "${meta.id}_derep.fasta.gz"
    """
}

process CLUSTERING {

    label 'process_high'

    conda "bioconda::vsearch=2.21.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    path(reads)
    output:
    path("otutab.txt"), emit: otutab
    path("consensus.fasta"), emit: otuseq
    path("clustering.log"), emit: log

    """
    cat ${reads} > concat_reads.fasta.gz

    vsearch \
        --cluster_fast concat_reads.fasta.gz \
        --sizein \
        --consout consensus.fasta \
        --otutabout otutab.txt \
        --threads $task.cpus \
        --id ${params.clusteringid} \
        --iddef ${params.clusteringiddef} 2> tmp.log

    
    awk 'BEGIN {print "sequences\tclusters"}; NR == 5 {printf \$4"\t"}; NR == 12 {print \$2}' tmp.log > clustering.log

    rm concat_reads.fasta.gz

    """
}

process GET_REF {

    label 'process_low'

    conda "r::r-tidyverse=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1' :
        'biocontainers/r-tidyverse:1.2.1' }"

    output:
    path("refdb.rds")

    """
    #!/usr/bin/env Rscript

    options(timeout = 1200)

    if ("${params.refdb}" == "silva"){

        download.file(
            "http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData",
            "SILVA_SSU_r138_2019.RData",
            quiet = TRUE
        )

        load("SILVA_SSU_r138_2019.RData", ex <- new.env())

        tmp <- get(names(ex)[1],envir=ex)

        saveRDS(tmp,file="refdb.rds")

    } else if ("${params.refdb}" == "pr2") {

        download.file(
            "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU.decipher.trained.rds",
            "refdb.rds",
            quiet = TRUE
        )

    }
    
    """
}

process ASSIGN_IDTAXA {

    label 'process_low'

    conda "bioconda::bioconductor-decipher=2.28.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-decipher:2.28.0--r43ha9d7317_0' :
        'biocontainers/bioconductor-decipher:2.28.0--r43ha9d7317_0' }"

    input:
    path(fasta_seq)
    path(refdb)

    output:
    path("taxonomy_${params.refdb}.tsv")

    """
    #!/usr/bin/env Rscript

    ##########################################
    # import sequences
    ##########################################

    sequences <- Biostrings::readDNAStringSet("${fasta_seq}", format = "fasta")
    names(sequences) <- gsub("^centroid=|;.+\$", "", names(sequences))

    ##########################################
    # assign with idtaxa
    ##########################################

    trainingSet <- readRDS("${refdb}")

    idtaxa_res <- DECIPHER::IdTaxa(
        sequences,
        trainingSet,
        strand = "top",
        threshold = ${params.idtaxa_thresh},
        minDescend = 0.9,
        processors = 1
    )

    taxonomy <- vapply(
        idtaxa_res,
        function(x) paste(x[["taxon"]], collapse = ";"),
        character(1)
    )

    confidence <- vapply(
        idtaxa_res,
        function(x) paste(round(x[["confidence"]], digits = 1), collapse = ";"),
        character(1)
    )


    ##########################################
    # assemble into one table
    ##########################################

    res <- data.frame(
        sequence = names(taxonomy),
        taxonomy = sub("^Root;", "", taxonomy),
        confidence = sub("^[^;]+;", "",confidence)
        )

    write.table(
        res,
        file = "taxonomy_${params.refdb}.tsv",
        quote = FALSE,
        sep = "\t",
        row.names = FALSE
    )

    """
}

process FINAL_TABLE {

    conda "r::r-tidyverse=1.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1' :
        'biocontainers/r-tidyverse:1.2.1' }"

    input:
    path(otutab)
    path(taxonomy)

    output:
    path("final_table_${params.refdb}.tsv")

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    otutab <- read.delim("${otutab}")
    names(otutab)[1] <- "sequence"

    taxonomy <- read.delim("${taxonomy}")

    right_join(taxonomy, otutab) %>%
        write_tsv("final_table_${params.refdb}.tsv")
    """

}

process GATHER_LOGS {

    input:
    path(trim_log)
    path(filter_log)

    output:
    path("workflow_stats.log")

    """
    (
        echo "sample in_reads with_primer filtered"
        join <(sort ${trim_log}) <(sort ${filter_log})
    ) > workflow_stats.log
    """

}

workflow {

    Channel.fromPath(params.manifest)
        .splitCsv(sep: " ", header: ['sampleId', 'raw_fastq'])
        .map { row -> tuple([id:row.sampleId], file(row.raw_fastq)) }
        .set { manifest }

    PRIMER_TRIMMING(manifest)
        .set { trimmed }

    QUALITY_FILTERING(trimmed.fastqs)
        .set { qual_filt }

    PREP_FOR_CLUST(qual_filt.fastqs)
        .set{ dereplicated }

    dereplicated
        .map { info, reads -> reads}
        .collect()
        .set { all_reads }

    CLUSTERING(all_reads)
        .set{ cluster_res }

    GET_REF().set{ refbd }

    cluster_res.otuseq
        .splitFasta( by: params.fastachunks, file: true )
        .set{ splitted_fasta }

    ASSIGN_IDTAXA(splitted_fasta, refbd)
        .collectFile(keepHeader: true, skip: 1)
        .set{ taxo_res }

    FINAL_TABLE(cluster_res.otutab, taxo_res)

    GATHER_LOGS(trimmed.log.collectFile(), qual_filt.log.collectFile())

}