process primer_trimming {
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("trimmed.fastq.gz")

    """
    PRIMER_F=${params.primer_f}
    PRIMER_R=${params.primer_r}

    cutadapt \
        --discard-untrimmed \
        -O \${#PRIMER_F} \
        -g \${PRIMER_F} \
        -e 0.2 \
        --revcomp \
        ${reads} 2> ${meta.id}_trimmed.log | \
        cutadapt \
            --discard-untrimmed \
            -O \${#PRIMER_R} \
            -a \${PRIMER_R} \
            -e 0.2 \
            - 2>> ${meta.id}_trimmed.log |
        gzip > trimmed.fastq.gz
    """
}

process quality_filtering {
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("filtered.fastq.gz")

    """
    vsearch --fastq_filter ${reads} \
        --fastq_maxee_rate 0.2 \
        --fastq_qmax 90 \
        --fastq_minlen ${params.min_length} \
        --fastq_maxlen ${params.max_length} \
        --fastqout filtered.fastq.gz \
        --fasta_width 0
    """
}

process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastqc=0.12.1"
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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.${entry.extension}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
        $args \\
        --threads $task.cpus \\
        $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}

process PREP_FOR_CLUST {
    
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

    cpus params.clusteringcpus

    input:
    path(reads)
    output:
    path("otutab.txt"), emit: otutab
    path("consensus.fasta"), emit: otuseq

    """
    cat ${reads} > concat_reads.fasta.gz

    vsearch \
        --cluster_fast concat_reads.fasta.gz \
        --sizein \
        --consout consensus.fasta \
        --otutabout otutab.txt \
        --threads ${params.clusteringcpus} \
        --id 0.97 \
        --iddef 0

    rm concat_reads.fasta.gz

    """
}

process GET_REF {
    output:
    path("refdb.rds")

    """
    #!/usr/bin/env Rscript

    if ("${params.refdb}" == "silva"){
        system("wget http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.RData")

        load("SILVA_SSU_r138_2019.RData", ex <- new.env())
        tmp <- get(names(ex)[1],envir=ex)
        saveRDS(tmp,file="refdb.rds")
    } else if ("${params.refdb}" == "pr2") {
        system("wget https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU.decipher.trained.rds")
        file.rename("pr2_version_5.0.0_SSU.decipher.trained.rds", "refdb.rds")
    }
    
    """
}

process ASSIGN_IDTAXA {

    memory 20.GB

    input:
    path(fasta_seq)
    path(refdb)

    output:
    path("taxonomy.tsv")

    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(Biostrings)

    ##########################################
    # import sequences
    ##########################################

    sequences <- readDNAStringSet("${fasta_seq}", format = "fasta")
    names(sequences) <- str_remove_all(names(sequences), "^centroid=|;.+\$")

    ##########################################
    # assign with idtaxa
    ##########################################

    trainingSet <- readRDS("${refdb}")

    idtaxa_res <- DECIPHER::IdTaxa(
        sequences,
        trainingSet,
        strand = "top",
        threshold = 50,
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

    res <- tibble(
        sequence = names(taxonomy),
        taxonomy = str_remove(taxonomy, "^Root;"),
        confidence = str_remove(confidence, "^[^;]+;"),
        )

    write_tsv(res, "taxonomy.tsv")
    """
}

process FINAL_TABLE {
    
    input:
    path(otutab)
    path(taxonomy)

    output:
    path("final_table.tsv")

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    otutab <- read_tsv("${otutab}")
    names(otutab)[1] <- "sequence"

    taxonomy <- read_tsv("${taxonomy}")

    right_join(taxonomy, otutab) %>%
        write_tsv("final_table.tsv")
    """

}

workflow {

    Channel.fromPath(params.manifest)
        .splitCsv(sep: " ", header: ['sampleId', 'raw_fastq'])
        .map { row -> tuple([id:row.sampleId], file(row.raw_fastq)) }
        .set { manifest }

    primer_trimming(manifest)
        .set { trimmed }

    FASTQC(trimmed)

    quality_filtering(trimmed)
        .set { qual_filt }

    PREP_FOR_CLUST(qual_filt)
        .set{ dereplicated }

    dereplicated
        .map { info, reads -> reads}
        .collect()
        .set { all_reads }

    CLUSTERING(all_reads)
        .set{ cluster_res }

    GET_REF().set{ refbd }

    cluster_res.otuseq
        .splitFasta( by: 1000, file: true )
        .set{ splitted_fasta }

    ASSIGN_IDTAXA(splitted_fasta, refbd)
        .collectFile(keepHeader: true, skip: 1)
        .set{ taxo_res }

    FINAL_TABLE(cluster_res.otutab, taxo_res)

}