process Manta {
    container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'

    input:
    tuple path(input), path(index)
    path(fasta)
    path(fai)

    script:
    def input_files = input.collect{"--bam ${it}"}.join(' ')
    """
    configManta.py \\
        ${input_files} \\
        --reference $fasta \\
        --runDir manta

    python manta/runWorkflow.py -m local -j $task.cpus

    mv manta/results/variants/candidateSmallIndels.vcf.gz \\
        ${prefix}.candidate_small_indels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \\
        ${prefix}.candidate_small_indels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz \\
        ${prefix}.candidate_sv.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi \\
        ${prefix}.candidate_sv.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz \\
        ${prefix}.diploid_sv.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi \\
        ${prefix}.diploid_sv.vcf.gz.tbi
    """
}

params.input = "${baseDir}/data/samplesheet.aligned.csv"

workflow {
    fasta = Channel.fromPath("s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta")
    fai = Channel.fromPath("s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai")

    Channel.fromPath(params.input)
    | splitCsv(header: true)
    | map { row -> [file(row.bam), file(row.bai)] }
    | toSortedList
    | transpose
    | map { [it] }
    | collect
    | set { bams }

    Manta(bams, fasta, fai)
}

