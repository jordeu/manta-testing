process Manta {
    container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
    cpus 10
    memory '8GB'

    input:
    tuple path(input, stageAs: 'sample.*.bam'), path(index, stageAs: 'sample.*.bam.bai')
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
        example.candidate_small_indels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \\
        example.candidate_small_indels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz \\
        example.candidate_sv.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi \\
        example.candidate_sv.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz \\
        example.diploid_sv.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi \\
        example.diploid_sv.vcf.gz.tbi
    """
}

params.input = "${baseDir}/data/samplesheet.aligned.csv"

workflow {
    fasta = Channel.fromPath("s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")
    fai = Channel.fromPath("s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa.fai")

    Channel.fromPath(params.input)
    | splitCsv(header: true)
    | map { row -> [file(row.bam), file(row.bai)] }
    | toSortedList
    | transpose
    | toSortedList
    | set { bams }

    Manta(bams, fasta, fai)
}

