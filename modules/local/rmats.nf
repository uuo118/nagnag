// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process RMATS {
    label 'process_high'


    conda (params.enable_conda ? "bioconda::rmats=4.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rmats:4.1.2--py37haf75f70_1':
        'quay.io/biocontainers/rmats:4.1.2--py37haf75f70_1' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    val bams
    path gtf
    val readlength


    output:
    path "*MATS.JCEC.txt"         , emit: jcec
    path "*MATS.JC.txt"           , emit: jc
    path "summary.txt"            , emit: summary
    path "comparison.txt"         , emit: comparison
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    def args = task.ext.args ?: ''
    def lib = bams[0][0].single_end ? "single" : "paired"
    bam1 = file("${workDir}/bam1.txt")
    bam1.append(bams[0][1].join(","))
    // ${workDir}/rmats/${bams[0][0].id}vs${bams[1][0].id}
    if( bams.size() > 1){
        bam2 = file("${workDir}/bam2.txt")
        bam2.append(bams[1][1].join(","))

        """
        mv ${workDir}/bam1.txt \${PWD}
        mv ${workDir}/bam2.txt \${PWD}
        echo "${bams[0][0].id}_vs_${bams[1][0].id}" > comparison.txt
        rmats.py \\
            $args \\
            --b1 \${PWD}/bam1.txt\\
            --b2 \${PWD}/bam2.txt \\
            --gtf $gtf \\
            --nthread $task.cpus \\
            --od \${PWD}/ \\
            --tmp \${PWD}/tmp \\
            -t $lib \\
            --readLength $readlength

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rmats: \$(echo \$(rmats.py --version 2>&1) | sed 's/^.*rmats //; s/Using.*\$//' ))
        END_VERSIONS
        """
        } else {
        """
        echo "${bams[0][0].id}_only" > comparison.txt
        mv ${workDir}/bam1.txt \${PWD}
        rmats.py \\
            $args \\
            --b1 \${PWD}/bam1.txt\\
            --nthread $task.cpus \\
            --gtf $gtf \\
            --od \${PWD}/ \\
            --tmp \${PWD}/tmp \\
            -t $lib \\
            --readLength $readlength \\
            --statoff

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rmats: \$(echo \$(rmats.py --version 2>&1) | sed 's/^.*rmats //; s/Using.*\$//' ))
        END_VERSIONS
        """
        }
}


