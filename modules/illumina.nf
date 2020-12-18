process readTrimming {
    /**
    * Trims paired fastq using trim_galore (https://github.com/FelixKrueger/TrimGalore)
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output trimgalore_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*_val_{1,2}.fq.gz', mode: 'copy'

    cpus 2

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz")) optional true

    script:
    """
    if [[ \$(gunzip -c ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      exit 0
    else
      trim_galore --paired $forward $reverse
    fi
    """
}

process filterResidualAdapters {
    /**
    * Discard reads that contain residual adapter sequences that indicate trimming may have failed
    * @input tuple(sampleName, path(forward), path(reverse))
    * @output untrim_filter_out tuple(sampleName, path("*_val_1.fq.gz"), path("*_val_2.fq.gz"))
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: '*{1,2}_posttrim_filter.fq.gz', mode: 'copy'

    cpus 2

    input:
    tuple(sampleName, path(forward), path(reverse))

    output:
    tuple(sampleName, path("*1_posttrim_filter.fq.gz"), path("*2_posttrim_filter.fq.gz")) optional true

    script:
    """
    if [[ \$(gunzip -c ${forward} | head -n4 | wc -l) -eq 0 ]]; then
      exit 0
    else
      filter_residual_adapters.py --input_R1 $forward --input_R2 $reverse 
    fi
    """
}

process indexReference {
    /**
    * Indexes reference fasta file in the scheme repo using bwa.
    */

    tag { ref }

    input:
        path(ref)

    output:
        tuple path('ref.fa'), path('ref.fa.*')

    script:
        """
        ln -s ${ref} ref.fa
        bwa index ref.fa
        """
}

process readMapping {
    /**
    * Maps trimmed paired fastq using BWA (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to BAM, sort and index sorted BAM (http://www.htslib.org/doc/samtools.html)
    * @input 
    * @output 
    */

    tag { sampleName }

    label 'largecpu'

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.sorted{.bam,.bam.bai}", mode: 'copy'

    input:
        tuple sampleName, path(forward), path(reverse), path(ref), path("*")

    output:
        tuple(sampleName, path("${sampleName}.sorted.bam"), path("${sampleName}.sorted.bam.bai"))

    script:
        """
        bwa mem -t ${task.cpus} ${ref} ${forward} ${reverse} | \
        samtools sort -o ${sampleName}.sorted.bam
        samtools index ${sampleName}.sorted.bam
        """
}

process trimPrimerSequences {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped{.bam,.bam.bai}", mode: 'copy'
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.mapped.primertrimmed.sorted{.bam,.bam.bai}", mode: 'copy'

    input:
    tuple sampleName, path(bam), path(bam_index), path(bedfile)

    output:
    tuple sampleName, path("${sampleName}.mapped.bam"), path("${sampleName}.mapped.bam.bai"), emit: mapped
    tuple sampleName, path("${sampleName}.mapped.primertrimmed.sorted.bam"), path("${sampleName}.mapped.primertrimmed.sorted.bam.bai" ), emit: ptrim

    script:
    if (params.allowNoprimer){
        ivarCmd = "ivar trim -e"
    } else {
        ivarCmd = "ivar trim"
    }
        """
        samtools view -F4 -o ${sampleName}.mapped.bam ${bam}
        samtools index ${sampleName}.mapped.bam
        ${ivarCmd} -i ${sampleName}.mapped.bam -b ${bedfile} -m ${params.illuminaKeepLen} -q ${params.illuminaQualThreshold} -f ${params.primer_pairs_tsv} -p ivar.out
        samtools sort -o ${sampleName}.mapped.primertrimmed.sorted.bam ivar.out.bam
        samtools index ${sampleName}.mapped.primertrimmed.sorted.bam
        """
}

process callVariants {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.variants.tsv", mode: 'copy'

    input:
    tuple(sampleName, path(bam), path(bam_index), path(ref))

    output:
    tuple sampleName, path("${sampleName}.variants.tsv")

    script:
        """
        samtools faidx ${ref}
        samtools mpileup -A -d 0 --reference ${ref} -B -Q 0 ${bam} |\
        ivar variants -r ${ref} -m ${params.ivarMinDepth} -p ${sampleName}.variants -q ${params.ivarMinVariantQuality} -t ${params.ivarMinFreqThreshold}
        """
}

process makeConsensus {

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.fa", mode: 'copy'

    input:
        tuple(sampleName, path(bam), path(bam_index))

    output:
        tuple(sampleName, path("${sampleName}.primertrimmed.consensus.fa"))

    script:
        """
        samtools mpileup -aa -A -B -d ${params.mpileupDepth} -Q0 ${bam} | \
        ivar consensus -t ${params.ivarFreqThreshold} -m ${params.ivarMinDepth} \
        -n N -p ${sampleName}.primertrimmed.consensus
        """
}

process cramToFastq {
    /**
    * Converts CRAM to fastq (http://bio-bwa.sourceforge.net/)
    * Uses samtools to convert to CRAM, to FastQ (http://www.htslib.org/doc/samtools.html)
    * @input
    * @output
    */

    input:
        tuple sampleName, file(cram)

    output:
        tuple sampleName, path("${sampleName}_1.fastq.gz"), path("${sampleName}_2.fastq.gz")

    script:
        """
        samtools collate -u ${cram} -o tmp.bam
        samtools fastq -1 ${sampleName}_1.fastq.gz -2 ${sampleName}_2.fastq.gz tmp.bam
        rm tmp.bam
        """
}

process alignConsensusToReference {
    /**
    * Aligns consensus sequence against reference using mafft. Uses the --keeplength
    * flag to guarantee that all alignments remain the same length as the reference.
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.aln.fa", mode: 'copy'

    input:
        tuple(sampleName, path(consensus), path(reference))

    output:
        tuple(sampleName, path("${sampleName}.primertrimmed.consensus.aln.fa"))

    script:
        // Convert multi-line fasta to single line
        awk_string = '/^>/ {printf("\\n%s\\n", $0); next; } { printf("%s", $0); }  END { printf("\\n"); }'
        """
        mafft \
          --preservecase \
          --keeplength \
          --add \
          ${consensus} \
          ${reference} \
          > ${sampleName}.with_ref.multi_line.alignment.fa
        awk '${awk_string}' ${sampleName}.with_ref.multi_line.alignment.fa > ${sampleName}.with_ref.single_line.alignment.fa
        tail -n 2 ${sampleName}.with_ref.single_line.alignment.fa > ${sampleName}.primertrimmed.consensus.aln.fa
        """
}

process trimUTRFromAlignment {
    /**
    * Trim the aligned consensus to remove 3' and 5' UTR sequences.
    */

    tag { sampleName }

    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.primertrimmed.consensus.aln.utr_trimmed.fa", mode: 'copy'

    input:
        tuple(sampleName, path(alignment))

    output:
        tuple(sampleName, path("${sampleName}.primertrimmed.consensus.aln.utr_trimmed.fa"))

    script:
    awk_string = '/^>/ { printf("%s\\n", $0); next; } { printf("%s", $0); } END { printf("\\n"); }'
        """
        echo -e "\$(head -n 1 ${alignment} | cut -c 2-):266-29674" > non_utr.txt
        samtools faidx ${alignment}
        samtools faidx -r non_utr.txt ${alignment} > ${sampleName}.primertrimmed.consensus.aln.utr_trimmed.multi_line.fa
        awk '${awk_string}' ${sampleName}.primertrimmed.consensus.aln.utr_trimmed.multi_line.fa > ${sampleName}.primertrimmed.consensus.aln.utr_trimmed.fa
        """
}
