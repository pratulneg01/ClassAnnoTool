nextflow.enable.dsl=2
params.input = null

workflow {
    input_file = file(params.input)

    cutadapt_result = run_cutadapt(input_file, file("bin/cutadapt.sh"))

    cutadapt_result
        .flatten()
        .filter { it.name.endsWith('.fastq') }
        .set { fastq_output }

    fastqc_result = run_fastqc(fastq_output, file("bin/fastqc.sh"))

    fastqc_result
        .flatten()
        .filter { it.name.endsWith('.zip') }
        .set { fastqc_zips }

    fastq_output
        .combine(fastqc_zips)
        .set { filtered_inputs }

    filtered_fastq = filter_fastq(filtered_inputs, file("bin/filter_fastq_by_fastqc.py"))
    fastas_from_fastq = fastq_to_fasta(filtered_fastq, file("bin/fastq_to_fasta.sh"))

    cutadapt_result
        .flatten()
        .filter { it.name.endsWith('.fasta') }
        .set { fasta_output }

    cleaned_fastas = clean_fasta(fasta_output, file("bin/clean_fasta.py"))

    all_fastas = cleaned_fastas.mix(fastas_from_fastq)

    orf_results = find_orfs(all_fastas, file("bin/orffinder.sh"))

    orf_results
        .flatten()
        .filter { it.toString().endsWith('.fasta') }
        .set { blastn_input }

    blastn_results = blastn_remote(blastn_input, file("bin/blastn.sh"))

    orf_results
        .flatten()
        .filter { it.name.endsWith('_orfs.fasta') }
        .set { orfs_fasta }

    (annotated_fasta, unannotated_fasta, _, _) = classify_orfs_by_blast(
    blastn_results,
    orfs_fasta,
    file("bin/classify_orf_by_blast.py"))

    translated = translate_unannotated(unannotated_fasta, file("bin/translate_to_aa.py"))

    blastp_output = blastp_local(translated, file("database"), file("bin/local_blastp.sh"))

    (unannotated_prot_fasta, _) = classify_proteins_by_blast(
    translated,
    blastp_output,
    file("bin/classify_prot_by_blast.py"))


    extract_features_result = extract_features(
    unannotated_fasta,
    file("bin/extract_caracteristics.py"))

    analysis_results = analyze_features_ml(
    extract_features_result,
    file("bin/features_clasification.R")
)

}

process run_cutadapt {
    input:
    path input_file
    path cutadapt_script

    output:
    path "results/cutadapt/*"

    publishDir "results/cutadapt", mode: 'copy'

    script:
    """
    bash ${cutadapt_script} ${input_file}
    """
}

process run_fastqc {
    input:
    path fastq_file
    path fastqc_script

    output:
    path "results/fastqc/*"

    publishDir "results/fastqc", mode: 'copy'

    script:
    """
    bash ${fastqc_script} ${fastq_file}
    """
}

process filter_fastq {
    input:
    tuple path(fastq_file), path(fastqc_zip)
    path filtering_script

    output:
    path "results/fastqc/*_filtered.fastq" 

    publishDir "results/fastqc", mode: 'copy'

    script:
    """
    base=\$(basename ${fastq_file} .fastq)
    mkdir -p results/fastqc
    python3 ${filtering_script} ${fastqc_zip} ${fastq_file} results/fastqc/\${base}_filtered.fastq
    """
}

process fastq_to_fasta {
    input:
    path fastq_file 
    path script

    output:
    path "results/fastq_to_fasta/*.fasta"

    publishDir "results/fastq_to_fasta", mode: 'copy'

    script:
    """
    base=\$(basename ${fastq_file} .fastq)
    mkdir -p results/fastq_to_fasta
    bash ${script} ${fastq_file} results/fastq_to_fasta/\${base}.fasta
    """
}

process clean_fasta {
    conda 'envs/biopython.yaml'

    input:
    path fasta_file
    path cleaning_script

    output:
    path "results/fasta_cleaned/*.fasta"

    publishDir "results/fasta_cleaned", mode: 'copy'

    script:
    """
    mkdir -p results/fasta_cleaned
    output_name=\$(basename ${fasta_file} .fasta)_cleaned.fasta
    python3 ${cleaning_script} ${fasta_file} results/fasta_cleaned/\$output_name
    """
}

process find_orfs {
    input:
    path fasta_file
    path script

    output:
    path "results/orfs/*"

    publishDir "results/orfs", mode: 'copy'

    script:
    """
    mkdir -p results/orfs
    base=\$(basename ${fasta_file} .fasta)
    bash ${script} ${fasta_file}

    # Mover el archivo resultante a results/orfs
    mv results/orfs/orfipy_\${base}.fasta_out/\${base}_orfs.fasta results/orfs/
    """
}

process blastn_remote {
    input:
    path fasta_file
    path script

    output:
    path "results/blastn/*.txt"

    publishDir "results/blastn", mode: 'copy'

    script:
    """
    mkdir -p results/blastn
    base=\$(basename ${fasta_file} .fasta)
    bash ${script} ${fasta_file} results/blastn/\${base}_blastn.txt
    """
}

process classify_orfs_by_blast {

    input:
    path blast_file
    path fasta_file
    path script

    output:
    path "annotated.fasta"
    path "unannotated.fasta"
    path "annotated_ids.txt"
    path "unannotated_ids.txt"

    publishDir "results/blast_class", mode: 'copy'

    script:
    """
    python3 ${script} ${blast_file} ${fasta_file} --identity 80 --coverage 70 --evalue 1e-5
    """
}

process translate_unannotated {

    input:
    path unannotated_fasta
    path script

    output:
    path "unannotated_translated.fasta"

    publishDir "results/protein_sequence", mode: 'copy'

    script:
    """
    python3 ${script} ${unannotated_fasta} unannotated_translated.fasta
    """
}

process blastp_local {
    input:
    path fasta_file
    path db_dir
    path script

    output:
    path "results/blastp/protein_blastp.txt"

    publishDir "results/blastp", mode: 'copy'

    script:
    """
    mkdir -p results/blastp
    cp -r ${db_dir}/* .  # Copiar todos los archivos de la DB al work dir
    bash $script $fasta_file swissprot results/blastp/protein_blastp.txt
    """
}

process classify_proteins_by_blast {
    input:
    path translated_fasta
    path blastp_result
    path script

    output:
    path "annotated_proteins.fasta"
    path "unannotated_proteins.fasta"

    publishDir "results/protein_class", mode: 'copy'

    script:
    """
    python3 ${script} ${translated_fasta} ${blastp_result} --identity 70 --evalue 1e-3
    """
}

process extract_features {
    input:
    path unannotated_prot_fasta
    path script

    output:
    path "results/characteristics/features.txt"

    publishDir "results/characteristics", mode: 'copy'

    script:
    """
    mkdir -p results/characteristics
    python3 ${script} ${unannotated_prot_fasta} results/characteristics/features.txt
    """
}

process analyze_features_ml {
    input:
    path features_txt
    path script

    output:
    path "results/analysis/*"

    publishDir "results/analysis", mode: 'copy'

    script:
    """
    mkdir -p results/analysis
    Rscript $script --features $features_txt --output results/analysis
    """
}