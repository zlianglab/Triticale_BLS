import subprocess as sp
import argparse
import os
import sys
import glob

def get_args():
    parser = argparse.ArgumentParser(
        description = 'RNASeq pipeline',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--annot', '-a',
                        metavar = 'GFF',
                        type = str,
                        required = True,
                        help = 'Annotation file (GFF3/GTF)'

    )

    parser.add_argument('--f1', '-1',
                        metavar = 'FQ1',
                        type = str,
                        required = True,
                        help = 'Fastq1'
                        
    )

    parser.add_argument('--f2', '-2',
                        metavar = 'FQ2',
                        type = str,
                        required = True,
                        help = 'Fastq2'
                        
    )

    parser.add_argument('--index', '-i',
                        metavar = 'INDEX',
                        type = str,
                        required = True,
                        help = 'Hisat2 Index files'
                        )

    parser.add_argument('--out', '-o',
                        metavar = 'OUT',
                        type = str,
                        required = True,
                        help = 'Output directory'
                        
    )

    parser.add_argument('--thread', '-t',
                        type = int,
                        default = 8,
                        help = 'Threads for fastp/hisat2/samtools'

    )

    return parser.parse_args()

def run_command(command):
    try:
        sp.run(command, check = True, shell = True)
    except sp.CalledProcessError as e:
        print(f'Error executing command: {e.cmd}')
        print(e.output)
        exit(1)

def fastp_trim(fq1, fq2, out1, out2, json, html, thread):
    command = (
        f'fastp -i {fq1} -I {fq2} -o {out1} -O {out2} -j {json} -h {html} -w {thread}'
    )
    run_command(command)

def hisat2_index(ref, index_out, thread):
    index_base = os.path.join(index_out, 'genome')
    if not os.path.exists(index_base + '.1.ht2'):
        command = f'hisat2-build -p {thread} {ref} {index_base}'
        run_command(command)
    else:
        print('Index already exists, skipping build.')
    return index_base

def hisat2_align(fq1, fq2, index, sam, thread):
    command = (
        f'hisat2 -x {index} -1 {fq1} -2 {fq2} -p {thread} | '
        f'samtools sort -@ {thread} -o {sam}'
    )
    run_command(command)

def index_bam(sam, thread):
    command = (
        f'samtools index -@ {thread} {sam}'
    )
    run_command(command)

def htseq_count(sam, annotation, out):
    command = (
        f'htseq-count -f bam -r pos -s no -t exon -i gene_id --nonunique=all '
        f'{sam} {annotation} > {out}'
    )
    run_command(command)

def main():
    args = get_args()

    annotation = args.annot
    out = args.out
    fastq1 = args.f1
    fastq2 = args.f2
    thread = args.thread
    index = args.index

    fastq_out = os.path.join(out, 'fastq')
    sam_out = os.path.join(out, 'sam')
    count_out = os.path.join(out, 'count')
    # index_out = os.path.join(out, 'index')

    for d in (fastq_out, sam_out, count_out):
        os.makedirs(d, exist_ok = True)

    # Hisat2 indexing
    index_base = hisat2_index(ref, index_out, thread)
    
    # Setting up trimmed fastq directories
    trim1 = os.path.join(fastq_out, f'{os.path.basename(fastq1).replace('_S_1.fq.gz', '')}_1.trim.fq.gz')
    trim2 = os.path.join(fastq_out, f'{os.path.basename(fastq2).replace('_S_2.fq.gz', '')}_2.trim.fq.gz')
    json = os.path.join(fastq_out, f'{os.path.basename(fastq1).replace('_S_1.fq.gz', '')}.fastp.json')
    html = os.path.join(fastq_out, f'{os.path.basename(fastq1).replace('_S_1.fq.gz', '')}.fastp.html')

    # trimming using fastp
    fastp_trim(fastq1, fastq2, trim1, trim2, json, html, thread)

    # Setting up sam directory
    sam = os.path.join(sam_out, f'{os.path.basename(fastq1).replace('_S_1.fq.gz', '')}.sorted.bam')

    # hisat2 align to sam files
    hisat2_align(trim1, trim2, index, sam, thread)

    # samtools indexing
    index_bam(sam, thread)

    # Setting up counts directory
    count_txt = os.path.join(count_out, f'{os.path.basename(fastq1).replace('_S_1.fq.gz', '')}_counts.txt')

    # htseq count
    htseq_count(sam, annotation, count_txt)

if __name__ == '__main__':
    main()