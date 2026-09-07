#!/mmfs1/home/mdfahad.hasan/miniconda3/envs/kallisto/bin/python

from pyfaidx import Fasta
# from idxgenome import pyfasta_genome
import sys
import os
import argcomplete
import argparse
import operator

### python ext_seq_from_annotation.py -genome ../../references/Barley/MorexV3/HvulgareMorex_702_V3.fa \\ 
### -gtf ../../references/Barley/MorexV3/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.chr.58.gtf \\ 
### -genes ../../zliang/projects/barley_zhaohui/DEG/filtered_set/filteredC1.csv.pos

def find_longest_transcript(trsdict):
    primary_GeneToTr = {}
    primary_TrToGene = {}
    for g in trsdict:
        sorted_x = sorted(trsdict[g].items(), key=operator.itemgetter(1), reverse=True)
        pgene = g
        ptranscript = sorted_x[0][0]
        if pgene not in primary_GeneToTr:
            primary_GeneToTr[pgene] = ptranscript
        if ptranscript not in primary_TrToGene:
            primary_TrToGene[ptranscript] = pgene
    return primary_GeneToTr, primary_TrToGene # gene <==> transcript

def sequence_extractor(genomeidx, chrom, start, stop, strand):
    outseq = ""
    seq = genomeidx[chrom][start:stop]
    if strand == "+":
        outseq = seq
    else:
        outseq = seq.reverse.complement
    return str(outseq)

#def parse_TSS_gff(gtf, mygenome, upext=2000, downext=500): # gtf file should be downloaded from ensembl
#def parse_TSS_gff(gtf, mygenome, upext=1000, downext=0):
def parse_TSS_gff(gtf, mygenome, upext=1000, downext=0):
    with open(gtf, 'r') as gfile:
        for line in gfile:
            if line.startswith("#"):
                continue
            new = line.strip().split("\t")
            if new[2] != "mRNA":
                continue
            attrs = dict(kv.split('=', 1) for kv in new[8].split(';') if '=' in kv)
            if attrs.get('tag') != 'Ensembl_canonical':
                continue

            chrom  = new[0]
            strand = new[6]
            start  = int(new[3])
            stop   = int(new[4])
            gene   = attrs['Parent'].split(':')[-1]

            if strand == "+":
                newstart = start - upext - 1
                newstop  = start + downext - 1
            else:
                newstart = stop - downext
                newstop  = stop + upext

            name  = ">" + gene + ":" + chrom + "_" + str(newstart) + "_" + str(newstop) + "_" + strand
            myseq = sequence_extractor(mygenome, chrom, newstart, newstop, strand)
            sys.stdout.write(name + "\n" + myseq + "\n")

def main():
    parser = argparse.ArgumentParser(description="extract TSS sequence")
    parser.add_argument('-gtf', required=False, help='gtf annotation file downloaded from Ensembl')
    parser.add_argument('-genome', required=True, help='genome file downloaded from Ensembl/Phytozome')
    # Create subparsers for gene and region
    subparsers = parser.add_subparsers(dest='genes', required=False, help='Select either "gene" or "region"')

    # Subparser for "gene"
    gene_parser = subparsers.add_parser('gene', help='Options for gene')
    gene_parser.add_argument('--gene-file', required=True, type=str, help='Path to the gene file')

    # Subparser for "region"
    region_parser = subparsers.add_parser('region', help='Options for region')
    region_parser.add_argument('--region-file', required=True, type=str, help='Path to the region file')

    # argcomplete.autocomplete(parser)
    args = parser.parse_args()
    genomefile = args.genome
    gtffile = args.gtf
    # genes = args.genes
    mygenome = Fasta(args.genome)
    print (mygenome)
    if args.genes:
        if args.genes == "gene":
            genefile = args.gene_file
            candidate_genes = set([])
            with open(genefile, 'r') as fh:
                for line in fh:
                    new = line.strip()
                    candidate_genes.add(new)
            parse_TSS_gff(gtffile, mygenome, filename=genefile, genelst=candidate_genes, allgene=False)
        elif args.genes == "region": # chr1.222.A.T
            regionfile = args.region_file
            with open(regionfile, 'r') as fh, open(regionfile+'.out', 'w') as out:
                for line in fh:
                    new = line.strip()
                    n = new.split('.')
                    chrom = n[0]
                    start = int(n[1])-100
                    stop = int(n[1])+100
                    seq = sequence_extractor(mygenome, chrom, start, stop, "+")
                    out.write(">"+chrom+':'+str(start)+'-'+str(stop)+'\n')
                    out.write(seq+'\n')
    else:
        parse_TSS_gff(gtffile, mygenome)

if __name__ == "__main__":
    main()
