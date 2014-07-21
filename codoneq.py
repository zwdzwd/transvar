


def main(args):

    codon1, codon2 = args.codons
    name2gene = parse_annotation(args.annotation_file)

    gene_name, codon_pos = codon1.split('.p')
    print gene_name, codon_pos
    codon_pos = int(codon_pos)
    gene = name2gene[gene_name]
    locs1 = set()
    for i, trans in enumerate(gene.transcripts):
        aa_pos = trans.aa_pos2nuc_pos(codon_pos)
        if aa_pos: 
            print "transcript %d \tcodon: %s" % (i, "-".join(map(str,aa_pos)))
            locs1.add(tuple(aa_pos))

    gene_name, codon_pos = codon2.split('.p')
    print gene_name, codon_pos
    codon_pos = int(codon_pos)
    gene = name2gene[gene_name]
    locs2 = set()
    for i, trans in enumerate(gene.transcripts):
        aa_pos = trans.aa_pos2nuc_pos(codon_pos)
        if aa_pos: 
            print "transcript %d \tcodon: %s" % (i, "-".join(map(str,aa_pos)))
            locs2.add(tuple(aa_pos))

    if locs1 & locs2:
        print "Genomic location might be the same."
    else:
        print "Genomic location shouldn't be the same."


def add_parser_codoneq(subparsers):

    parser = subparsers.add_parser("codoneq", help="test whether two codon annotations may be the same nucleotide position", epilog="Example: prog -c MET.p1010 MET.p992 -a hg19.map")
    parser.add_argument("-c", nargs=2, help="two codons to test. Format: [gene_name].p[codon_location], e.g., MET.p1010", dest='codons')
    parser.add_argument('-a', metavar='annotation',
                        required = True,
                        dest='annotation_file', 
                        help='protein annotation file')
    parser.set_defaults(func=main)
