def main(args):

    name2gene = parse_annotation(args.annotation_file)
    # print name2gene['TP53'].transcripts[0].aa_pos2nuc_pos(348), name2gene['TP53'].transcripts[0].cds_beg, name2gene['TP53'].transcripts[0].cds_end, name2gene['TP53'].transcripts[0].exons
    # print name2gene['TP53'].transcripts[1].aa_pos2nuc_pos(393), name2gene['TP53'].transcripts[1].cds_beg, name2gene['TP53'].transcripts[1].cds_end, name2gene['TP53'].transcripts[1].exons
    with open(args.codon_list) as f:
        for line in f:
            gene_name, codon_pos = line.strip().split(':')

            # if line.strip() != "TP53:281":
                # continue

            gene = name2gene[gene_name]
            codon_pos = int(codon_pos)
            loc2trans = {}
            locs = set()
            for i, trans in enumerate(gene.transcripts):
                nucpos = trans.aa_pos2nuc_pos(codon_pos)
                if nucpos:
                    locs.add(tuple(nucpos))

                    if tuple(nucpos) in loc2trans:
                        loc2trans[tuple(nucpos)].append(trans)
                    else:
                        loc2trans[tuple(nucpos)] = [trans]

                    # print trans, trans.cds_beg, trans.cds_end, nucpos, trans.exons

            aa_poses = set()
            aa_pos2trans = {}
            for loc in locs:
                for i, trans in enumerate(gene.transcripts):
                    result = trans.nuc_pos2aa_pos(loc[0])
                    if not result:
                        continue
                    aa_pos, codon_loc = result
                    # print i, trans, "nucleotide loc: ", loc, "aa loc: ", result

                    # if aa_pos <=0:
                    #     print trans, trans.cds_beg, trans.cds_end, loc, aa_pos
                    #     print trans.exons
                    #     import sys
                    #     sys.exit(1)

                    if aa_pos != codon_pos:
                        # aa_poses.add("%s:%d:%d" % (gene_name, aa_pos, codon_loc))
                        aa_poses.add("%s.p%d" % (gene_name, aa_pos))

                        aa_pos_str = "%s.p%d" % (gene_name, aa_pos)
                        if aa_pos_str in aa_pos2trans:
                            aa_pos2trans[aa_pos_str].append((trans))
                        else:
                            aa_pos2trans[aa_pos_str] = [(aa_pos2trans, trans)]


            print "%s.p%s\t%s\t%d\t%s" % (gene_name, codon_pos, gene.strand(), len(gene.transcripts), '\t'.join(aa_poses))


def add_parser_codonsearch(subparsers):

    parser = subparsers.add_parser('codonsearch', help='search alternative codonpositions due to different transcript usage')
    parser.add_argument('-c', dest='codon_list', help='a list of codons to search alternatives')
    parser.add_argument('-a',
                        metavar='annotation',
                        required = True,
                        dest='annotation_file', 
                        help='protein annotation file')
    parser.set_defaults(func=main)
