"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Wanding Zhou, Tenghui Chen, Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

from .transcripts import *

def parse_ucsc_refgene(map_file, name2gene):

    """ start 1-based, end 1-based """

    cnt = 0
    for line in opengz(map_file):
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        if fields[13] != 'cmpl' or fields[14] != 'cmpl':
            continue
        gene_name = fields[12].upper()
        if gene_name in name2gene:
            gene = name2gene[gene_name]
        else:
            gene = Gene(name=gene_name)
            name2gene[gene_name] = gene
        t = Transcript()
        t.name = fields[1]
        t.chrm = normalize_chrm(fields[2])
        t.strand = fields[3]
        t.beg    = int(fields[4])+1
        t.end    = int(fields[5])
        t.cds_beg = int(fields[6])+1
        t.cds_end = int(fields[7])
        t.source = 'UCSC_refGene'
        ex_begs, ex_ends = fields[9], fields[10]

        for ex_beg, ex_end in zip([int(x)+1 for x in ex_begs.strip(',').split(',')],
                                  list(map(int, ex_ends.strip(',').split(',')))):
            t.exons.append((ex_beg, ex_end))
            
        t.exons.sort() # keep exons sorted
        gene.tpts.append(t)
        t.gene = gene
        cnt += 1

    err_print('loaded %d transcripts from UCSC refgene.' % cnt)

    return

def parse_ucsc_refgene_customized(map_file, name2gene):

    """ start 1-based, end 1-based """
    cnt = 0
    for line in open(map_file):
        fields = line.strip().split()
        gene_name = fields[0].upper()
        if gene_name in name2gene:
            gene = name2gene[gene_name]
        else:
            gene = Gene(name=gene_name)
            name2gene[gene_name] = gene

        t = Transcript()
        t.chrm = normalize_chrm(fields[1])
        t.strand = fields[2]
        t.beg    = int(fields[3])
        t.end    = int(fields[4])
        t.seq    = fields[-1]
        t.cds_beg = int(fields[5])
        t.cds_end = int(fields[6])
        t.source = 'custom'
        t.name = '.'
        ex_begs, ex_ends = fields[8], fields[9]

        for ex_beg, ex_end in zip(list(map(int, ex_begs.split(','))),
                                  list(map(int, ex_ends.split(',')))):
            t.exons.append((ex_beg, ex_end))
            
        t.exons.sort() # keep exons sorted
        gene.tpts.append(t)
        t.gene = gene
        cnt += 1

    err_print('loaded %d transcripts from customized table.' % cnt)

    return

def parse_refseq_gff(gff_fn, name2gene):

    id2ent = {}
    gff_fh = opengz(gff_fn)
    reg = None
    cnt = 0
    for line in gff_fh:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        # print line.strip()
        info = dict([_.split('=') for _ in fields[8].split(';')])
        if fields[2] == 'region':
            if 'chromosome' in info:
                reg = Region(info['chromosome'], int(fields[3]), int(fields[4]))
            if 'map' in info and info['map']=='unlocalized':
                reg.unlocalized = False
            # else:
            # reg = None
        elif (reg and fields[2] == 'gene' and
              ('pseudo' not in info or info['pseudo'] != 'true')):

            if reg.unlocalized:
                continue

            gene_name = info['Name'].upper()
            if gene_name in name2gene:
                g = name2gene[gene_name]
                if hasattr(g, '_gene_id') and g._gene_id != info['ID']:
                    continue   # if a gene_name appears twice, then all the subsequent occurrences are all ignored.
            else:
                g = Gene(name=gene_name)
                name2gene[gene_name] = g
            g._gene_id = info['ID']
            g.beg = int(fields[3])
            g.end = int(fields[4])
            id2ent[info['ID']] = g
            if 'Dbxref' in info:
                g.dbxref = info['Dbxref']

        elif (fields[2] in ['mRNA', 'ncRNA', 'rRNA', 'tRNA']
              and 'Parent' in info and info['Parent'] in id2ent):

            if reg.unlocalized:
                continue
            
            if fields[2] == 'mRNA':
                fields[2] = 'protein_coding'
            if fields[2] == 'ncRNA' and 'ncrna_class' in info:
                fields[2] = info['ncrna_class']
            t = Transcript(transcript_type=fields[2])
            t.chrm = normalize_chrm(reg.name)
            t.strand = fields[6]
            t.beg = int(fields[3])
            t.end = int(fields[4])
            t.name = info['Name'] if 'Name' in info else info['product']
            t.gene = id2ent[info['Parent']]
            t.gene.tpts.append(t)
            t.source = 'RefSeq'
            id2ent[info['ID']] = t
            cnt += 1
            
        elif fields[2] == 'exon' and info['Parent'] in id2ent:

            if reg.unlocalized:
                continue

            t = id2ent[info['Parent']]
            if (isinstance(t, Gene)):
                g = t
                if not hasattr(g, 'gene_t'):
                    g.gene_t = Transcript()
                    g.tpts.append(g.gene_t)
                    g.gene_t.chrm = normalize_chrm(reg.name)
                    g.gene_t.strand = fields[6]
                    g.gene_t.gene = g
                    g.gene_t.beg = g.beg
                    g.gene_t.end = g.end
                    g.gene_t.source = 'RefSeq'
                    cnt += 1
                t = g.gene_t
            t.exons.append((int(fields[3]), int(fields[4])))
        elif fields[2] == 'CDS' and info['Parent'] in id2ent:

            if reg.unlocalized:
                continue

            t = id2ent[info['Parent']]
            if (isinstance(t, Gene)):
                g = t
                if not hasattr(g, 'gene_t'):
                    g.gene_t = Transcript()
                    g.tpts.append(g.gene_t)
                    g.gene_t.chrm = normalize_chrm(reg.name)
                    g.gene_t.strand = fields[6]
                    g.gene_t.gene = g
                    g.gene_t.beg = g.beg
                    g.gene_t.end = g.end
                    g.gene_t.source = 'RefSeq'
                    cnt += 1
                t = g.gene_t
            t.cds.append((int(fields[3]), int(fields[4])))

    err_print("loaded %d transcripts from RefSeq GFF3 file." % cnt)

def parse_ensembl_gtf(gtf_fn, name2gene):
    """
    This parses the new GTF after or equal to hg19.
    The function does not handle hg18.
    gtf file is gffv2
    parser does not assume order in the GTF file
    """

    gtf_fh = opengz(gtf_fn)
    id2ent = {}
    cnt = 0
    for line in gtf_fh:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
        # info = dict([_.split('=') for _ in fields[8].split(';')])
        if fields[2] == 'gene':
            gene_id = info['gene_id']
            if gene_id not in id2ent:
                id2ent[gene_id] = Gene(gene_type=info['gene_biotype'])
            g = id2ent[gene_id]
            if 'gene_name' in info:
                g.name = info['gene_name'].upper()
            else:
                g.name = gene_id
            if g.name not in name2gene: name2gene[g.name] = g
            g.beg = int(fields[3])
            g.end = int(fields[4])
            
        elif fields[2] == 'transcript':

            # there exits two transcript format in ensembl gtf
            # the old version has no 'transcript_biotype'
            # the equivalent transcript_biotype is fields[1]
            tid = info['transcript_id']
            if tid not in id2ent: 
                transcript_type = info['transcript_biotype'] if 'transcript_biotype' in info else fields[1]
                id2ent[tid] = Transcript(transcript_type=transcript_type)
            t = id2ent[tid]
            t.chrm = normalize_chrm(fields[0])
            t.strand = fields[6]
            t.beg = int(fields[3])
            t.end = int(fields[4])
            t.name = info['transcript_id']
            gene_id = info['gene_id']
            if gene_id not in id2ent:
                id2ent[gene_id] = Gene(gene_type=info['gene_biotype'])
            t.gene = id2ent[gene_id]
            t.gene.tpts.append(t)
            t.source = 'Ensembl'
            cnt += 1
        elif fields[2] == 'exon':
            tid = info['transcript_id']
            if tid not in id2ent:
                transcript_type = info['transcript_biotype'] if 'transcript_biotype' in info else fields[1]
                id2ent[tid] = Transcript(transcript_type=transcript_type)
            t = id2ent[tid]
            t.exons.append((int(fields[3]), int(fields[4])))
        elif fields[2] == 'CDS':
            tid = info['transcript_id']
            if tid not in id2ent:
                transcript_type = info['transcript_biotype'] if 'transcript_biotype' in info else fields[1]
                id2ent[tid] = Transcript(transcript_type=transcript_type)
            t = id2ent[tid]
            t.cds.append((int(fields[3]), int(fields[4])))

    err_print("loaded %d transcripts from Ensembl GTF file." % cnt)

def parse_ensembl_gtf_hg18(gtf_fn, name2gene):

    """
    This parses the old ensembl GTF before or equal to hg18.
    The function does not handle hg19 or later.
    """

    gtf_fh = opengz(gtf_fn)
    tid2transcript = {}
    cnt = 0
    for line in gtf_fh:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
        if fields[2] == "exon":
            if info['transcript_id'] in tid2transcript:
                t = tid2transcript[info['transcript_id']]
            else:
                t = Transcript(transcript_type=fields[1])
                t.chrm = normalize_chrm(fields[0])
                t.strand = fields[6]
                t.name = info['transcript_id']
                tid2transcript[t.name] = t
                if info['gene_name'] in name2gene:
                    g = name2gene[info['gene_name']]
                else:
                    g = Gene()
                    g.name = info['gene_name']
                    name2gene[g.name] = g
                t.gene = g
                g.tpts.append(t)
                t.source = 'Ensembl'
                cnt += 1
            t.exons.append((int(fields[3]), int(fields[4])))
        elif fields[2] == 'CDS':
            if info['transcript_id'] in tid2transcript:
                t = tid2transcript[info['transcript_id']]
            else:
                t = Transcript(transcript_type=fields[1])
                t.chrm = normalize_chrm(fields[0])
                t.strand = fields[6]
                t.name = info['transcript_id']
                tid2transcript[t.name] = t
                if info['gene_name'] in name2gene:
                    g = name2gene[info['gene_name']]
                else:
                    g = Gene()
                    g.name = info['gene_name']
                    name2gene[g.name] = g
                t.gene = g
                g.tpts.append(t)
                t.source = 'Ensembl'
                cnt += 1
            t.cds.append((int(fields[3]), int(fields[4])))

    for t in list(tid2transcript.values()):
        t.exons.sort()
        t.beg = t.exons[0][0]
        t.end = t.exons[-1][1]

    err_print("loaded %d transcripts from Ensembl GTF file." % cnt)

def parse_ccds_table(ccds_fn, name2gene):

    """ start 0-based end 0-based """

    ccds_fh = open(ccds_fn)
    ccds_fh.readline()
    cnt = 0
    for line in ccds_fh:
        fields = line.strip().split('\t')
        if fields[5] != 'Public':
            continue
        gene_name = fields[2].upper()
        if gene_name not in name2gene:
            name2gene[gene_name] = Gene(name=gene_name)

        g = name2gene[gene_name]
        t = Transcript()
        t.chrm = normalize_chrm(fields[0])
        t.strand = fields[6]
        t.cds_beg = int(fields[7])+1
        t.cds_end = int(fields[8])+1

        # without UTR information, take CDS boundary as the exon boundary
        t.beg = t.cds_beg
        t.end = t.cds_end

        t.name = fields[4]
        # note that CCDS do not annotate UTR, so all the exons are equivalently cds
        t.exons = [(int(b)+1, int(e)+1) for b,e in re.findall(r"[\s\[](\d+)-(\d+)[,\]]", fields[9])]
        t.source = 'CDDS'
        t.gene = g
        g.tpts.append(t)
        cnt += 1

    err_print("loaded %d transcripts from CCDS table." % cnt)

def parse_ucsc_kg_table(kg_fn, alias_fn, name2gene):

    kg_fh = opengz(kg_fn)
    id2aliases = {}
    if alias_fn:
        alias_fh = opengz(alias_fn)
        for line in alias_fh:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if fields[0] in id2aliases:
                id2aliases[fields[0]].append(fields[1])
            else:
                id2aliases[fields[0]] = [fields[1]]

    cnt = 0
    for line in kg_fh:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        g = None
        if fields[0] in id2aliases:
            for alias in id2aliases[fields[0]]:
                if alias in name2gene:
                    g = name2gene[alias]
            if not g:
                g = Gene(name=fields[0])
            for alias in id2aliases[fields[0]]:
                name2gene[alias] = g
        else:
            if fields[0] in name2gene:
                g = name2gene[fields[0]]
            else:
                g = Gene(name=fields[0])
            name2gene[fields[0]] = g

        t = Transcript()
        t.name = fields[0]
        t.chrm = normalize_chrm(fields[1])
        t.strand = fields[2]
        t.beg = int(fields[3])
        t.end = int(fields[4])
        t.cds_beg = int(fields[5])
        t.cds_end = int(fields[6])
        t.source = 'UCSC_knownGene'
        ex_begs, ex_ends = fields[8], fields[9]
        for ex_beg, ex_end in zip(list(map(int, ex_begs.strip(',').split(','))),
                                  list(map(int, ex_ends.strip(',').split(',')))):
            t.exons.append((ex_beg, ex_end))
        t.exons.sort()
        g.tpts.append(t)
        t.gene = g
        cnt += 1

    err_print("loaded %d transcripts from UCSC knownGene table." % cnt)

def parse_gencode_gtf(gencode_fn, name2gene):

    id2ent = {}
    gencode_fh = opengz(gencode_fn)
    cnt = 0
    for line in gencode_fh:
        # if cnt > 1000:
        #     break
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
        if fields[2] == 'gene':
            gene_name = info['gene_name'].upper()
            gid = info['gene_id']
            if gene_name in name2gene:
                g = name2gene[gene_name]
                id2ent[gid] = g
            else:
                if gid not in id2ent:
                    id2ent[gid] = Gene(name=gene_name, gene_type=info['gene_type'])
                g = id2ent[gid]
                name2gene[gene_name] = g
            g.beg = int(fields[3])
            g.end = int(fields[4])
            # if info['gene_type'] == 'pseudogene':
            #     g.pseudo = True
        elif fields[2] == 'transcript':
            tid = info['transcript_id']
            if tid not in id2ent:
                id2ent[tid] = Transcript(transcript_type=info['transcript_type'])
                
            t = id2ent[tid]
            t.chrm = normalize_chrm(fields[0])
            t.strand = fields[6]
            t.beg = int(fields[3])
            t.end = int(fields[4])
            t.name = tid
            gid = info['gene_id']
            if gid not in id2ent:
                id2ent[gid] = Gene(gene_type=info['gene_type'])
            t.gene = id2ent[gid]
            t.gene.tpts.append(t)
            t.source = 'GENCODE'
            id2ent[t.name] = t
            cnt += 1
        elif fields[2] == 'exon':
            tid = info['transcript_id']
            if tid not in id2ent:
                id2ent[tid] = Transcript(transcript_type=info['transcript_type'])
            t = id2ent[tid]
            t.exons.append((int(fields[3]), int(fields[4])))
        elif fields[2] == 'CDS':
            tid = info['transcript_id']
            if tid not in id2ent:
                id2ent[tid] = Transcript(transcript_type=info['transcript_type'])
            if tid not in id2ent: id2ent[tid] = Transcript()
            t = id2ent[tid]
            t.cds.append((int(fields[3]), int(fields[4])))

    err_print("loaded %d transcripts from GENCODE GTF file." % cnt)

def parse_aceview_transcripts(aceview_gff_fn, name2gene):

    id2tpt = {}
    aceview_fh = opengz(aceview_gff_fn)
    for line in aceview_fh:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        if len(fields) < 9: continue # the old transcript definition (hg18) from AceView is a bit corrupted.
        info = dict(re.findall(r'\s*(\S+) (\S+);', fields[8]))
        if fields[2] == 'CDS':
            gene_name = info['gene_id'].upper()
            if gene_name in name2gene:
                g = name2gene[gene_name]
            else:
                g = Gene(name=gene_name)
                name2gene[gene_name] = g

            if info['transcript_id'] in id2tpt:
                t = id2tpt[info['transcript_id']]
            else:
                t = Transcript()
                t.chrm = normalize_chrm(fields[0])
                t.strand = fields[6]
                t.name = info['transcript_id']
                id2tpt[t.name] = t
                t.gene = g
                g.tpts.append(t)
                t.source = 'AceView'

            t.cds.append((int(fields[3]), int(fields[4])))

        elif fields[2] == 'exon':
            gene_name = info['gene_id'].upper()
            if gene_name in name2gene:
                g = name2gene[gene_name]
            else:
                g = Gene(name=gene_name)
                name2gene[gene_name] = g

            if info['transcript_id'] in id2tpt:
                t = id2tpt[info['transcript_id']]
            else:
                t = Transcript()
                t.chrm = normalize_chrm(fields[0])
                t.strand = fields[6]
                t.name = info['transcript_id']
                id2tpt[t.name] = t
                t.gene = g
                g.tpts.append(t)
                t.source = 'AceView'

            t.exons.append((int(fields[3]), int(fields[4])))

    # skip transcripts without CDS, e.g., LOC391566.aAug10-unspliced
    for tid, t in id2tpt.items():
        if t.cds and t.exons:
            t.exons.sort()
            t.beg = t.exons[0][0]
            t.end = t.exons[-1][1]
        else:
            t.gene.tpts.remove(t)

    err_print("loaded %d transcripts from AceView GFF file." % len(id2tpt))

def parse_uniprot_mapping(fn):

    tid2uniprot = {}
    for line in opengz(fn):
        fields = line.strip().split('\t')
        if fields[2] != fields[0]:
            # tid2uniprot[fields[2]] = fields[0]
            if fields[0] in tid2uniprot:
                tid2uniprot[fields[0]].append(fields[2])
            else:
                tid2uniprot[fields[0]] = [fields[2]]

    err_print('loaded %d transcript with UniProt mapping.' % len(tid2uniprot))

    return tid2uniprot


def parse_annotation(args):

    name2gene = {}
    if args.ensembl:
        if args.refversion in ['hg18']:
            parse_ensembl_gtf_hg18(args.ensembl, name2gene)
        else:
            parse_ensembl_gtf(args.ensembl, name2gene)
    
    if args.refseq:
        parse_refseq_gff(args.refseq, name2gene)

    if args.ccds:
        parse_ccds_table(args.ccds, name2gene)

    if args.gencode:
        parse_gencode_gtf(args.gencode, name2gene)
        # try:
        #     import pysam
        #     args.ffhs['GENCODE'] = pysam.Tabixfile(args.gencode)
        # except:
        #     err_print("Cannot import non-coding features (may need pysam).")

    if args.ucsc:
        parse_ucsc_refgene(args.ucsc, name2gene)

    # if args.custom:
    #     parse_ucsc_refgene_customized(args.custom, name2gene)

    if args.kg:
        parse_ucsc_kg_table(args.kg, args.alias, name2gene)

    if args.aceview:
        parse_aceview_transcripts(args.aceview, name2gene)

    # remove genes without transcripts
    names_no_tpts = []
    for name, gene in name2gene.items():
        # print gene, len(gene.tpts)
        if not gene.tpts:
            names_no_tpts.append(name)
    for name in names_no_tpts:
        del name2gene[name]
    err_print('loaded %d genes.' % len(name2gene))

    # index transcripts in a gene
    thash = THash()
    genes = set(name2gene.values())
    for g in genes:
        for t in g.tpts:
            t.exons.sort()
            if not (hasattr(t, 'cds_beg') and hasattr(t, 'cds_end')):
                if t.cds:
                    t.cds.sort()
                    t.cds_beg = t.cds[0][0]
                    t.cds_end = t.cds[-1][1]
                elif hasattr(t,'beg') and hasattr(t,'end'):
                    t.cds_beg = t.beg
                    t.cds_end = t.end
                else:
                    t.cds_beg = t.exons[0][0]
                    t.cds_end = t.exons[-1][1]

            thash.insert(t)
            if len(t.exons) == 0:     # if no exon, use cds
                t.exons = t.cds[:]

        g.std_tpt = g.longest_tpt()

    if args.uniprot:
        tid2uniprot = parse_uniprot_mapping(args.uniprot)
        name2protein = {}
        for name, gene in name2gene.items():
            for tpt in gene.tpts:
                if tpt.name in tid2uniprot:
                    uniprot = tid2uniprot[tpt.name]
                    if uniprot not in name2protein:
                        name2protein[uniprot] = Gene(uniprot)
                    prot = name2protein[uniprot]
                    prot.tpts.append(tpt)
        return name2protein, thash
    else:
        return name2gene, thash

def parser_add_annotation(parser):

    parser.add_argument("--longest", action="store_true",
                        help="consider only longest transcript")
    parser.add_argument("--longestcoding", action="store_true",
                        help="consider only protein-coding transcript with longest cds")
    parser.add_argument('--refversion', nargs='?', default=None,
                        help='reference version (hg18, hg19, hg38 etc) (config key: refversion)')
    parser.add_argument('--reference', nargs='?', default='_DEF_',
                        help='indexed reference fasta (with .fai) (config key: reference)')
    parser.add_argument('--ensembl', nargs='?', default=None, const='_DEF_',
                        help='Ensembl GTF transcript annotation (config key: ensembl)')
    parser.add_argument('--gencode', nargs='?', default=None, const='_DEF_',
                        help='GENCODE GTF transcript annotation (config key: gencode)')
    parser.add_argument('--kg', nargs='?', default=None, const='_DEF_',
                        help='UCSC knownGene transcript annotation (config key: known_gene)')
    parser.add_argument('--alias', nargs='?', default=None, const='_DEF_',
                        help='UCSC knownGene aliases (without providing aliases, only the knownGene id can be searched (config key: known_gene_alias)')
    parser.add_argument('--ucsc', nargs='?', default=None, const='_DEF_',
                        help='UCSC transcript annotation table (config key: ucsc')
    parser.add_argument('--refseq', nargs='?', default=None, const='_DEF_',
                        help='RefSeq transcript annotation (config key: refseq)')
    parser.add_argument('--ccds', nargs='?', default=None, const='_DEF_',
                        help='CCDS transcript annotation table (config key: ccds)')
    parser.add_argument('--aceview', nargs='?', default=None, const='_DEF_',
                        help='AceView GFF transcript annotation (config key: aceview)')
    # parser.add_argument
    # parser.add_argument('--custom', nargs='?', default=None, const='_DEF_',
    #                     help='A customized transcript table with sequence (config key: custom)')
    parser.add_argument('--uniprot', nargs='?', default=None, const='_DEF_',
                        help='use uniprot ID rather than gene id (config key: uniprot)')
    parser.add_argument('--mem', action='store_true', help='for processing large input, preload indices')
    parser.add_argument('--sql', action='store_true', help='SQL mode')
    parser.add_argument('--prombeg', type=int, default=1000, help='promoter starts from n1 bases upstream of transcription start site (default: n1=1000)')
    parser.add_argument('--promend', type=int, default=0, help='promoter ends extends to n2 bases downstream of transcription start site (default: n2=0)')

    return
