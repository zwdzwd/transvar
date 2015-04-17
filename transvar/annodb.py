from transcripts import *
import parser

class AnnoDB():
    
    def __init__(self, args, config):

        if args.refversion:
            self.rv = args.refversion
        elif 'refversion' in config.defaults():
            self.rv = config.get('DEFAULT', 'refversion')
        else:
            err_warn('please specify reference version, either in transvar.cfg, or as argument --refversion')
            sys.exit(1)
        
        replace_defaults(args, config)
        
        faidx.init_refgenome(args.reference if args.reference else None)
        self.session = None
        self.name2gene = None
        self.thash = None
        # args.ffhs = {}
        # if args.dbsnp:
        #     import pysam
        #     args.dbsnp_fh = pysam.Tabixfile(args.dbsnp)
        # else:
        #     args.dbsnp_fh = None
        
        if args.sql:
            import sqlmodel
            self.sqlmodel = sqlmodel
            self.session = sqlmodel.sessionmaker(bind=sqlmodel.engine, autoflush=False)()
            refversion_ids = self.session.query(sqlmodel.DRefVersion).filter_by(name=args.refversion).all()
            if not refversion_ids:
                raise Exception("No reference version named: %s" % args.refversion)
            self.refversion_id = refversion_ids[0].id
            self.source = []
            if args.ensembl:
                self.source.append('Ensembl')
            if args.ccds:
                self.source.append('CCDS')
            if args.refseq:
                self.source.append('RefSeq')
            if args.gencode:
                self.source.append('GENCODE')
            if args.aceview:
                self.source.append('AceView')
            if args.ucsc:
                self.source.append('UCSC')
            if not self.source:
                self.source.append('Ensembl')
        else:
            self.name2gene, self.thash = parser.parse_annotation(args)

        self.config = config
        self.args = args
        self.resources = {}
        self.init_resource()

    def init_resource(self):

        for rname in ['dbsnp']:
            if self.config.has_option(self.rv, 'dbsnp'):
                import tabix
                self.resources['dbsnp'] = tabix.open(self.config.get(self.rv, 'dbsnp'))


    def _query_dbsnp_(self, chrm, beg, end, ref=None, alt=None):

        dbsnps = []
        if 'dbsnp' in self.resources:
            if beg == end and (alt is None or len(alt)==1): # SNV
                ret = self.resources['dbsnp'].query(normalize_chrm_dbsnp(chrm), int(beg), int(end))
                for fields in ret:
                    if int(fields[1]) != int(beg):
                        continue
                    alts = fields[4].split(',')
                    if ref is not None and ref != fields[3]:
                        continue
                    if alt is not None and alt not in alts:
                        continue
                    if alt is None:
                        for alt in alts:
                            dbsnps.append('%s(%s:%s%s>%s)' % (fields[2], chrm, fields[1], fields[3], alt))
                    else:
                        dbsnps.append('%s(%s:%s%s>%s)' % (fields[2], chrm, fields[1], fields[3], alt))
            else:
                ret = self.resources['dbsnp'].query(normalize_chrm_dbsnp(chrm), int(beg)-1, int(end))
                for fields in ret:
                    if int(fields[1]) != int(beg)-1:
                        continue
                    alts = [_[1:] for _ in fields[4].split(',')]
                    if ref is not None and ref != fields[3][1:]:
                        continue
                    if alt is not None and alt not in alts:
                        continue
                    if alt is None:
                        for alt in alts:
                            dbsnps.append('%s(%s:%s%s>%s)' % (fields[2], chrm, fields[1], fields[3], fields[3][0]+alt))
                    else:
                        dbsnps.append('%s(%s:%s%s>%s)' % (fields[2], chrm, fields[1], fields[3], fields[3][0]+alt))
        return dbsnps

    def query_dbsnp_range(self, r, beg, end, alt):

        dbsnps = self._query_dbsnp_(r.chrm, beg, end, alt=alt)
        if dbsnps:
            r.append_info('dbsnp='+','.join(dbsnps))
        
    def query_dbsnp_codon(self, r, codon, taa_alt):

        dbsnps = []
        for tnuc_altseq in reverse_codon_table[taa_alt]:
            subs = []
            for i, (tnuc_refbase, tnuc_altbase) in enumerate(zip(codon.seq, tnuc_altseq)):
                if tnuc_refbase != tnuc_altbase:
                    subs.append((i, tnuc_altbase))
                    
            if len(subs) == 1:
                i, tnuc_altbase = subs[0]
                gnuc_pos = codon.tloc(i)
                gnuc_alt = tnuc_altbase if codon.strand == '+' else complement(tnuc_altbase)
                dbsnps.extend(self._query_dbsnp_(r.chrm, gnuc_pos, gnuc_pos, alt=gnuc_alt))

            if len(subs) == 2:
                i1 = subs[0][0]
                i2 = subs[1][0]
                if i1+1 == i2 and abs(codon.tloc(i1)-codon.tloc(i2)) == 1:
                    gnuc_beg = min(codon.tloc(i1), codon.tloc(i2))
                    gnuc_end = max(codon.tloc(i1), codon.tloc(i2))
                    tnuc_alt = subs[0][1]+subs[1][1]
                    gnuc_alt = tnuc_alt if codon.strand == '+' else reverse_complement(tnuc_alt)
                    dbsnps.extend(self._query_dbsnp_(r.chrm, gnuc_beg, gnuc_end, gnuc_alt))

            if len(subs) == 3:
                i1 = subs[0][0]
                i2 = subs[2][0]
                if i1+2 == i2 and abs(codon.tloc(i1)-codon.tloc(i2)) == 2:
                    gnuc_beg = min(codon.tloc(i1), codon.tloc(i2))
                    gnuc_end = max(codon.tloc(i1), codon.tloc(i2))
                    tnuc_alt = codon.seq
                    gnuc_alt = codon.seq if codon.strand == '+' else reverse_complement(codon.seq)
                    dbsnps.extend(self._query_dbsnp_(r.chrm, gnuc_beg, gnuc_end, gnuc_alt))

        if dbsnps:
            r.append_info('dbsnp='+','.join(dbsnps))
        
    def query_dbsnp(self, r, pos, ref=None, alt=None):

        dbsnps = self._query_dbsnp_(r.chrm, pos, pos, ref, alt)
        if dbsnps:
            r.append_info('dbsnp='+','.join(dbsnps))

    def get_gene(self, name):

        if self.session:
            gs = self.session.query(self.sqlmodel.DGene).filter_by(name = name).all()
            if not gs: return None
            g = gs[0]
            gene = Gene(name)
            for t in g.transcripts:
                if t.source.name not in self.source:
                    continue
                if t.refversion_id != self.refversion_id:
                    continue

                tpt = Transcript()
                tpt.chrm = t.chrm.name
                tpt.transcript_type = t.transcript_type.name
                tpt.strand = '-' if t.strand == 1 else '+'
                tpt.name = t.name
                tpt.beg = t.beg
                tpt.end = t.end
                tpt.cds_beg = t.cds_beg
                tpt.cds_end = t.cds_end
                tpt.source = t.source.name
                for ex in t.exons:
                    tpt.exons.append((ex.beg, ex.end))
                tpt.gene = gene
                tpt.exons.sort()
                gene.tpts.append(tpt)

            return gene
        elif self.name2gene:
            if name in self.name2gene:
                return self.name2gene[name]
            else:
                return None
        else:
            return None
            # raise Exception("No valid source of transcript definition")

    def _sql_get_transcripts(self, chrm, beg, end=None, flanking=0):

        if not end:
            end = beg
            
        beg = int(beg)
        end = int(end)
        chrm = normalize_chrm(chrm)
        dchrms = self.session.query(self.sqlmodel.DChromosome).filter_by(name=chrm).all()
        if dchrms:
            dchrm = dchrms[0]
        else:
            return []

        dtranscripts = self.session.query(self.sqlmodel.DTranscript).filter(
            self.sqlmodel.DTranscript.chrm_id == dchrm.id,
            self.sqlmodel.DTranscript.beg - flanking <= end,
            self.sqlmodel.DTranscript.end + flanking >= beg,
        ).all()

        # db_features = self.session.query(self.sqlmodel.Feature).filter(
        #     self.sqlmodel.Feature.chrm_id == db_chrm.id,
        #     self.sqlmodel.Feature.beg-flanking <= end,
        #     self.sqlmodel.Feature.end+flanking >= beg,
        # ).all()
        # print db_features
        tpts = []
        name2gene = {}
        for dtranscript in dtranscripts:

            if dtranscript.source.name not in self.source:
                continue

            tpt = Transcript()
            tpt.chrm = dtranscript.chrm.name
            tpt.strand = '-' if dtranscript.strand == 1 else '+'
            tpt.name = dtranscript.name
            tpt.beg = dtranscript.beg
            tpt.end = dtranscript.end
            tpt.cds_beg = dtranscript.cds_beg
            tpt.cds_end = dtranscript.cds_end
            tpt.transcript_type = dtranscript.transcript_type.name
            gene_name = dtranscript.gene.name
            if gene_name in name2gene:
                tpt.gene = name2gene[gene_name]
            else:
                tpt.gene = Gene(gene_name)
                name2gene[gene_name] = tpt.gene
            tpt.gene.tpts.append(tpt)
            tpt.source = dtranscript.source.name
            for ex in dtranscript.exons:
                tpt.exons.append((ex.beg, ex.end))
            tpt.exons.sort()
            tpts.append(tpt)

        return tpts

    def get_transcripts(self, chrm, beg, end=None, flanking=0):
        if self.session:
            return self._sql_get_transcripts(chrm, beg, end, flanking)
        elif self.thash:
            return [t for t in self.thash.get_transcripts(chrm, beg, end, flanking)]
        else:
            raise Exception("No valid source of transcript definition")

    def get_closest_transcripts_upstream(self, chrm, pos):
        if self.session:
            raise Exception("Not implemented.")
        elif self.thash:
            return self.thash.get_closest_transcripts_upstream(chrm, pos)
        else:
            raise Exception("No valid source of transcript definition")

    def get_closest_transcripts_downstream(self, chrm, pos):
        if self.session:
            raise Exception("Not implemented.")
        elif self.thash:
            return self.thash.get_closest_transcripts_downstream(chrm, pos)
        else:
            raise Exception("No valid source of transcript definition")

    def get_closest_transcripts(self, chrm, beg, end):
        """ closest transcripts upstream and downstream """
        return (self.get_closest_transcripts_upstream(chrm, beg),
                self.get_closest_transcripts_downstream(chrm, end))
