""" faidx python code adapted from Allen Yu 
http://www.allenyu.info/item/24-quickly-fetch-sequence-from-samtools-faidx-indexed-fasta-sequences.html """
import sys
class RefGenome:

    def __init__(self, fasta_file):       
        self.faidx = {}
         
        self.fasta_file=fasta_file
        try:
            self.fasta_handle=open(fasta_file)
        except IOError:
            print "Reference sequence doesn't exist"
 
        try:
            self.faidx_handle=open(fasta_file+".fai")
        except IOError:
            print "samtools faidx file doesn't exist for reference"   
        self.load_faidx()
         
    # Function to cache fasta index in dictionary
    # faidx format contains the following columns:
    ##.the name of the sequence
    ##.the length of the sequence
    ##.the offset of the first base in the file
    ##.the number of bases in each fasta line
    ##.the number of bytes in each fasta line   
    def load_faidx(self):

        for line in self.faidx_handle:
            line=line.strip()
            cols=line.split('\t')
            chrom = cols[0]
            slen,offset,blen,bytelen=[int(i) for i in cols[1:]]
            self.faidx[chrom]=(slen,offset,blen,bytelen)
     
    # Function to fetch sequence from an indexed fasta
    # *chrom--Chromosome name (str)
    # *start--Start position (1-based) (int)
    # *end--End position (1-based) (int)
    # *keepN--Keep ambiguous bases in sequence (boolean)       
    def fetch_sequence(self, chrom, start, end):

        # Fetch a sequence from start to end in 1-based coordinates
        seq=""
        
        if chrom not in self.faidx:
            if chrom.startswith('chr') and chrom[3:] in self.faidx:
                chrom = chrom[3:]
            elif 'chr'+chrom in self.faidx:
                chrom = 'chr'+chrom
            else:
                sys.stderr.write('Chromosome %s not found in reference\n' % chrom)
                return None

        slen,offset,blen,bytelen=self.faidx[chrom]
        start = start-1 #To 0-base
        # Sanity check of start and end position
        if start<0:
            raise IndexError('Sequence window out of bound--Chr: %s\tStart:%d\tEnd:%s' % (chrom,start+1,end))
        elif start>=end:
            raise ValueError('Start position %d is larger than end position %d' % (start+1,end))
        elif end>slen and start-(end-slen)>=0: #end is out of bound, adjust the window towards start
            end=slen
            start=start-(end-slen)
        elif end>slen:
            raise IndexError('Sequence window out of bound--Chr: %s\tStart:%d\tEnd:%s' % (chrom,start+1,end))
         
        self.fasta_handle.seek(offset+start/blen*bytelen+start%blen)
         
        while len(seq)<end-start:
            line=self.fasta_handle.readline()
            line=line[:-1] #Remove newline symbols
            seq=seq+line
            
        #chomp off extra bases
        return seq[:end-start]
        
    def __exit__(self, type, value, traceback):
        self.fasta_handle.close()
        self.faidx_handle.close()

def init_refgenome(r=None):
    global refgenome
    refgenome = RefGenome(r) if r else None
