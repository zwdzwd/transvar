******************************
Protein level annotation
******************************

Protein level inputs are handled by the `panno` subcommand.

Protein sites
####################

To use uniprot id as protein name, one must first download the uniprot id map by

.. code:: bash

   transvar config --download_idmap

Then one could use protein id instead of gene name by applying the `--uniprot` option to TransVar. For example,


.. code:: bash

   $ transvar panno --ccds -i 'Q5VUM1:47' --uniprot


::

   Q5VUM1:47	CCDS4972 (protein_coding)	C6ORF57	+
      chr6:g.71289191_71289193/c.139_141/p.47S	inside_[cds_in_exon_2]
      protein_sequence=S;cDNA_sequence=TCC;gDNA_sequence=TCC;source=CCDS

TransVar use a keyword extension `ref` in `Q5VUM1:p.47refS` to differentiate from the synonymous mutation `Q5VUM1:p.47S`. The former notation specifies that the reference protein sequence is `S` while the later specifies the target protein sequence is `S`.

Protein motif
####################

For example, one can find the genomic location of a DRY motif in protein P28222 by issuing the following command,

.. code:: bash

   $ transvar panno -i 'P28222:p.146_148refDRY' --uniprot --ccds

::

   P28222:p.146_148refDRY	CCDS4986 (protein_coding)	HTR1B	-
      chr6:g.78172677_78172685/c.436_444/p.D146_Y148	inside_[cds_in_exon_1]
      protein_sequence=DRY;cDNA_sequence=GACCGCTAC;gDNA_sequence=GTAGCGGTC;source=C
      CDS

One can also use wildcard `x` (lowercase) in the motif.

.. code:: bash

   $ transvar panno -i 'HTR1B:p.365_369refNPxxY' --ccds --seqmax 30

::

   HTR1B:p.365_369refNPxxY	CCDS4986 (protein_coding)	HTR1B	-
      chr6:g.78172014_78172028/c.1093_1107/p.N365_Y369	inside_[cds_in_exon_1]
      protein_sequence=NPIIY;cDNA_sequence=AACCCCATAATCTAT;gDNA_sequence=ATAGATTATG
      GGGTT;source=CCDS

Protein region
####################


.. code:: bash

   $ transvar panno --ccds -i 'ABCB11:p.200_400'

outputs

::

   ABCB11:p.200_400	CCDS46444 (protein_coding)	ABCB11	-
      chr2:g.169833195_169851872/c.598_1200/p.T200_K400	inside_[cds_in_exons_[6,7,8,9,10,11]]
      protein_sequence=TRF..DRK;cDNA_sequence=ACA..AAA;gDNA_sequence=TTT..TGT;sourc
      e=CCDS


Protein variants
#################################

Single amino acid substitution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mutation formats acceptable in TransVar are ```PIK3CA:p.E545K``` or without reference or alternative amino acid identity, e.g., ```PIK3CA:p.545K``` or ```PIK3CA:p.E545```. TransVar takes native HGVS format inputs and outputs. The reference amino acid is used to narrow the search scope of candidate transcripts. The alternative amino acid is used to infer nucleotide change which results in the amino acid.


.. code:: bash

   $ transvar panno -i PIK3CA:p.E545K --ensembl

outputs

::

   PIK3CA:p.E545K	ENST00000263967 (protein_coding)	PIK3CA	+
      chr3:g.178936091G>A/c.1633G>A/p.E545K	inside_[cds_in_exon_10]
      CSQN=Missense;reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_vari
      ants=chr3:g.178936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G
      >A);aliases=ENSP00000263967;source=Ensembl

One may encounter **ambiguous cases** where the multiple substitutions exist in explaining the amino acid change. For example,

.. code:: bash

   $ transvar panno -i ACSL4:p.R133R --ccds

outputs

::

   ACSL4:p.R133R	CCDS14548 (protein_coding)	ACSL4	-
      chrX:g.108926078G>T/c.399C>A/p.R133R	inside_[cds_in_exon_2]
      CSQN=Synonymous;reference_codon=CGC;candidate_codons=AGG,AGA,CGA,CGG,CGT;cand
      idate_snv_variants=chrX:g.108926078G>C,chrX:g.108926078G>A;candidate_mnv_vari
      ants=chrX:g.108926078_108926080delGCGinsCCT,chrX:g.108926078_108926080delGCGi
      nsTCT;source=CCDS

In those cases, TransVar prioritizes all the candidate base changes by minimizing the edit distance between the reference codon sequence and the target codon sequence. One of the optimal base changes is arbitrarily chosen as the default and all the candidates are included in the appended `CddMuts` entry.

Ambiguous amino acid code
++++++++++++++++++++++++++++

TransVar instantiates input of ambiguous amino acid code such as ('B', for "Asx", which stands for "Asp" or "Asn") to more specific amino acid. Even if the reference amino acid is a subset of the ambiguous alternative amino acid, TransVar assume a mutation on the nucleotide level (can still deduce synonymous mutations):

.. code:: bash

   $ transvar panno -i 'APC:p.D326B' --ccds

::

   APC:p.D326B	CCDS4107 (protein_coding)	APC	+
      chr5:g.112154705G>A/c.976G>A/p.D326N	inside_[cds_in_exon_9]
      CSQN=Missense;reference_codon=GAT;candidate_codons=AAC,AAT,GAC;candidate_snv_
      variants=chr5:g.112154707T>C;candidate_mnv_variants=chr5:g.112154705_11215470
      7delGATinsAAC;source=CCDS

Here input alternative amino acids is B (D or N). After TransVar processing, a 'N' is derived (though a D is equally likely, as shown in the candidates).

Insertion
^^^^^^^^^^^^


.. code:: bash

   $ transvar panno --ccds -i 'AATK:p.P1331_A1332insTP'

::

   AATK:p.P1331_A1332insTP	CCDS45807 (protein_coding)	AATK	-
      chr17:g.79093270_79093271insAGGTGT/c.3993_3994insACACCT/p.T1330_P1331dupTP	inside_[cds_in_exon_13]
      CSQN=InFrameInsertion;left_align_protein=p.A1326_P1327insPT;unalign_protein=p
      .T1330_P1331dupTP;left_align_gDNA=g.79093270_79093271insAGGTGT;unalign_gDNA=g
      .79093270_79093271insAGGTGT;left_align_cDNA=c.3993_3994insACACCT;unalign_cDNA
      =c.3993_3994insACACCT;16_CandidatesOmitted;source=CCDS

Deletion
^^^^^^^^^^

.. code:: bash

   $ transvar panno --ccds -i 'AADACL4:p.W263_I267delWRDAI'

::

   AADACL4:p.W263_I267delWRDAI	CCDS30590 (protein_coding)	AADACL4	+
      chr1:g.12726310_12726324del15/c.788_802del15/p.W263_I267delWRDAI	inside_[cds_in_exon_4]
      CSQN=InFrameDeletion;left_align_gDNA=g.12726308_12726322del15;unaligned_gDNA=
      g.12726309_12726323del15;left_align_cDNA=c.786_800del15;unalign_cDNA=c.787_80
      1del15;left_align_protein=p.W263_I267delWRDAI;unalign_protein=p.W263_I267delW
      RDAI;imprecise;source=CCDS

Block substitution
^^^^^^^^^^^^^^^^^^^^

.. code:: bash

   $ transvar panno --ccds -i 'ABCC3:p.Y556_V557delinsRRR'

::

   ABCC3:p.Y556_V557delinsRRR	CCDS32681 (protein_coding)	ABCC3	+
      chr17:g.48745254_48745259delinsAGGAGGAGG/c.1666_1671delinsAGGAGGAGG/p.Y556_V557delinsRRR	inside_[cds_in_exon_13]
      CSQN=MultiAAMissense;216_CandidatesOmitted;source=CCDS

Sometimes block substitution comes from in-frame deletion on the nucleotide level.

.. code:: bash

   $ transvar panno -i 'MAP2K1:p.F53_Q58delinsL' --ensembl

::

   MAP2K1:p.F53_Q58delinsL	ENST00000307102 (protein_coding)	MAP2K1	+
      chr15:g.66727443_66727457del15/c.159_173del15/p.F53_Q58delinsL	inside_[cds_in_exon_2]
      CSQN=MultiAAMissense;left_align_gDNA=g.66727443_66727457del15;unaligned_gDNA=
      g.66727443_66727457del15;left_align_cDNA=c.159_173del15;unalign_cDNA=c.159_17
      3del15;candidate_alternative_sequence=CTT/CTG/CTA/CTC/TTA/TTG;aliases=ENSP000
      00302486;source=Ensembl


Frame-shift variants
^^^^^^^^^^^^^^^^^^^^^^^

Frame-shift variants can be results of either insertion or deletion. In the cases where both are plausible the variants are prioritized by the length of the insertion/deletion. Mutations of smallest variants are given as the most likely inference. Other candidates are in given in the `candidates` field.

.. code:: bash

   $ transvar panno --refseq -i 'PTEN:p.T319fs*1' --max-candidates 2

::

   PTEN:p.T319fs*1	NM_000314 (protein_coding)	PTEN	+
      chr10:g.89720803_89720804insTA/c.954_955insTA/p.T319fs*1	inside_[cds_in_exon_8]
      CSQN=Frameshift;left_align_cDNA=c.954_955insTA;left_align_gDNA=g.89720803_897
      20804insTA;candidates=g.89720803_89720804insTG/c.954_955insTG/g.89720803_8972
      0804insTG/c.954_955insTG,g.89720804_89720807delACTT/c.955_958delACTT/g.897207
      99_89720802delTACT/c.950_953delTACT;1_CandidatesOmitted;dbxref=GeneID:5728,HG
      NC:9588,MIM:601728;aliases=NP_000305;source=RefSeq

In this example, both deletion `c.950_953delTACT` and insertion `c.954_955insTA` are possible. Both insertion involves fewer nucleotides and is chosen as the most likely inference. Deletion is given in the `candidates` tag.

The `candidates` field shows the right-aligned genomic, right-aligned cDNA, left-aligned genomic and left-aligned cDNA identifiers separated by `/`.

.. code:: bash

   $ transvar panno --ccds -i 'A1BG:p.G132fs*2' --max-candidates 1

::

   A1BG:p.G132fs*2	CCDS12976 (protein_coding)	A1BG	-
      chr19:g.58863868delC/c.395delG/p.G132fs*2	inside_[cds_in_exon_4]
      CSQN=Frameshift;left_align_cDNA=c.394delG;left_align_gDNA=g.58863867delC;cand
      idates=g.58863873delG/c.393delC/g.58863869delG/c.389delC;13_CandidatesOmitted
      ;source=CCDS

Frameshift variants can be difficult since there might be too many valid underlying nucleotide variants. 
Suppose we have a relatively long insertion,

.. code:: bash

   $ transvar ganno -i 'chr11:g.32417908_32417909insACCGTACA' --ccds

::

   chr11:g.32417908_32417909insACCGTACA	CCDS55750 (protein_coding)	WT1	-
      chr11:g.32417908_32417909insACCGTACA/c.456_457insTGTACGGT/p.A153Cfs*70	inside_[cds_in_exon_6]
      CSQN=Frameshift;left_align_gDNA=g.32417908_32417909insACCGTACA;unalign_gDNA=g
      .32417908_32417909insACCGTACA;left_align_cDNA=c.456_457insTGTACGGT;unalign_cD
      NA=c.456_457insTGTACGGT;source=CCDS
   chr11:g.32417908_32417909insACCGTACA	CCDS55751 (protein_coding)	WT1	-
      chr11:g.32417908_32417909insACCGTACA/c.507_508insTGTACGGT/p.A170Cfs*70	inside_[cds_in_exon_7]
      CSQN=Frameshift;left_align_gDNA=g.32417908_32417909insACCGTACA;unalign_gDNA=g
      .32417908_32417909insACCGTACA;left_align_cDNA=c.507_508insTGTACGGT;unalign_cD
      NA=c.507_508insTGTACGGT;source=CCDS
   chr11:g.32417908_32417909insACCGTACA	CCDS7878 (protein_coding)	WT1	-
      chr11:g.32417908_32417909insACCGTACA/c.1143_1144insTGTACGGT/p.A382Cfs*70	inside_[cds_in_exon_7]
      CSQN=Frameshift;left_align_gDNA=g.32417908_32417909insACCGTACA;unalign_gDNA=g
      .32417908_32417909insACCGTACA;left_align_cDNA=c.1143_1144insTGTACGGT;unalign_
      cDNA=c.1143_1144insTGTACGGT;source=CCDS
   chr11:g.32417908_32417909insACCGTACA	CCDS44561 (protein_coding)	WT1	-
      chr11:g.32417908_32417909insACCGTACA/c.1092_1093insTGTACGGT/p.A365Cfs*70	inside_[cds_in_exon_6]
      CSQN=Frameshift;left_align_gDNA=g.32417908_32417909insACCGTACA;unalign_gDNA=g
      .32417908_32417909insACCGTACA;left_align_cDNA=c.1092_1093insTGTACGGT;unalign_
      cDNA=c.1092_1093insTGTACGGT;source=CCDS
   chr11:g.32417908_32417909insACCGTACA	CCDS44562 (protein_coding)	WT1	-
      chr11:g.32417908_32417909insACCGTACA/c.1143_1144insTGTACGGT/p.A382Cfs*70	inside_[cds_in_exon_7]
      CSQN=Frameshift;left_align_gDNA=g.32417908_32417909insACCGTACA;unalign_gDNA=g
      .32417908_32417909insACCGTACA;left_align_cDNA=c.1143_1144insTGTACGGT;unalign_
      cDNA=c.1143_1144insTGTACGGT;source=CCDS

But now suppose we only know its protein identifier and forget about the original identifier. Using `panno`, we can get roughly how the original identifier look like:

.. code:: bash

   $ transvar panno -i 'WT1:p.A170Cfs*70' --ccds --max-candidates 2

would return more than 80 underlying variants. In this case the argument `--max-candidates` (default to 10) controls the maximum number of candidates output.

::

   WT1:p.A170Cfs*70	CCDS55751 (protein_coding)	WT1	-
      chr11:g.32417908_32417909insTTGGGGCA/c.507_508insTGCCCCAA/p.A170Cfs*70	inside_[cds_in_exons_[7,8,9]]
      CSQN=Frameshift;left_align_cDNA=c.507_508insTGCCCCAA;left_align_gDNA=g.324179
      08_32417909insTTGGGGCA;candidates=g.32417908_32417909insTTGNNNCA/c.507_508ins
      TGNNNCAA/g.32417908_32417909insTTGNNNCA/c.507_508insTGNNNCAA,g.32417908_32417
      909insGTGNNNCA/c.507_508insTGNNNCAC/g.32417908_32417909insGTGNNNCA/c.507_508i
      nsTGNNNCAC;80_CandidatesOmitted;source=CCDS

Sometimes the alternative amino acid can be missing

.. code:: bash

   $ transvar panno -i ADAMTSL1:p.I396fs*30 --ccds --max-candidates 2

::

   ADAMTSL1:p.I396fs*30	CCDS6485 (protein_coding)	ADAMTSL1	+
      chr9:g.18680360_18680361insG/c.1187_1188insG/p.I396fs*30	inside_[cds_in_exon_11]
      CSQN=Frameshift;left_align_cDNA=c.1187_1188insG;left_align_gDNA=g.18680360_18
      680361insG;candidates=g.18680359dupA/c.1186dupA/g.18680358_18680359insA/c.118
      5_1186insA,g.18680359_18680360insC/c.1186_1187insC/g.18680359_18680360insC/c.
      1186_1187insC;11_CandidatesOmitted;source=CCDS
   ADAMTSL1:p.I396fs*30	CCDS47954 (protein_coding)	ADAMTSL1	+
      chr9:g.18680360_18680361insG/c.1187_1188insG/p.I396fs*30	inside_[cds_in_exon_11]
      CSQN=Frameshift;left_align_cDNA=c.1187_1188insG;left_align_gDNA=g.18680360_18
      680361insG;candidates=g.18680359dupA/c.1186dupA/g.18680358_18680359insA/c.118
      5_1186insA,g.18680359_18680360insC/c.1186_1187insC/g.18680359_18680360insC/c.
      1186_1187insC;11_CandidatesOmitted;source=CCDS


TransVar can also take protein identifiers such as  as input. For example,

.. code:: bash

   $ transvar panno --refseq -i 'NP_006266.2:p.G240Afs*50'


::

   NP_006266.2:p.G240Afs*50	NM_006275 (protein_coding)	SRSF6	+
      chr20:g.42089387delG/c.719delG/p.G240Afs*50	inside_[cds_in_exon_6]
      CSQN=Frameshift;left_align_cDNA=c.718delG;left_align_gDNA=g.42089386delG;cand
      idates=g.42089385delA/c.717delA/g.42089382delA/c.714delA;dbxref=GeneID:6431,H
      GNC:10788,HPRD:09054,MIM:601944;aliases=NP_006266;source=RefSeq

The output gives the exact details of the mutation on the DNA levels, properly right-aligned. The `candidates` fields also include other equally-likely mutation identifiers. `candidates` have the format `[right-align-gDNA]/[right-align-cDNA]/[left-align-gDNA]/[left-align-cDNA]` for each hit and `,` separation between hits. 

Similar applies when the underlying mutation is an insertion. TransVar can infer insertion sequence of under 3 base pairs long. For example,

.. code:: bash

   $ transvar panno -i 'AASS:p.I355Mfs*10' --ccds --max-candidates 1

::

   AASS:p.I355Mfs*10	CCDS5783 (protein_coding)	AASS	-
      chr7:g.121753753_121753754insTC/c.1064_1065insGA/p.I355Mfs*10	inside_[cds_in_exon_9]
      CSQN=Frameshift;left_align_cDNA=c.1064_1065insGA;left_align_gDNA=g.121753753_
      121753754insTC;candidates=g.121753753_121753754insGC/c.1064_1065insGC/g.12175
      3753_121753754insGC/c.1064_1065insGC;3_CandidatesOmitted;source=CCDS

When the alternative becomes a stop codon, frameshift mutation becomes a nonsense mutation:

.. code:: bash

   $ transvar panno -i 'APC:p.I1557*fs*3' --ccds

returns a nonsense mutation

::

   APC:p.I1557*fs*3	CCDS4107 (protein_coding)	APC	+
      chr5:g.112175960_112175962delATTinsTAA/c.4669_4671delATTinsTAA/p.I1557*	inside_[cds_in_exon_15]
      CSQN=Nonsense;reference_codon=ATT;candidate_codons=TAA,TAG,TGA;candidate_mnv_
      variants=chr5:g.112175960_112175962delATTinsTAG,chr5:g.112175960_112175962del
      ATTinsTGA;source=CCDS


Whole transcript
###################

TransVar provides an easy way to investigate a whole transcript by supplying the gene id.

.. code:: bash

   $ transvar panno -i 'Dnmt3a' --refseq

outputs the basic information of transcripts of the protein, in an intuitive way,

::

   Dnmt3a	XM_005264176 (protein_coding)	DNMT3A	-
      chr2:g.25451421_25537541/c.1_2739/p.M1_*913	whole_transcript
      promoter=chr2:25537541_25538541;#exons=23;cds=chr2:25457148_25536853
   Dnmt3a	XM_005264175 (protein_coding)	DNMT3A	-
      chr2:g.25451421_25537354/c.1_2739/p.M1_*913	whole_transcript
      promoter=chr2:25537354_25538354;#exons=23;cds=chr2:25457148_25536853
   Dnmt3a	XM_005264177 (protein_coding)	DNMT3A	-
      chr2:g.25451421_25475145/c.1_2070/p.M1_*690	whole_transcript
      promoter=chr2:25475145_25476145;#exons=18;cds=chr2:25457148_25471091
   Dnmt3a	NM_175629 (protein_coding)	DNMT3A	-
      chr2:g.25455830_25565459/c.1_2739/p.M1_*913	whole_transcript
      promoter=chr2:25565459_25566459;#exons=23;cds=chr2:25457148_25536853
   Dnmt3a	NM_022552 (protein_coding)	DNMT3A	-
      chr2:g.25455830_25564784/c.1_2739/p.M1_*913	whole_transcript
      promoter=chr2:25564784_25565784;#exons=23;cds=chr2:25457148_25536853
   Dnmt3a	NM_153759 (protein_coding)	DNMT3A	-
      chr2:g.25455830_25475184/c.1_2172/p.M1_*724	whole_transcript
      promoter=chr2:25475184_25476184;#exons=19;cds=chr2:25457148_25475066
   Dnmt3a	NM_175630 (protein_coding)	DNMT3A	-
      chr2:g.25504321_25565459/c.1_501/p.M1_*167	whole_transcript
      promoter=chr2:25565459_25566459;#exons=4;cds=chr2:25505257_25536853


Search alternative codon identifiers
########################################

An identifier is regarded as an alternative if the underlying codon overlap with the one from the original identifier.
Example: to search alternative identifiers of CDKN2A.p.58 (without knowing reference allele),

.. code:: bash

   $ transvar codonsearch --ccds -i CDKN2A:p.58

::

   origin_id	alt_id	chrm	codon1	codon2	transcripts_choice
   CDKN2A:p.58	CDKN2A.p.73	chr9	21971184-21971185-21971186
      21971182-21971183-21971184	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
   CDKN2A:p.58	CDKN2A.p.72	chr9	21971184-21971185-21971186
      21971185-21971186-21971187	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]

The pair of transcript id listed corresponds to the transcripts based on which, the original and alternative identifiers are defined. Multiple pairs of transcript definitions are appended following a `,`.

Example: to search alternative identifiers of DHODH:G152R (knowing reference allele `G`, alternative allele here will be ignored),

.. code:: bash

   $ transvar codonsearch -i DHODH:G152R --refseq

outputs

::

   origin_id	alt_id	chrm	codon1	codon2	transcripts_choice
   DHODH:G152R	DHODH.p.G124	chr16	72050942-72050943-72050944
      72050942-72050943-72050944	NM_001361[RefSeq]/XM_005255827[RefSeq]
   DHODH:G152R	DHODH.p.G16	chr16	72050942-72050943-72050944
      72050942-72050943-72050944	NM_001361[RefSeq]/XM_005255828[RefSeq]
   DHODH:G152R	DHODH.p.G9	chr16	72050942-72050943-72050944
      72050942-72050943-72050944	NM_001361[RefSeq]/XM_005255829[RefSeq]

TransVar outputs genomic positions of codons based on original transcript (4th column in the output) and alternative transcript (5th column in the output). The potential transcript usages are also appended.

Example: to run `transvar codonsearch` to **batch process** a list of mutation identifiers.

.. code:: bash

   $ transvar codonsearch -l example/input_table2 --ccds -m 1 -o 1

Example input table

::

   origin_id	alt_id	chrm	codon1	codon2	transcripts_choice
   CDKN2A:p.61	CDKN2A.p.76	chr9	21971175-21971176-21971177
      21971173-21971174-21971175	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
   CDKN2A:p.61	CDKN2A.p.75	chr9	21971175-21971176-21971177
      21971176-21971177-21971178	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
   CDKN2A:p.69	CDKN2A.p.84	chr9	21971151-21971152-21971153
      21971149-21971150-21971151	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
   CDKN2A:p.69	CDKN2A.p.83	chr9	21971151-21971152-21971153
      21971152-21971153-21971154	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
   CDKN2A:p.69	CDKN2A.p.55	chr9	21971194-21971195-21971196
      21971193-21971194-21971195	CCDS6511[CCDS]/CCDS6510[CCDS],CCDS6511[CCDS]/CCDS56565[CCDS]
   CDKN2A:p.69	CDKN2A.p.54	chr9	21971194-21971195-21971196
      21971196-21971197-21971198	CCDS6511[CCDS]/CCDS6510[CCDS],CCDS6511[CCDS]/CCDS56565[CCDS]
   ERBB2:p.755	ERBB2.p.725	chr17	37880219-37880220-37880221
      37880219-37880220-37880221	CCDS32642[CCDS]/CCDS45667[CCDS]
   ERBB2:p.755	ERBB2.p.785	chr17	37881024-37881025-37881026
      37881024-37881025-37881026	CCDS45667[CCDS]/CCDS32642[CCDS]

outputs

::

   origin_id	alt_id	chrm	codon1
      codon2	transcripts_choice
   CDKN2A:p.61	CDKN2A.p.76	chr9	21971175-21971176-21971177
      21971173-21971174-21971175	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
   CDKN2A:p.61	CDKN2A.p.75	chr9	21971175-21971176-21971177
      21971176-21971177-21971178	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
   CDKN2A:p.69	CDKN2A.p.54	chr9	21971194-21971195-21971196
      21971196-21971197-21971198	CCDS6511[CCDS]/CCDS6510[CCDS],CCDS6511[CCDS]/CCDS56565[CCDS]
   CDKN2A:p.69	CDKN2A.p.55	chr9	21971194-21971195-21971196
      21971193-21971194-21971195	CCDS6511[CCDS]/CCDS6510[CCDS],CCDS6511[CCDS]/CCDS56565[CCDS]
   CDKN2A:p.69	CDKN2A.p.83	chr9	21971151-21971152-21971153
      21971152-21971153-21971154	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
   CDKN2A:p.69	CDKN2A.p.84	chr9	21971151-21971152-21971153
      21971149-21971150-21971151	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
   ERBB2:p.755	ERBB2.p.785	chr17	37881024-37881025-37881026
      37881024-37881025-37881026	CCDS45667[CCDS]/CCDS32642[CCDS]
   ERBB2:p.755	ERBB2.p.725	chr17	37880219-37880220-37880221
      37880219-37880220-37880221	CCDS32642[CCDS]/CCDS45667[CCDS]

The third column indicates the potential transcript usage for the alternative identifier. Each transcript usage is denoted by <listing transcript>/<actual transcript>. Different potential choices are separated by ','.

Infer potential codon identity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example: to check if MET.p1010 and MET.p992 may be refering to one mutation due to different usage of transcripts,

.. code:: bash

   $ transvar codonsearch --refseq -i MET:p.1010

gives

::

   origin_id	alt_id	chrm	codon1	codon2	transcripts_choice
   MET:p.1010	MET.p.991	chr7	116411932-116411933-116411934
      116411932-116411933-116411934	XM_005250353[RefSeq]/NM_001127500[RefSeq]
   MET:p.1010	MET.p.973	chr7	116411932-116411933-116411934
      116411932-116411933-116411934	XM_005250353[RefSeq]/NM_000245[RefSeq]
   MET:p.1010	MET.p.543	chr7	116411932-116411933-116411934
      116411932-116411933-116411934	XM_005250353[RefSeq]/XM_005250354[RefSeq]
   MET:p.1010	MET.p.1029	chr7	116411989-116411990-116411991
      116411989-116411990-116411991	NM_001127500[RefSeq]/XM_005250353[RefSeq]
   MET:p.1010	MET.p.992	chr7	116411989-116411990-116411991
      116411989-116411990-116411991	NM_001127500[RefSeq]/NM_000245[RefSeq]
   MET:p.1010	MET.p.562	chr7	116411989-116411990-116411991
      116411989-116411990-116411991	NM_001127500[RefSeq]/XM_005250354[RefSeq]
   MET:p.1010	MET.p.1047	chr7	116412043-116414935-116414936
      116412043-116414935-116414936	NM_000245[RefSeq]/XM_005250353[RefSeq]
   MET:p.1010	MET.p.1028	chr7	116412043-116414935-116414936
      116412043-116414935-116414936	NM_000245[RefSeq]/NM_001127500[RefSeq]
   MET:p.1010	MET.p.580	chr7	116412043-116414935-116414936
      116412043-116414935-116414936	NM_000245[RefSeq]/XM_005250354[RefSeq]

Since MET.p.992 is in the list, the two identifiers might be due to the same genomic mutation.
