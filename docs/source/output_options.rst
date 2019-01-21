******************************
cDNA level annotation
******************************

Annotation from cDNA level is handled by the `canno` subcommand.

VCF-like output
#################


With `--gseq` transvar appends genomic sequence information as additional columns with pos, ref, alt following the VCF convention (i.e., indels are left-aligned)

.. code:: bash

   $ transvar canno -i 'MRE11A:c.592_593delGTinsTA' --ensembl --gseq

::

   MRE11A:c.592_593delGTinsTA	ENST00000323929 (protein_coding)	MRE11A	-
      chr11:g.94209521_94209522delinsTA/c.592_593delinsTA/p.V198*	inside_[cds_in_exon_7]
      CSQN=Missense;codon_cDNA=592-593-594;aliases=ENSP00000325863;source=Ensembl	c
      hr11	94209520	TAC	TTA
   MRE11A:c.592_593delGTinsTA	ENST00000323977 (protein_coding)	MRE11A	-
      chr11:g.94209521_94209522delinsTA/c.592_593delinsTA/p.V198*	inside_[cds_in_exon_7]
      CSQN=Missense;codon_cDNA=592-593-594;aliases=ENSP00000326094;source=Ensembl	c
      hr11	94209520	TAC	TTA
   MRE11A:c.592_593delGTinsTA	ENST00000393241 (protein_coding)	MRE11A	-
      chr11:g.94209521_94209522delinsTA/c.592_593delinsTA/p.V198*	inside_[cds_in_exon_7]
      CSQN=Missense;codon_cDNA=592-593-594;aliases=ENSP00000376933;source=Ensembl	c
      hr11	94209520	TAC	TTA
   MRE11A:c.592_593delGTinsTA	ENST00000540013 (protein_coding)	MRE11A	-
      chr11:g.94209521_94209522delinsTA/c.592_593delinsTA/p.V198*	inside_[cds_in_exon_7]
      CSQN=Missense;codon_cDNA=592-593-594;aliases=ENSP00000440986;source=Ensembl	c
      hr11	94209520	TAC	TTA


Another example of deletion

.. code:: bash

   $ transvar ganno -i "chr2:g.234183368_234183379del" --ccds --gseq --seqmax 200

::

   chr2:g.234183368_234183379del	CCDS2502.2 (protein_coding)	ATG16L1	+
      chr2:g.234183368_234183379delACTCATCCTGGT/c.841_852delACTCATCCTGGT/p.T281_G284delTHPG	inside_[cds_in_exon_8]
      CSQN=InFrameDeletion;left_align_gDNA=g.234183367_234183378delTACTCATCCTGG;una
      ligned_gDNA=g.234183368_234183379delACTCATCCTGGT;left_align_cDNA=c.840_851del
      TACTCATCCTGG;unalign_cDNA=c.841_852delACTCATCCTGGT;left_align_protein=p.T281_
      G284delTHPG;unalign_protein=p.T281_G284delTHPG;source=CCDS	chr2	234183366	ATA
      CTCATCCTGG	A
   chr2:g.234183368_234183379del	CCDS2503.2 (protein_coding)	ATG16L1	+
      chr2:g.234183368_234183379delACTCATCCTGGT/c.898_909delACTCATCCTGGT/p.T300_G303delTHPG	inside_[cds_in_exon_9]
      CSQN=InFrameDeletion;left_align_gDNA=g.234183367_234183378delTACTCATCCTGG;una
      ligned_gDNA=g.234183368_234183379delACTCATCCTGGT;left_align_cDNA=c.897_908del
      TACTCATCCTGG;unalign_cDNA=c.898_909delACTCATCCTGGT;left_align_protein=p.T300_
      G303delTHPG;unalign_protein=p.T300_G303delTHPG;source=CCDS	chr2	234183366	ATA
      CTCATCCTGG	A
   chr2:g.234183368_234183379del	CCDS54438.1 (protein_coding)	ATG16L1	+
      chr2:g.234183368_234183379delACTCATCCTGGT/c.409_420delACTCATCCTGGT/p.T137_G140delTHPG	inside_[cds_in_exon_5]
      CSQN=InFrameDeletion;left_align_gDNA=g.234183367_234183378delTACTCATCCTGG;una
      ligned_gDNA=g.234183368_234183379delACTCATCCTGGT;left_align_cDNA=c.408_419del
      TACTCATCCTGG;unalign_cDNA=c.409_420delACTCATCCTGGT;left_align_protein=p.T137_
      G140delTHPG;unalign_protein=p.T137_G140delTHPG;source=CCDS	chr2	234183366	ATA
      CTCATCCTGG	A

