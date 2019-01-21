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
