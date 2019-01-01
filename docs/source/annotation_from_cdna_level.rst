******************************
cDNA level annotation
******************************

Annotation from cDNA level is handled by the `canno` subcommand.

cDNA region
#############


.. code:: bash

   $ transvar canno --ccds -i 'ABCB11:c.1198-8_1202'

outputs

::

   ABCB11:c.1198-8_1202	CCDS46444 (protein_coding)	ABCB11	-
      chr2:g.169833193_169833205GGTTTCTGGAGTG/c.1198-8_1202CACTCCAGAAACC/p.400_401KP	from_[cds_in_exon_11]_to_[intron_between_exon_10_and_11]
      C2=acceptor_splice_site_on_exon_11_at_chr2:169833198_included;source=CCDS

cDNA variant
##############

Single Nucleotide Variation (SNV)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

TransVar infers nucleotide mutation through ``PIK3CA:c.1633G>A``. Note that nucleotide identity follows the natural sequence, i.e., if transcript is interpreted on the reverse-complementary strand, the base at the site needs to be reverse-complemented too.

.. code:: bash

   $ transvar canno --ccds -i 'PIK3CA:c.1633G>A'

outputs

::

   PIK3CA:c.1633G>A	CCDS43171 (protein_coding)	PIK3CA	+
      chr3:g.178936091G>A/c.1633G>A/p.E545K	inside_[cds_in_exon_9]
      CSQN=Missense;dbsnp=rs104886003(chr3:178936091G>A);reference_codon=GAG;altern
      ative_codon=AAG;source=CCDS

The SNV can be in the intronic region, e.g.,

.. code:: bash

   $ transvar canno --ccds -i 'ABCB11:c.1198-8C>A'

outputs

::

   ABCB11:c.1198-8C>A	CCDS46444 (protein_coding)	ABCB11	-
      chr2:g.169833205G>T/c.1198-8C>A/.	inside_[intron_between_exon_10_and_11]
      CSQN=IntronicSNV;source=CCDS

Or in the 5'-UTR region, e.g.,

.. code:: bash

   $ transvar canno -i 'KCNJ11:c.-134G>T' --ensembl

::

   KCNJ11:c.-134G>T	ENST00000339994 (protein_coding)	KCNJ11	-
      chr11:g.17409772C>A/c.1-134G>T/.	inside_[5-UTR;noncoding_exon_1]
      CSQN=5-UTRSNV;dbsnp=rs387906398(chr11:17409772C>A);aliases=ENSP00000345708;so
      urce=Ensembl

Or in the 3'-UTR region, e.g.,

.. code:: bash

   $ transvar canno -i 'MSH2:c.*95C>T' --refseq

::

   MSH2:c.*95C>T	NM_000251 (protein_coding)	MSH2	+
      chr2:g.47710183C>T/c.*95C>T/.	inside_[3-UTR;noncoding_exon_16]
      CSQN=3-UTRSNV;dbsnp=rs587779062(chr2:47710183C>T);dbxref=GeneID:4436,HGNC:732
      5,HPRD:00389,MIM:609309;aliases=NP_000242;source=RefSeq
   MSH2:c.*95C>T	NM_001258281 (protein_coding)	MSH2	+
      chr2:g.47710183C>T/c.*95C>T/.	inside_[3-UTR;noncoding_exon_17]
      CSQN=3-UTRSNV;dbsnp=rs587779062(chr2:47710183C>T);dbxref=GeneID:4436,HGNC:732
      5,HPRD:00389,MIM:609309;aliases=NP_001245210;source=RefSeq
   MSH2:c.*95C>T	XM_005264333 (protein_coding)	MSH2	+
      chr2:g.47710183C>T/c.*95C>T/.	inside_[3-UTR;noncoding_exon_15]
      CSQN=3-UTRSNV;dbsnp=rs587779062(chr2:47710183C>T);dbxref=GeneID:4436,HGNC:732
      5,HPRD:00389,MIM:609309;aliases=XP_005264390;source=RefSeq

      
insertion
^^^^^^^^^^^^

An insertion may result in: 1) a pure insertion of amino acids; 2) a block substitution of amino acids, when insertion occur after 1st or 2nd base in a codon; or 3) a frame-shift. Following HGVS nomenclature, TransVar labels the first different amino acid and the length of the peptide util stop codon, assuming no change in the splicing.

Example: to annotate an **in-frame, in-phase insertion**,

.. code:: bash

   $ transvar canno --ccds -i 'ACIN1:c.1932_1933insATTCAC'

::

   ACIN1:c.1932_1933insATTCAC	CCDS9587 (protein_coding)	ACIN1	-
      chr14:g.23548785_23548786insGTGAAT/c.1932_1933insATTCAC/p.R644_S645insIH	inside_[cds_in_exon_6]
      CSQN=InFrameInsertion;left_align_gDNA=g.23548785_23548786insGTGAAT;unalign_gD
      NA=g.23548785_23548786insGTGAAT;left_align_cDNA=c.1932_1933insATTCAC;unalign_
      cDNA=c.1932_1933insATTCAC;left_align_protein=p.R644_S645insIH;unalign_protein
      =p.R644_S645insIH;phase=0;source=CCDS
   ACIN1:c.1932_1933insATTCAC	CCDS53889 (protein_coding)	ACIN1	-
      chr14:g.23548157_23548158insGTGAAT/c.1932_1933insATTCAC/p.P644_V645insIH	inside_[cds_in_exon_6]
      CSQN=InFrameInsertion;left_align_gDNA=g.23548157_23548158insGTGAAT;unalign_gD
      NA=g.23548157_23548158insGTGAAT;left_align_cDNA=c.1932_1933insATTCAC;unalign_
      cDNA=c.1932_1933insATTCAC;left_align_protein=p.P644_V645insIH;unalign_protein
      =p.P644_V645insIH;phase=0;source=CCDS
   ACIN1:c.1932_1933insATTCAC	CCDS55905 (protein_coding)	ACIN1	-
      chr14:g.23548785_23548786insGTGAAT/c.1932_1933insATTCAC/p.R644_S645insIH	inside_[cds_in_exon_6]
      CSQN=InFrameInsertion;left_align_gDNA=g.23548785_23548786insGTGAAT;unalign_gD
      NA=g.23548785_23548786insGTGAAT;left_align_cDNA=c.1932_1933insATTCAC;unalign_
      cDNA=c.1932_1933insATTCAC;left_align_protein=p.R644_S645insIH;unalign_protein
      =p.R644_S645insIH;phase=0;source=CCDS

``Phase = 0,1,2`` indicates whether the insertion happen after the 3rd, 1st or 2nd base of a codon, respectively. An insertion *in phase* refers to one with ``Phase=0``.

Example: to annotate an **out-of-phase, in-frame insertion**,

.. code:: bash

   $ transvar canno --ccds -i 'ACIN1:c.1930_1931insATTCAC'

::

   ACIN1:c.1930_1931insATTCAC	CCDS9587 (protein_coding)	ACIN1	-
      chr14:g.23548792_23548793insTGTGAA/c.1930_1931insATTCAC/p.S643_R644insHS	inside_[cds_in_exon_6]
      CSQN=InFrameInsertion;left_align_gDNA=g.23548787_23548788insGTGAAT;unalign_gD
      NA=g.23548787_23548788insGTGAAT;left_align_cDNA=c.1925_1926insTTCACA;unalign_
      cDNA=c.1930_1931insATTCAC;left_align_protein=p.R642_S643insSH;unalign_protein
      =p.S643_R644insHS;phase=1;source=CCDS
   ACIN1:c.1930_1931insATTCAC	CCDS53889 (protein_coding)	ACIN1	-
      chr14:g.23548162_23548163insAATGTG/c.1930_1931insATTCAC/p.P643_P644insHS	inside_[cds_in_exon_6]
      CSQN=InFrameInsertion;left_align_gDNA=g.23548159_23548160insGTGAAT;unalign_gD
      NA=g.23548159_23548160insGTGAAT;left_align_cDNA=c.1927_1928insCACATT;unalign_
      cDNA=c.1930_1931insATTCAC;left_align_protein=p.P643_P644insHS;unalign_protein
      =p.P643_P644insHS;phase=1;source=CCDS
   ACIN1:c.1930_1931insATTCAC	CCDS55905 (protein_coding)	ACIN1	-
      chr14:g.23548792_23548793insTGTGAA/c.1930_1931insATTCAC/p.S643_R644insHS	inside_[cds_in_exon_6]
      CSQN=InFrameInsertion;left_align_gDNA=g.23548787_23548788insGTGAAT;unalign_gD
      NA=g.23548787_23548788insGTGAAT;left_align_cDNA=c.1925_1926insTTCACA;unalign_
      cDNA=c.1930_1931insATTCAC;left_align_protein=p.R642_S643insSH;unalign_protein
      =p.S643_R644insHS;phase=1;source=CCDS

Reverse annotation can result in different identifiers after left/right alignments, e.g., 

.. code:: bash

   $ transvar canno --ccds -i 'AATK:c.3976_3977insCGCCCA'

results in

::

   AATK:c.3976_3977insCGCCCA	CCDS45807 (protein_coding)	AATK	-
      chr17:g.79093282_79093287dupTGGGCG/c.3988_3993dupACGCCC/p.T1330_P1331dupTP	inside_[cds_in_exon_13]
      CSQN=InFrameInsertion;left_align_gDNA=g.79093270_79093271insGGGCGT;unalign_gD
      NA=g.79093282_79093287dupTGGGCG;left_align_cDNA=c.3976_3977insCGCCCA;unalign_
      cDNA=c.3976_3977insCGCCCA;left_align_protein=p.A1326_P1327insPT;unalign_prote
      in=p.A1326_P1327insPT;phase=1;source=CCDS

Note how insertion switch to duplication when 5'flanking is identical. This conforms to HGVS recommendation to replace insertion notation with duplication when possible.

Example: to annotate a **frame-shift insertion**, frameshift mutations have not alternative alignments. Hence only cDNA and gDNA have left alignment and unalignment reports.

.. code:: bash

   $ transvar canno --ccds -i 'AAAS:c.1225_1226insG'

results in

::

   AAAS:c.1225_1226insG	CCDS8856 (protein_coding)	AAAS	-
      chr12:g.53702093dupC/c.1225dupG/p.E409Gfs*17	inside_[cds_in_exon_13]
      CSQN=Frameshift;left_align_gDNA=g.53702089_53702090insC;unalign_gDNA=g.537020
      89_53702090insC;left_align_cDNA=c.1221_1222insG;unalign_cDNA=c.1225dupG;sourc
      e=CCDS
   AAAS:c.1225_1226insG	CCDS53797 (protein_coding)	AAAS	-
      chr12:g.53701842_53701843insC/c.1225_1226insG/p.L409Rfs*54	inside_[cds_in_exon_13]
      CSQN=Frameshift;left_align_gDNA=g.53701842_53701843insC;unalign_gDNA=g.537018
      42_53701843insC;left_align_cDNA=c.1225_1226insG;unalign_cDNA=c.1225_1226insG;
      source=CCDS

Example: to annotate an **intronic insertion**,

.. code:: bash

   $ transvar canno --ccds -i 'ADAM33:c.991-3_991-2insC'

outputs

::

   ADAM33:c.991-3_991-2insC	CCDS13058 (protein_coding)	ADAM33	-
      chr20:g.3654151dupG/c.991-3dupC/.	inside_[intron_between_exon_10_and_11]
      CSQN=IntronicInsertion;left_align_gDNA=g.3654145_3654146insG;unalign_gDNA=g.3
      654145_3654146insG;left_align_cDNA=c.991-9_991-8insC;unalign_cDNA=c.991-3dupC
      ;source=CCDS

In the case of intronic insertions, amino acid identifier is not applicable, represented in a `.`. But cDNA and gDNA identifier are right-aligned according to their natural order, respecting HGVS nomenclature.

Insertion could occur to *splice sites*. TransVar identifies such cases and report splice site and repress translation of protein change.

.. code:: bash

   $ transvar canno --ccds -i 'ADAM33:c.991_992insC'

results in

::

   ADAM33:c.991_992insC	CCDS13058 (protein_coding)	ADAM33	-
      chr20:g.3654142_3654143insG/c.991_992insC/.	inside_[cds_in_exon_11]
      CSQN=SpliceAcceptorInsertion;left_align_gDNA=g.3654142_3654143insG;unalign_gD
      NA=g.3654142_3654143insG;left_align_cDNA=c.991_992insC;unalign_cDNA=c.991_992
      insC;C2=acceptor_splice_site_on_exon_11_at_chr20:3654144_affected;source=CCDS

deletion
^^^^^^^^^^

Similar to insertions, deletion can be in-frame or frame-shift. The consequence of deletion to amino acid sequence may appear a simple deletion or a block substitution (in the case where in-frame deletion is out of phase, i.e., partially delete codons).

Example: to annotate an **in-frame deletion**,

.. code:: bash

   $ transvar canno --ccds -i 'A4GNT:c.694_696delTTG'

::

   A4GNT:c.694_696delTTG	CCDS3097 (protein_coding)	A4GNT	-
      chr3:g.137843435_137843437delACA/c.694_696delTTG/p.L232delL	inside_[cds_in_exon_2]
      CSQN=InFrameDeletion;left_align_gDNA=g.137843433_137843435delCAA;unaligned_gD
      NA=g.137843433_137843435delCAA;left_align_cDNA=c.692_694delTGT;unalign_cDNA=c
      .694_696delTTG;left_align_protein=p.L232delL;unalign_protein=p.L232delL;sourc
      e=CCDS

Example: to annotate a **in-frame, out-of-phase deletion**,

.. code:: bash

   $ transvar canno --ccds -i 'ABHD15:c.431_433delGTG'

::

   ABHD15:c.431_433delGTG	CCDS32602 (protein_coding)	ABHD15	-
      chr17:g.27893552_27893554delCAC/c.431_433delGTG/p.C144_V145delinsF	inside_[cds_in_exon_1]
      CSQN=MultiAAMissense;left_align_gDNA=g.27893552_27893554delCAC;unaligned_gDNA
      =g.27893552_27893554delCAC;left_align_cDNA=c.431_433delGTG;unalign_cDNA=c.431
      _433delGTG;source=CCDS

Example: to annotate a **frame-shift deletion**,

.. code:: bash

   $ transvar canno --ccds -i 'AADACL3:c.374delG'

::

   AADACL3:c.374delG	CCDS41252 (protein_coding)	AADACL3	+
      chr1:g.12785494delG/c.374delG/p.C125Ffs*17	inside_[cds_in_exon_3]
      CSQN=Frameshift;left_align_gDNA=g.12785494delG;unaligned_gDNA=g.12785494delG;
      left_align_cDNA=c.374delG;unalign_cDNA=c.374delG;source=CCDS

Example: to annotate a **deletion that span from intronic to coding region**, protein prediction is suppressed due to loss of splice site.

.. code:: bash

   $ transvar canno --ccds -i 'ABCB11:c.1198-8_1199delcactccagAA'

::

   ABCB11:c.1198-8_1199delcactccagAA	CCDS46444 (protein_coding)	ABCB11	-
      chr2:g.169833196_169833205delTTCTGGAGTG/c.1198-8_1199delCACTCCAGAA/.	from_[cds_in_exon_11]_to_[intron_between_exon_10_and_11]
      CSQN=SpliceAcceptorDeletion;left_align_gDNA=g.169833196_169833205delTTCTGGAGT
      G;unaligned_gDNA=g.169833196_169833205delTTCTGGAGTG;left_align_cDNA=c.1198-8_
      1199delCACTCCAGAA;unalign_cDNA=c.1198-8_1199delCACTCCAGAA;C2=acceptor_splice_
      site_on_exon_11_at_chr2:169833198_lost;source=CCDS


block substitution
^^^^^^^^^^^^^^^^^^^^

Example: to annotate a block substitution in **coding region**,

.. code:: bash

   $ transvar canno --ccds -i 'A1CF:c.508_509delinsTT'

::

   A1CF:c.508_509delinsTT	CCDS7241 (protein_coding)	A1CF	-
      chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
      CSQN=Missense;codon_cDNA=508-509-510;source=CCDS
   A1CF:c.508_509delinsTT	CCDS7242 (protein_coding)	A1CF	-
      chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
      CSQN=Missense;codon_cDNA=508-509-510;source=CCDS
   A1CF:c.508_509delinsTT	CCDS7243 (protein_coding)	A1CF	-
      chr10:g.52595953_52595954delinsAA/c.508_509delinsTT/p.G170F	inside_[cds_in_exon_4]
      CSQN=Missense;codon_cDNA=508-509-510;source=CCDS

Block substitution does not necessarily results in block substitution in amino acid. For example, the following substitution results in a deletion, where protein alternative alignment should be reported.

.. code:: bash

   $ transvar canno --ccds -i 'CSRNP1:c.1212_1224delinsGGAGGAGGAA'

::

   CSRNP1:c.1212_1224delinsGGAGGAGGAA	CCDS2682 (protein_coding)	CSRNP1	-
      chr3:g.39185092_39185104delinsTTCCTCCTCC/c.1212_1224delinsGGAGGAGGAA/p.E411delE	inside_[cds_in_exon_4]
      CSQN=InFrameDeletion;begin_codon_cDNA=1210-1211-1212;end_codon_cDNA=1222-1223
      -1224;left_align_protein=p.E405delE;unalign_protein=p.E408delE;source=CCDS

Likewise, block substitution could occur to **intronic region**,

.. code:: bash

   $ transvar canno --ccds -i 'A1CF:c.1460+2_1460+3delinsCC'


::

   A1CF:c.1460+2_1460+3delinsCC	CCDS7241 (protein_coding)	A1CF	-
      chr10:g.52570797_52570798delinsGG/c.1460+2_1460+3delinsCC/.	inside_[intron_between_exon_9_and_10]
      CSQN=IntronicBlockSubstitution;source=CCDS

When block substitution occurs **across splice site**, TransVar put a tag in the info fields and does not predict amino acid change.

.. code:: bash

   $ transvar canno --ccds -i 'A1CF:c.1459_1460+3delinsCC'


::

   A1CF:c.1459_1460+3delinsCC	CCDS7241 (protein_coding)	A1CF	-
      chr10:g.52570797_52570801delinsGG/c.1459_1460+3delinsCC/.	from_[intron_between_exon_9_and_10]_to_[cds_in_exon_9]
      CSQN=SpliceDonorBlockSubstitution;C2=donor_splice_site_on_exon_9_at_chr10:525
      70799_lost;source=CCDS


With `--gseq` transvar appends genomic sequence information as additional columns

.. code:: bash

   $ transvar canno -i 'MRE11A:c.592_593delGTinsTA' --ensembl --gseq

::

   MRE11A:c.592_593delGTinsTA	ENST00000323929 (protein_coding)	MRE11A	-
      chr11:g.94209521_94209522delinsTA/c.592_593delinsTA/p.V198*	inside_[cds_in_exon_7]
      CSQN=Missense;codon_cDNA=592-593-594;aliases=ENSP00000325863;source=Ensembl	c
      hr11	94209521	94209522	AC	TA
   MRE11A:c.592_593delGTinsTA	ENST00000323977 (protein_coding)	MRE11A	-
      chr11:g.94209521_94209522delinsTA/c.592_593delinsTA/p.V198*	inside_[cds_in_exon_7]
      CSQN=Missense;codon_cDNA=592-593-594;aliases=ENSP00000326094;source=Ensembl	c
      hr11	94209521	94209522	AC	TA
   MRE11A:c.592_593delGTinsTA	ENST00000393241 (protein_coding)	MRE11A	-
      chr11:g.94209521_94209522delinsTA/c.592_593delinsTA/p.V198*	inside_[cds_in_exon_7]
      CSQN=Missense;codon_cDNA=592-593-594;aliases=ENSP00000376933;source=Ensembl	c
      hr11	94209521	94209522	AC	TA
   MRE11A:c.592_593delGTinsTA	ENST00000540013 (protein_coding)	MRE11A	-
      chr11:g.94209521_94209522delinsTA/c.592_593delinsTA/p.V198*	inside_[cds_in_exon_7]
      CSQN=Missense;codon_cDNA=592-593-594;aliases=ENSP00000440986;source=Ensembl	c
      hr11	94209521	94209522	AC	TA


duplication
^^^^^^^^^^^^^^^

Duplication can be thought of as special insertion where the inserted sequence is identical to the sequence flanking the breakpoint.
Similar to insertion, the annotation of duplication may possess alternative alignment.

Example: to annotate a duplication coding region,

.. code:: bash

   $ transvar canno --ccds -i 'CHD7:c.1669_1674dup'

::

   CHD7:c.1669_1674dup	CCDS47865 (protein_coding)	CHD7	+
      chr8:g.61693564_61693569dupCCCGTC/c.1669_1674dup/p.P558_S559dupPS	inside_[cds_in_exon_2]
      CSQN=InFrameInsertion;left_align_gDNA=g.61693561_61693562insTCCCCG;unalign_gD
      NA=g.61693562_61693567dupTCCCCG;left_align_cDNA=c.1668_1669insTCCCCG;unalign_
      cDNA=c.1669_1674dupTCCCCG;left_align_protein=p.H556_S557insSP;unalign_protein
      =p.S557_P558dupSP;phase=0;source=CCDS

Example: a duplication on the nucleotide level may lead to frame-shift or block substitution on the amino acid level,

.. code:: bash

   $ transvar canno --ccds -i 'CHD7:c.1668_1669dup'

::

   CHD7:c.1668_1669dup	CCDS47865 (protein_coding)	CHD7	+
      chr8:g.61693561_61693562dupTT/c.1668_1669dup/p.S557Ffs*8	inside_[cds_in_exon_2]
      CSQN=Frameshift;left_align_gDNA=g.61693560_61693561insTT;unalign_gDNA=g.61693
      561_61693562dupTT;left_align_cDNA=c.1667_1668insTT;unalign_cDNA=c.1668_1669du
      pTT;source=CCDS

Example: to annotate a duplication in intronic region,

.. code:: bash

   $ transvar canno --ccds -i 'CHD7:c.1666-5_1666-3dup'


::

   CHD7:c.1666-5_1666-3dup	CCDS47865 (protein_coding)	CHD7	+
      chr8:g.61693554_61693556dupCTC/c.1666-5_1666-3dup/.	inside_[intron_between_exon_1_and_2]
      CSQN=IntronicInsertion;left_align_gDNA=g.61693553_61693554insCTC;unalign_gDNA
      =g.61693554_61693556dupCTC;left_align_cDNA=c.1666-6_1666-5insCTC;unalign_cDNA
      =c.1666-5_1666-3dupCTC;source=CCDS


