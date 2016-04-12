********************
Setup and Customize
********************

Use environment variable to direct download and configuration
##############################################################

TRANSVAR_CFG
^^^^^^^^^^^^^^

store the path to transvar.cfg

.. code:: bash

   export TRANSVAR_CFG=path_to_transvar.cfg

If not specified, TransVar will use **[installdir]/lib/transvar/transvar.cfg** directory or your local **~/.transvar.cfg** if the installation directory is inaccessible.

TRANSVAR_DOWNLOAD_DIR
^^^^^^^^^^^^^^^^^^^^^^^^

store the path to the directory where auto-download of annotation and reference go

.. code:: bash

   export TRANSVAR_DOWNLOAD_DIR=path_to_transvar_download_directory

If not specified, TransVar will use **[installdir]/lib/transvar/transvar.download** directory or your local **~/.transvar.download** if the installation directory is inaccessible.

Install and specify reference genome assembly
###############################################

Download from TransVar database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For some genome assembly (currently hg18, hg19, hg38, mm9 and mm10) we provide download via

.. code:: bash

   transvar config --download_ref --refversion [reference name]

See ``transvar config -h`` for all choices of ``[reference name]``).

Manual download and index
^^^^^^^^^^^^^^^^^^^^^^^^^^

For other genome assemblies, one could manually download the genome and index it manually by, 

.. code:: bash

   transvar index --reference [fasta]

Under the hood, TransVar uses the "samtools faidx". So one could use any existing faidx indices without a glitch.
Once downloaded and indexed, the genome can be used through the "--reference" option followed by path to the genome.

To set the default location of genome file for a reference version, say, to ./hg19.fa,

.. code:: bash

   transvar config -k reference -v ./hg19.fa --refversion hg19

will create in transvar.cfg an entry

::
   
   [hg19]
   reference = hg19.fa

so that there is no need to specify the location of reference on subsequent usages.

Install and specify transcript annotations
############################################

Download from TransVar database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One could automatically download transcript annotations via E.g., 

.. code:: bash

   transvar config --download_anno --refversion hg19

which download annotation from TransVar database to **[installdir]/lib/transvar/transvar.download** directory or your local **~/.transvar.download** if the installation directory is inaccessible. See **transvar config -h** for all version names.
These will also create default mappings under the corresponding reference version section of **transvar.cfg** like

::
   
   [hg19]
   ucsc = /home/wzhou1/download/hg19.ucsc.txt.gz

Download from Ensembl ftp
^^^^^^^^^^^^^^^^^^^^^^^^^^

One also has the option of downloading from Ensembl collection.

.. code:: bash

   transvar config --download_ensembl --refversion mus_musculus

Without specifying the refversion, user will be prompted a collection of options to choose from.

Know Current configuration
###########################

One can read the transvar.cfg file for the information. Alternatively one may run

.. code:: bash

   transvar current

which returns information about the setup regarding to the current reference selection, including the location of the reference file and database file.

::
   
   Current reference version: mm10
   reference: /home/wzhou/genomes_link/mm10/mm10.fa
   Available databases:
   refseq: /home/wzhou/tools/transvar/transvar/transvar.download/mm10.refseq.gff.gz
   ccds: /home/wzhou/tools/transvar/transvar/transvar.download/mm10.ccds.txt
   ensembl: /home/wzhou/tools/transvar/transvar/transvar.download/mm10.ensembl.gtf.gz

specifying ``--refversion`` displays the information under that reference version (without changing the default reference version setup).

Use Additional Resources
##################################

TransVar uses optional additional resources for annotation.

dbSNP
^^^^^^^

For example, one could annotate SNP with dbSNP id by downloading the dbSNP files.
This can be done by

.. code:: bash

   transvar config --download_dbsnp

TransVar automatically download dbSNP file which correspoding to the current default reference version (as set in **transvar.cfg**). This also sets the entry in **transvar.cfg**.
With dbSNP file downloaded, TransVar automatically looks for dbSNP id when performing annotation.

.. code:: bash

   transvar panno -i 'A1CF:p.A309A' --ccds

::

   A1CF:p.A309A	CCDS7243 (protein_coding)	A1CF	-
      chr10:g.52576004T>G/c.927A>C/p.A309A	inside_[cds_in_exon_7]
      CSQN=Synonymous;reference_codon=GCA;candidate_codons=GCC,GCG,GCT;candidate_sn
      v_variants=chr10:g.52576004T>C,chr10:g.52576004T>A;dbsnp=rs201831949(chr10:52
      576004T>G);source=CCDS

Note that in order to use dbSNP, one must download the dbSNP database through

.. code:: bash

   transvar config --download_dbsnp

or by configure the ``dbsnp`` slot in the configure file via

.. code:: bash

   transvar config -k dbsnp -v [path to dbSNP VCF]

Manually set path for dbSNP file must have the file tabix indexed.

Control the length of reference sequence
##########################################

TransVar reduces the reference sequence in a deletion to its length when the deleted reference sequence is too long. For example

.. code:: bash

   $ transvar ganno -i 'chr14:g.101347000_101347023del' --ensembl

outputs

::

   chr14:g.101347000_101347023del	ENST00000534062 (protein_coding)	RTL1	-
      chr14:g.101347000_101347023del24/c.4074+29_4074+52del24/.	inside_[3-UTR;noncoding_exon_1]
      CSQN=3-UTRDeletion;left_align_gDNA=g.101347000_101347023del24;unaligned_gDNA=
      g.101347000_101347023del24;left_align_cDNA=c.4074+29_4074+52del24;unalign_cDN
      A=c.4074+29_4074+52del24;aliases=ENSP00000435342;source=Ensembl

where the deletion sequence is reduced to its length (`del24`). The `--seqmax` option changes the length threshold (default:10) when this behavior occur. When `--seqmax` is given a negative number, the threshold is lifted such that the reference sequence is always reported regardless of its length, i.e.,

.. code:: bash

   $ transvar ganno -i 'chr14:g.101347000_101347023del' --ensembl --seqmax -1

outputs the full reference sequence:

::

   chr14:g.101347000_101347023del	ENST00000534062 (protein_coding)	RTL1	-
      chr14:g.101347000_101347023delTTGGGGTGAGAAATAGAGGGGACT/c.4074+29_4074+52delAGTCCCCTCTATTTCTCACCCCAA/.	inside_[3-UTR;noncoding_exon_1]
      CSQN=3-UTRDeletion;left_align_gDNA=g.101347000_101347023delTTGGGGTGAGAAATAGAG
      GGGACT;unaligned_gDNA=g.101347000_101347023delTTGGGGTGAGAAATAGAGGGGACT;left_a
      lign_cDNA=c.4074+29_4074+52delAGTCCCCTCTATTTCTCACCCCAA;unalign_cDNA=c.4074+29
      _4074+52delAGTCCCCTCTATTTCTCACCCCAA;aliases=ENSP00000435342;source=Ensembl
