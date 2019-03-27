***************************
Using non-canonical IDs
***************************

TransVar provides the use of non-canonical IDs by the means of ID mapping. This is achieved by providing an ID mapping file.

Create ID Mapping File
#######################

One can create an ID Mapping file by indexing a tab-delimited file with "synonym"(noncanonical ID) in the first column and canonical ID in the second column. The content of such tab-delimited file looks like

.. code:: text

	 MLL2 	KMT2D

   
And to create a ID mapping index

.. code:: bash

   transvar index --idmap [file_name] -o test.idmap_idx

Now you can use `--idmap` option in annotation to get annotation of non-canonical ID mapping

.. code:: bash

   transvar panno -i 'MLL2:p.Asp5492Asn' --ensembl --idmap test.idmap_idx

::

   MLL2:p.Asp5492Asn       ENST00000301067 (protein_coding)        KMT2D   -
   chr12:g.49415873C>T/c.16474G>A/p.D5492N inside_[cds_in_exon_53]
   CSQN=Missense;reference_codon=GAC;candidate_codons=AAC,AAT;candidate_mnv_varian
   ts=chr12:g.49415871_49415873delGTCinsATT;aliases=ENSP00000301067;source=Ensembl

You can see now TransVar can identify MLL2 which is a noncanonical ID in addition to the standard KMT2D.

Inherently, if you name the generated ID mapping to `[path_to_transvardb].XXX.idmap_idx`. You can use the shortcut of `--idmap XXX` as long as the annotation transcript database is provided. For example, `--idmap HGNC` when used with `--ensembl path_to_ensembl.transvardb` will also look for a ID mapping file of name `path_to_ensembl.transvardb.HGNC.idmap_idx`.
