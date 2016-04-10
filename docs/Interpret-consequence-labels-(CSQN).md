For each genetic variant, TransVar assigns a consequence label with `CSQN` tag. The consequence label sometimes explains the behaviour of the output, e.g., the missing of protein level representation due to the loss of splice site.

The consequence label is in the following alphabet:

| label | interpretation |
|-------|----------------|
| Synonymous | Variation in protein-coding sequence results in the same protein sequence |
| Missense | Single or multiple amino acid substitution to protein-coding gene (1 to 1) |
| MultiAAMissense | In-frame multiple amino acid replacement (m to n, either m>1 or n>1) |
| Nonsense | Introduction of stop codon by single, multiple amino acid substitution or in-frame insertions/deletions |
|--- Coding Start/Stop ---|-----|
| CdsStartLoss | Loss of coding start |
| CdsStopLoss | Loss of coding stop |
|--- Coding Insertion/Deletion ---|-----|
| Frameshift | Frameshift mutation to a protein-coding gene |
| InFrameDeletion | In-frame deletion to a protein-coding gene |
| InFrameInsertion | In-frame insertion to a protein-coding gene |
|--- Intronic ---|-----|
| IntronicSNV | Intronic single nucleotide variation |
| IntronicDeletion | Intronic deletion |
| IntronicInsertion | Intronic insertion |
| IntronicBlockSubstitution | Intronic block substitution |
|--- Intergenic ---|-----|
| IntergenicSNV | Intergenic single nucleotide variation |
| IntergenicDeletion | Intergenic deletion |
| IntergenicInsertion | Intergenic insertion |
| IntergenicBlockSubstitution | Intergenic block substitution |
|--- Splice site ---|-----|
| SpliceDonorDeletion | Deletion occurs to splice donor |
| SpliceAcceptorDeletion | Deletion occurs to splice acceptor |
| SpliceDonorSNV | Single-nucleotide variation at splice donor |
| SpliceAcceptorSNV | Single-nucleotide variation at splice acceptor |
| SpliceDonorBlockSubstitution | Block substitution occurs at splice donor |
| SpliceAcceptorBlockSubstitution | Block substitution occurs at splice acceptor |
| SpliceDonorInsertion | Insertion at splice donor |
| SpliceAcceptorInsertion | Insertion at splice acceptor |
|--- Others ---|-----|
| Unclassified | Unclassified due to mostly mal-formated input |
