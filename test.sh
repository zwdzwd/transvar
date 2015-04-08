#!/bin/bash
function er() { echo "$@"; "$@"; }
function er2() { echo "======="; echo "$@"; echo "======="; echo "\n";}

er transvar panno --ccds -i 'Q5VUM1:47' --uniprot
er2 "Q5VUM1:47       CCDS4972.1 (protein_coding)     C6ORF57 +       chr6:g.71289191_71289193/c.139_141/p.47S        cds_in_exon_2   protein_sequence=S;cDNA_sequence=TCC;gDNA_sequence=TCC"

er transvar panno -i 'P28222:p.146_148refDRY' --uniprot --ccds
er2 "
P28222:p.146_148refDRY  CCDS4986.1 (protein_coding)     HTR1B   -       chr6:g.78172677_78172685/c.436_444/p.D146_Y148  cds_in_exon_1   protein_sequence=DRY;cDNA_sequence=GACCGCTAC;gDNA_sequence=GTAGCGGTC"

er transvar revanno -i 'HTR1B.p.365_369refNPxxY' --ccds
er2 "HTR1B.p.365_369refNPxxY 6       78172014-78172028       CCDS4986.1      HTR1B (-, coding)       6:g.78172014_78172028/c.1093_1107/p.365_369     PRefSeq=NPIIY;NRefSeq=AAC..TAT;RefSeq=ATA..GTT"


er transvar revanno -i PIK3CA:E545K --ensembl
er2 'PIK3CA:E545K    3       178936091-178936092-178936093   ENST00000263967 PIK3CA (+, coding, exon 10)     3:g.178936091G>A/c.1633G>A/p.E545K      NCodonSeq=GAG;NCddSeqs=AAG,AAA;CddMNVMuts=3:g.178936091_178936093GAG>AAA'

er transvar revanno -i 'A1CF:p.A309A' --ccds
er2 'A1CF:p.A309A    10    52576004-52576005-52576006    CCDS7243.1    A1CF (-, coding, exon 7)    10:g.52576004T>G/c.927A>C/p.A309A    NCodonSeq=GCA;NCddSeqs=GCC,GCG,GCT;CddSNVMuts=10:g.52576004T>C,10:g.52576004T>A; DBSNP=rs201831949(10:52576004T>G)'

er transvar revanno --ccds -i 'PIK3CA:c.1633G>A'
er2 'PIK3CA:c.1633G>A     3     178936091-178936092-178936093   CCDS43171.1	PIK3CA (+ coding)  3:G178936091A/c.1633G>A/p.E545K NCodonSeq=GAG;NAltCodonSeq=AAG'

er transvar revanno --ccds -i 'ACIN1:c.1932_1933insATTCAC'
er2 'ACIN1:c.1932_1933insATTCAC    14     14:23548785-(ins)-23548786    CCDS55905.1    ACIN1 (-, coding)    14:23548785_23548786insGTGAAT/c.1932_1933insATTCAC/p.R644_S645insIH    NatInsSeq=ATTCAC;RefInsSeq=GTGAAT;Phase=0'

er transvar revanno --ccds -i 'ACIN1:c.1930_1931insATTCAC'
er2 'ACIN1:c.1930_1931insATTCAC   14     14:23548787-(ins)-23548788      CCDS9587.1    ACIN1 (-, coding)       14:23548787_23548788insGTGAAT/c.1930_1931insATTCAC/p.S643_R644insHS    NatInsSeq=C(ATTCAC)GT;RefInsSeq=GTGAAT;Phase=1'

er transvar revanno --ccds -i 'AAAS:c.1225_1226insG'
er2 'AAAS:c.1225_1226insG    12   12:53702089-53702090    CCDS8856.1    AAAS (-, coding)    12:53702089_53702090insC/c.1225_1226insG/p.E409Gfs*17   NatInsSeq=G;RefInsSeq=C'

er transvar revanno --ccds -i 'ADAM33:c.991-3_991-2insC'
er2 'ADAM33:c.991-3_991-2insC   20   20:3654141-3654142-3654143-(3654145)-(ins)-(3654146)    CCDS13058.1     ADAM33 (-, intronic)     20:3654145_3654146insG/c.991-3_991-2insC/.   RefInsSeq=G;NatInsSeq=C'


er transvar revanno --ccds -i 'A4GNT:c.694_696delTTG'
er2 'A4GNT:c.694_696delTTG   3   137843433-137843435     CCDS3097.1    A4GNT (- coding)  3:137843433_137843435del/c.694_696del/p.L232del RefDelSeq=CAA;NatDelSeq=TTG'

er transvar revanno --ccds -i 'ABHD15.c.431_433delGTG'
er2 'ABHD15.c.431_433delGTG  17   27893552-27893554    CCDS32602.1     ABHD15 (- coding)  17:27893552_27893554del/c.431_433del/p.C144_V145delinsF RefDelSeq=CAC;NatDelSeq=GTG'

er transvar revanno --ccds -i 'AADACL3.c.374delG'
er2 'AADACL3.c.374delG    1    12785494-12785494    CCDS41252.1   AADACL3 (+ coding)    1:12785494_12785494del/c.374del/p.C125Ffs*17    RefDelSeq=G;NatDelSeq=G'

er transvar revanno --ccds -i 'ABCB11:c.1198-8_1199delcactccagAA'
er2 'ABCB11:c.1198-8_1199delcactccagAA       2       2:169833196-169833205    CCDS46444.1     ABCB11 (- coding & intronic)    2:169833196_169833205del/c.1198-8_1199del/p.K400Tfs*4    RefDelSeq=TTCTGGAGTG;NatDelSeq=CACTCCAGAA'

er transvar revanno --ccds -i 'A1CF:c.508_509CC>TT'
er2 'A1CF:c.508_509CC>TT     10      52595929-52595930       CCDS7241.1      A1CF (-, coding)        10:52595929_52595930GG>AA/c.508_509CC>TT/p.P170L        .
A1CF:c.508_509CC>TT     10      52595929-52595930       CCDS7242.1      A1CF (-, coding)        10:52595929_52595930GG>AA/c.508_509CC>TT/p.P170L        .'

er transvar revanno --ccds -i 'CSRNP1.c.1212_1224>GGAGGAGGAA'
er2 'CSRNP1.c.1212_1224>GGAGGAGGAA   3    39185092-39185104   CCDS2682.1   CSRNP1 (-, coding)   3:39185092_39185104TTCCTCCTCCTCC>TTCCTCCTCC/c.1212_1224GGAGGAGGAGGAA>GGAGGAGGAA/p.E4'

er transvar revanno --ccds -i 'A1CF:c.1460+2_1460+3TG>CC'
er2 'A1CF:c.1460+2_1460+3TG>CC    10    52570797-52570798   CCDS7241.1    A1CF (-, intronic)    10:52570797_52570798CA>GG/c.1460+2_1460+3TG>CC/.        .'

er transvar revanno --ccds -i 'A1CF:c.1459_1460+3ATGTG>CC'
er2 'A1CF:c.1459_1460+3ATGTG>CC    10   52570797-52570801   CCDS7241.1   A1CF (-, coding;intronic)	10:52570797_52570801CACAT>GG/c.1459_1460+3ATGTG>CC/.    CrossSplitSite'

er transvar revanno --ccds -i 'CHD7:c.1669_1674dup'
er2 'CHD7:c.1669_1674dup    8    61693562-61693567 (dup) CCDS47865.1     CHD7 (+, Coding)    8:61693562-61693567dupTCCCCG/c.1669_1674dup/p.S557_P558dupSP   RefDupSeq=TCCCCG;NatDupSeq=TCCCCG'

er transvar revanno --ccds -i 'CHD7:c.1668_1669dup'
er2 'CHD7:c.1668_1669dup    8    61693561-61693562 (dup) CCDS47865.1     CHD7 (+, Coding)    8:61693561-61693562dupTT/c.1668_1669dup/p.S557Ffs*8   RefDupSeq=TT;NatDupSeq=TT'

er transvar revanno --ccds -i 'CHD7:c.1666-5_1666-3dup'
er2 'CHD7:c.1666-5_1666-3dup 8   61693554-61693556 (dup) CCDS47865.1   CHD7 (+, Intronic)    8:61693554-61693556dupCTC/c.1666-5_1666-3dup/.  RefDupSeq=CTC;NatDupSeq=CTC'

er transvar revanno --ccds -i 'AATK.p.P1331_A1332insTP'
er2 'AATK    c.3993_3994insACGCCC   p.P1331_A1332insTP   17:79093270-79093271    17    79093268-79093273 (insertion)   CCDS45807.1     AATK (-, coding)   17:79093268-79093273ins6/c.3991-3996ins6/p.P1331_A1332insTP     Uncertain'

er transvar revanno --ccds -i 'AADACL4.p.W263_I267delWRDAI'
er2 'AADACL4   c.788_802del15  p.W263_I267delWRDAI   1:12726310-12726324     1       12726309-12726323 (deletion)    CCDS30590.1     AADACL4 (+, coding)     1:12726309-12726323/c.787-801/p.W263_I267delWRDAI     Uncertain'

er transvar revanno --ccds -i 'ABCC3:p.Y556_V557delinsRRR'
er2 'ABCC3:p.Y556_V557delinsRRR   17   48745254-48745259 (block substitution)  CCDS32681.1    ABCC3 (+, coding)    17:48745254-48745259TACGTG>AGGAGGAGG/c.1666-1671TACGTG>AGGAGGAGG/p.Y556_V557delinsRRR    CddNatAlt=AGG/AGA/CGA/CGC/CGG/CGT+AGG/AGA/CGA/CGC/CGG/CGT+AGG/AGA/CGA/CGC/CGG/CGT;Uncertain'

er transvar revanno --ccds -i 'A1BG.p.G132fs*2'
er2 'A1BG.p.G132fs*2 19      58863866-58863867-58863868      CCDS12976.1     A1BG (-, coding)    19:58863860-58863868/c.394-402/p.G132fs*2       RoughEstimateFromFrameShift'

er transvar codonsearch --ccds -i CDKN2A.p.58
er2 'CDKN2A.p.58    CDKN2A.p.72   9    21971184-21971185-21971186   21971185-21971186-21971187
    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p.58    CDKN2A.p.73   9    21971184-21971185-21971186   21971182-21971183-21971184
    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]'

####### forward annotation #####
er transvar anno --ccds -i 'chr3:178936091.G>A'
er2 'chr3:178936091.G>A   3   178936091-178936092-178936093   CCDS43171.1     PIK3CA (+, coding)      3:G178936091A/c.1633>/p.E545K   .'

er transvar anno -i "9:135782704C>G" --ccds
er2 '9:135782704C>G  9    135782704    CCDS6956.1    TSC1 (-, coding)    9:g.135782704C>G/c.1317G>C/p.L439L    CodonPos=135782704-135782705-135782706;NCodonSeq=CTG
9:135782704C>G  9    135782704    CCDS55350.1   TSC1 (-, coding)    9:g.135782704C>G/c.1164G>C/p.L388L    CodonPos=135782704-135782705-135782706;NCodonSeq=CTG'

er transvar anno --ccds -i 'chr3:g.178936091_178936192'
er2 'chr3:g.178936091_178936192   3    178936091-178936192    CCDS43171.1    PIK3CA (+, coding,intronic)   3:g.178936091_178936192/c.1633_1664+70/p.E545_R555    BEGCodon=178936091-178936092-178936093;ENDCodon=178936121-178936122-178936984'

er transvar anno -i '9:g.133750356_137990357' --ccds
er2 '9:g.133750356_137990357 9     133750356-137990357     BEG=CCDS35165.1,END=CCDS6986.1    4,240,002 bp covering 53 genes  9:g.133750356_137990357/./.    BEGreg=ABL1 (+, coding);BEGid=9:g.133750356A>/c.1244A>/p.H415;	ENDreg=OLFM1 (+, intronic);ENDid=9:g.137990357C>/c.622+6C>/.'

er transvar anno -i '9:g.133750356_1337503570' --ccds
er2 '9:g.133750356_1337503570    9    133750356-141213431    BEG=CCDS35165.1,END=.    7,463,076 bp covering 137 genes 9:g.133750356_141213431/./.    BEGreg=ABL1 (+, coding);BEGid=9:g.133750356A>/c.1244A>/p.H415;	ENDreg=Noncoding (up: 484,026 bp to EHMT1, down: 0 bp to 3-telomere);ENDid=././.'
