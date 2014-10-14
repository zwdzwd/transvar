#!/bin/bash
function er() { echo "$@"; "$@"; }
function er2() { echo "======="; echo "$@"; echo "======="; echo -e "\n";}

er transvar revanno --ccds -i 'Q5VUM1:47' --uniprot
er2 "Q5VUM1:47   6   71289191-71289193   CCDS4972.1  C6ORF57 (+, coding)    6:g.71289191_71289193/c.139_141/p.47    PRefSeq=S;NRefSeq=TCC;RefSeq=TCC"


er transvar revanno -i 'P28222.p.146_148refDRY' --uniprot --ccds
er2 "P28222.p.146_148refDRY  6       78172677-78172685       CCDS4986.1      HTR1B (-, coding)       6:g.78172677_78172685/c.436_444/p.146_148       PRefSeq=DRY;NRefSeq=GACCGCTAC;RefSeq=GTAGCGGTC"

er transvar revanno -i 'HTR1B.p.365_369refNPxxY' --ccds
er2 "HTR1B.p.365_369refNPxxY 6       78172014-78172028       CCDS4986.1      HTR1B (-, coding)       6:g.78172014_78172028/c.1093_1107/p.365_369     PRefSeq=NPIIY;NRefSeq=AAC..TAT;RefSeq=ATA..GTT"


er transvar revanno -i PIK3CA:E545K --ensembl
er2 'PIK3CA:E545K    3       178936091-178936092-178936093   ENST00000263967 PIK3CA (+, coding, exon 10)     3:g.178936091G>A/c.1633G>A/p.E545K      NCodonSeq=GAG;NCddSeqs=AAG,AAA;CddMNVMuts=3:g.178936091_178936093GAG>AAA'

er transvar revanno -i 'A1CF:p.A309A' --ccds --dbsnp
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
