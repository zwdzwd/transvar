
mkdir -p testout

transvar panno -l cosmic/cosmic_panno_frameshift -g 1 -m 4 --ccds | tee testout/panno_cosmic_frameshift

transvar panno -l data/cosmic/cosmic_panno_deletion -g 1 -m 4 --ccds | tee testout/panno_cosmic_deletion

transvar panno -l data/cosmic/cosmic_panno_insertion -g 1 -m 4 --ccds | tee testout/panno_cosmic_insertion

transvar panno -l data/cosmic/cosmic_panno_saav -g 1 -m 4 --ccds | tee testout/panno_cosmic_saav

transvar canno -l data/cosmic/cosmic_canno_snv -g 1 -m 3 --ccds | tee testout/canno_cosmic_snv

transvar canno -l data/cosmic/cosmic_canno_deletion -g 1 -m 3 --ccds | tee testout/canno_cosmic_deletion

transvar canno -l data/cosmic/cosmic_canno_insertion -g 1 -m 3 --ccds | tee testout/canno_cosmic_insertion





