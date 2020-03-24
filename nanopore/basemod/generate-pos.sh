#!/bin/sh
pos=$1
ORF="E M N ORF1ab ORF3a ORF6 ORF7a ORF7b ORF8 S nsp2trunc nsp3Ctrunc"
CONTROL_FAST5_DIR='/blaze/hyeshik/p/npworks/20200222-Vero-SCV2/subguppy-leader/IVT1-single/'
FULL_LENGTH_GROUP='SARSCoV2FullGenomeCorrected_000'
OUTPUTDIR='/qbio/nest/hyeshik/working/20200222-Vero-SCV2/basemod-figures'

for orf in $ORF; do
echo $orf
tombo plot genome_locations --fast5-basedirs flg-tombo/$orf \
    --control-fast5-basedirs $CONTROL_FAST5_DIR \
    --plot-standard-model \
    --corrected-group $FULL_LENGTH_GROUP \
    --genome-locations chrSCV:$pos \
    --pdf-filename $OUTPUTDIR/genomelocation-$pos-$orf.pdf
done
