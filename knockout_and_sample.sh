#!/bin/bash

n_samp=1000
n_update=10
obj_flux=3.074

for S in $(seq 0 5)
do
    echo ${S}
    #echo $n_samp
    #echo $n_update
    Rscript --vanilla knockout_model_setup.R GitHub/tnseq_modeling/data/grb_updated_mutans_model.lp bio00001 ${S} ko_model.lp $obj_flux
    python sample_model.py ko_model.lp Documents/jensn\ lab/tnseq\ fitting/sampling_data/ko_${S}_samples.csv $n_samp $n_update
done

