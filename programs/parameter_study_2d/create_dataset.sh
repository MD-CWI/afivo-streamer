#!/bin/bash

rm -f dataset_128.h5
python create_dataset.py 'npz_128/*_e.npz' dataset_128.h5
python create_dataset.py 'npz_128/*_phi.npz' dataset_128.h5
python create_dataset.py 'npz_128/*_electric_fld.npz' dataset_128.h5
python create_dataset.py 'npz_128/*_rhs.npz' dataset_128.h5
python create_dataset.py 'npz_128/*_N2_B3.npz' dataset_128.h5
python create_dataset.py 'npz_128/*_N2_C3.npz' dataset_128.h5
