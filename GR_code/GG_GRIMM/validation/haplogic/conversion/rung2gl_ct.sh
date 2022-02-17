#!/bin/bash

# download rel_dna_ser.txt and rel_ser_ser.txt from WMDA website

# generate rel_ser_dna.txt to expand serology to DNA
perl ds.pl

# convert CT validation dataset to GL String format
python3 g2gl_ct.py ./data/recip_dpre.in.txt ./data/don.gl.txt D
python3 g2gl_ct.py ./data/recip_dpre.in.txt ./data/pat.gl.txt R

gzip ./data/don.gl.txt
gzip ./data/pat.gl.txt