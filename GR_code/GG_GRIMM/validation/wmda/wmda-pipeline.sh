#!/usr/bin/env bash

CONF_DIR=../../conf
GRAPH_ENGINE="neo4j"
# GRAPH_ENGINE should be neo4j/networkx

# use WMDA haplotype frequencies
cd ../../imputation/graph_generation
make wmda CONFIG_FILE=$CONF_DIR/wmda/wmda-pat.json

# uninstall and reinstall imputegl
cd ../
pip uninstall -y imputegl

python setup.py install
pip install imputegl

# run imputation
cd ../validation
CONF_DIR=../conf
python runfile.py -c  $CONF_DIR/wmda/wmda-pat.json
python runfile.py -c  $CONF_DIR/wmda/wmda-don.json

if [ "$GRAPH_ENGINE" = "networkx" ]
then
  # generate networkx matching graph, perform matching, and compare results
  cd ../matching/search
  mkdir output
  mkdir -p ../graph_generation/networkx/perl/output/graph
  python graph_match_wmda.py ../../validation/output/don.umug.freqs ../../validation/output/pat.umug.freqs ../../imputation/graph_generation/data/wmda/set3.consensus.txt
else
  # generate neo4j matching graph
  cd ../matching/graph_generation/neo4j/perl
  mkdir -p output/graph
  ./generate_matchgraph.pl ../../../../validation/output/don.umug.freqs ../../../../validation/output/pat.umug.freqs
  ./bulk_load_neo4j.sh

  # wait for neo4j db to load
  sleep 10 #zzzz

  # run matching queries
  cd ../../../search
  mkdir -p output
  ./match_results.pl

  # compare with WMDA reference results
  ./compare_results.pl
fi
