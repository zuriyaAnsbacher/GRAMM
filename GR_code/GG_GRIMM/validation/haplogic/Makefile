# Project Directory
GRIMM_DIR=$(realpath ../../)

# Haplogic Patient and Donor Configuration
PAT_CONFIG_FILE=$(GRIMM_DIR)/conf/haplogic/6loci-pat-configuration.json
DON_CONFIG_FILE=$(GRIMM_DIR)/conf/haplogic/6loci-don-configuration.json

# Needed for perl graph match generator. Needs to match config.json
DON_IMPUTATION_OUT=$(GRIMM_DIR)/validation/output/don.umug.freqs
PAT_IMPUTATION_OUT=$(GRIMM_DIR)/validation/output/pat.umug.freqs

# Matching Configuration
PREDICTIONS_TO_CALCULATE=aggregate
FOLD_RESULTS_FILE=$(GRIMM_DIR)/matching/search/data/valid_fold_results.csv.gz
TRUE_RESULTS_PATH=$(GRIMM_DIR)/imputation/graph_generation/data/wmda/set3.consensus.txt

# Number of PROCESSES
# 0 -> Use num of cpus - 1
NUM_OF_PROCESSES=0

all: setup neo4j_grimm

setup: clean_all \
	hpf \
	install_imputegl \
	install_py2neo

networkx_grimm: parallel_impute \
	networkx_matching

neo4j_grimm: parallel_impute \
	neo4j_matching_graph \
	neo4j_match

clean_all:
	@echo "Removing all output directories"
	@echo "Removing $(GRIMM_DIR)/imputation/graph_generation/output"
	@rm -rf $(GRIMM_DIR)/imputation/graph_generation/output
	@echo "Removing $(GRIMM_DIR)/validation/output"
	@rm -rf $(GRIMM_DIR)/validation/output
	@echo "Removing $(GRIMM_DIR)/matching/graph_generation/perl/output"
	@rm -rf $(GRIMM_DIR)/matching/graph_generation/perl/output
	@echo "Removing $(GRIMM_DIR)/matching/search/output"
	@rm -rf $(GRIMM_DIR)/matching/search/output

hpf:
	@echo "generating hpf and imputation graph files"
	@echo $(shell date)
	@cd $(GRIMM_DIR)/imputation/graph_generation && $(MAKE) em2014 CONFIG_FILE=$(PAT_CONFIG_FILE)

install_imputegl:
	@cd $(GRIMM_DIR)/imputation && pip uninstall -y imputegl && python setup.py install && pip install imputegl

install_py2neo:
	@pip install py2neo

parallel_impute: parallel_impute_pat parallel_impute_don

parallel_impute_pat:
	@echo $(shell date)
	@echo "Patients imputation"
	@cd $(GRIMM_DIR)/validation && python parallel-imputation.py -c $(PAT_CONFIG_FILE) -p  $(NUM_OF_PROCESSES)
	@cd $(GRIMM_DIR)/validation && python reduce_loci.py -c $(PAT_CONFIG_FILE)

parallel_impute_don:
	@echo $(shell date)
	@echo "Donors imputation"
	@cd $(GRIMM_DIR)/validation && python parallel-imputation.py -c $(DON_CONFIG_FILE) -p $(NUM_OF_PROCESSES)
	@cd $(GRIMM_DIR)/validation && python reduce_loci.py -c $(DON_CONFIG_FILE)

networkx_matching:
	@echo $(shell date)
	@cd $(GRIMM_DIR)/matching/search && \
	mkdir -p output && \
  	mkdir -p ../graph_generation/networkx/perl/output/graph && \
	python3 graph_match_wmda.py $(DON_IMPUTATION_OUT) $(PAT_IMPUTATION_OUT) $(TRUE_RESULTS_PATH)

neo4j_matching_graph:
	@echo $(shell date)
	@cd $(GRIMM_DIR)/matching/graph_generation/neo4j/perl && mkdir -p output/graph
	@echo "Generating CSV graph files from" $(DON_IMPUTATION_OUT) and $(PAT_IMPUTATION_OUT) "imputation files"
	@cd $(GRIMM_DIR)/matching/graph_generation/neo4j/perl && ./generate_matchgraph.pl  $(DON_IMPUTATION_OUT) $(PAT_IMPUTATION_OUT)
	@echo "Loading generated graph into Neo4J"
	@cd $(GRIMM_DIR)/matching/graph_generation/neo4j/perl && ./bulk_load_neo4j.sh
	@echo "Waiting 10 seconds for Neo4J database to come up...................."
	@sleep 10

neo4j_match:	
	@echo $(shell date)
	@cd $(GRIMM_DIR)/matching/search && mkdir -p output && \
    	python parallel-pairwise_matching.py -p $(PREDICTIONS_TO_CALCULATE) -f $(FOLD_RESULTS_FILE)

print:
	@echo	"Project directory is:  $(GRIMM_DIR)"
	@echo   "========================="
	@echo   "JAVA_HOME is $(JAVA_HOME)"
	@echo   "========================="
	@echo   "NEO4J_HOME is $(NEO4J_HOME)"
	@echo   "========================="
	@echo   "checking python version:"
	@echo	"$(shell python --version)"
	@echo   "========================="
	@echo   "checking perl version: "
	@echo	"$(shell perl --version)"
	@echo   "========================="



