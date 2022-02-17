
# Haplogic validation pipeline

## Imputation input files

Imputation input files are found at `validation/haplogic/data`

Valid files are donors/patients who have no serological data or mismatch.

| file  | Description |
| :-------------- | :------|
| `6loci_valid_don.gl.csv.gz` | 6 loci compressed csv file with donors typing and multi-race info (anonymized ID, GL-formatted typing, race1, race2)| 
| `5loci_valid_don.gl.csv.gz` | 5 loci compressed csv file with donors typing and multi-race info (anonymized ID, GL-formatted typing, race1, race2)| 
| `6loci_valid_pat.gl.csv.gz` | 6 loci compressed csv file with patients typing and single-race info (anonymized ID, GL-formatted typing, race, race) | 
| `5loci_valid_pat.gl.csv.gz` | 5 loci compressed csv file with patients typing and single-race info (anonymized ID, GL-formatted typing, race, race) | 
| `sample/don_100.txt.gz`   | 100 sample donors |
| `sample/pat_100.txt.gz`   | 100 sample patients |

## Configuration files

Configuration files are found in `../../conf/haplogic/`

| file  | Description |
| :-------------- | :------|
| `6loci-pat-configuration.json` | 6 loci patients JSON config. file  | 
| `6loci-don-configuration.json` | 6 loci donors  JSON config. file | 
| `5loci-pat-configuration.json` | 5 loci patients JSON config. file  | 
| `5loci-don-configuration.json` | 5 loci donors  JSON config. file | 


## Make commands
The following table shows different commands to use run the whole haplogic validation pipeline or an individual step.

| Command         | Description |
| :-------------- | :---------- |
| `make all`  | run the whole dataset validation pipeline |
| `make hpf`  | generate hpf.csv file and imputation graph csv files |
| `make parallel_impute` | run parallel imputation  |
| `make matching_graph`  | generate and load Neo4j matching graph |
| `make match`          | run pairwise matching queries for 9/10 aggregate, 10/10 & 9/10 and SLUG 2/2  |
| `make print`     | prints `$JAVA_HOME` & `$NEO4J_HOME` & Perl, Python version|

Running previous commands will assume using the default values for input arguments. To run a custom settings,`Makefile` has several input arguments. The following table summeries these arguments.

| input argument  | Description | Default|
| :-------------- | :---------- | :------|
| `GRIMM_DIR`  | sets GRIMM directory | `../../` |
| `PAT_CONFIG_FILE`|JSON config. file for patients|`$(GRIMM_DIR)/conf/6loci-pat-configuration.json`|
| `DON_CONFIG_FILE`|JSON config. file for donors|`$(GRIMM_DIR)/conf/6loci-don-configuration.json`|
| `PAT_IMPUTATION_OUT`| patients imputation output file | `$(GRIMM_DIR)/validation/output/pat.umug.freqs`|
| `DON_IMPUTATION_OUT`| donors imputation output file | `$(GRIMM_DIR)/validation/output/don.umug.freqs`|
| `PREDICTIONS_TO_CALCULATE`| predictions to use.  | `aggregate` |
| `FOLD_RESULTS_FILE` | fold results to merge 9/19 predictions | `matching/search/data/anon57K_fold_results.csv.gz` |

# Using Docker
## Install Docker

- For Mac/Windows install [Docker Desktop](https://www.docker.com/products/docker-desktop)
  
  By default, Docker Desktop only allocates 2GB of memory to Docker. Change Memory/CPU allocation in "Prefrences -> Resources -> Advanced" settings.

- For Linux, install using your package manger:

  On Redhat/Centos/EC2: Install docker and start it.
  ```
  sudo yum install -y docker
  sudo service docker restart
  ```

  On Ubuntu: Install docker and start it.
  ```
  sudo apt-get install -y docker
  sudo service docker restart
  ```

## Checkout code and run haplogic validation

  To run 6 Locus Imputation
  ```
  git clone https://github.com/nmdp-bioinformatics/graph-imputation-match
  cd graph-imputation-match
  ```
  To run 6 Locus Imputation
  ```
  docker run -d --name grimm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx haplogic 6loci
  ```
  To run 5 Locus Imputation
  ```
  docker run -d --name grimm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx haplogic 5loci
  ```

  You may need to use `sudo docker` depending upon your install.

  To look at logs:
  ```
  docker logs -f grimm
  ```

  To capture logs to a file:
  ```
  docker logs -f grimm > grimm_run.log
  ```

## Imputing your own data

 1. gzip your csv file and put it in a path somewhere
2. Create `conf/my-configuration.json` from `conf/haplogic/5loci-pat-configuration.json`
	Change `"imputation_in_file"` and `"populations"` settings for your data.
 3. Run Imputation in Docker
	 
    Setup and Run imputation
         
    ```
    docker run --name grimm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx make setup parallel_impute_pat PAT_CONFIG_FILE=/grimm/conf/my-configuration.json
    ```

    If you want to rerun, you can avoid setup with 
         
    ```
    docker run --name grimm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx make install_imputegl parallel_impute_pat PAT_CONFIG_FILE=/grimm/conf/my-configuration.json
    ```


## Manually run the pipeline in Docker

The image for the development is available from dockerhub: https://hub.docker.com/repository/docker/nmdpbioinformatics/grimm-networkx

It includes Python, Perl, Java and Neo4j Database.

```
Project directory is:  /grimm
=========================
JAVA_HOME is /usr/lib/jvm/java-11-openjdk-amd64
=========================
NEO4J_HOME is /usr/share/neo4j
=========================
checking python version:
Python 3.8.2
=========================
checking perl version:
 This is perl 5, version 28, subversion 1 (v5.28.1) built for x86_64-linux-gnu-thread-multi (with 61 registered patches, see perl -V for more detail)  Copyright 1987-2018, Larry Wall  Perl may be copied only under the terms of either the Artistic License or the GNU General Public License, which may be found in the Perl 5 source kit.  Complete documentation for Perl, including FAQ lists, should be found on this system using man perl or perldoc perl.  If you have access to the Internet, point your browser at http://www.perl.org/, the Perl Home Page.
=========================
```

- Enter the container to access the development environment
  ```
  cd graph-imputation-match
  docker run -it -v $PWD:/grimm nmdpbioinformatics/grimm-networkx bash
  ```


## Steps to run the pipeline

1. Validate your environment is setup correctly by running `make print`.
```
(venv) [ec2-user@ip-10-223-0-94 haplogic]$ make print
Project directory is:  /home/ec2-user/graph-imputation-match
=========================
JAVA_HOME is /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.222.b10-0.amzn2.0.1.x86_64
=========================
NEO4J_HOME is /home/ec2-user/neo4j-community-3.5.11
=========================
checking python version:
Python 3.7.4
=========================
checking perl version:
This is perl 5, version 16, subversion 3 (v5.16.3) built for x86_64-linux-thread-multi (with 39 registered patches, see perl -V for more detail)  Copyright 1987-2012, Larry Wall  Perl may be copied only under the terms of either the Artistic License or the GNU General Public License, which may be found in the Perl 5 source kit.  Complete documentation for Perl, including FAQ lists, should be found on this system using man perl or perldoc perl.  If you have access to the Internet, point your browser at http://www.perl.org/, the Perl Home Page.
=========================
```


2. Run `make`. This will:
 - Use 6loci donor/patient files
 - produce `hpf.csv` files from all frequency files. 
 - install all required python libraries.
 - run imputation for patients and donors.
 - match using aggregate queries to produce a validation results file in `matching/search/output/validation_results.tsv.gz`
 
You can specify the number of parallel processes to use. eg. To run 2 parallel processes, `make NUM_OF_PROCESSES=2` 

## Run with 5 loci files

The following will run the full pipleine with the 100K sample in haplogic/data/{pat,don}_100.csv
```
make PAT_CONFIG_FILE=$PWD/../../conf/haplogic/5loci-pat-configuration.json\
     DON_CONFIG_FILE=$PWD/../../conf/haplogic/5loci-don-configuration.json
```

## Run 100K Sample files

The following will run the full pipleine with the 100K sample in haplogic/data/{pat,don}_100.csv
```
make PAT_CONFIG_FILE=$PWD/../../conf/haplogic/pat-100-configuration.json\
  DON_CONFIG_FILE=$PWD/../../conf/haplogic/don-100-configuration.json\
  FOLD_RESULTS_FILE=$PWD/../../matching/search/data/anon100_fold.csv.gz
```
