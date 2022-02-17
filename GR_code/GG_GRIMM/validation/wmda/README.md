# WMDA validation 

## Using Docker
### Install Docker

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

### Checkout code and run WMDA validation

  ```
  git clone https://github.com/nmdp-bioinformatics/graph-imputation-match
  cd graph-imputation-match
  ```

  Run WMDA validation
  ```
  docker run -d --name grimm -v $PWD:/grimm nmdpbioinformatics/grimm-networkx:0.0.1 wmda
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

## Without Docker

To run WMDA validation pipeline (make sure to set `JAVA_HOME`, `NEO4J_HOME` and required Perl libraries)

```
./wmda-pipeline.sh 
```
This will run imputation, matching queries and compare results with WMDA reference matching results. 
1. Matching results can be found at `/matching/search/output/mr.txt` 
2. comparison results can be found at `/matching/search/output/cmp.txt`

This script assumes that donor/patient input files are `data/don.gl.txt` & `data/pat.gl.txt` and there are two config. files `wmda-pat.json` & `wmda-don.json` 

## NOTE
Make sure donor/patient input files are comma-separated and have population information as follows 
`subject_id,typing_in_gl_format,CAU,CAU` 

If not, use `./code/input-converter.py` to convert wmda files (it accepts both .txt and .gz files).

```
cd code/
python  input-converter.py input_file
```

