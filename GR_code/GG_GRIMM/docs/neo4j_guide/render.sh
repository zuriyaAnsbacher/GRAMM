#!/bin/sh

if [ -z "${GUIDES}" ]; then
  echo "Error:  Location of neo4j-guides is undefined."
  echo "Usage example:"
  echo "  GUIDES=\${HOME}/git/neo4j-guides ./render.sh"
  exit 1
fi

if [ ! -f ${GUIDES}/run.sh ]; then
  echo "Error:  ${GUIDES}/run.sh not found."
  exit 1
fi

rm -rf output
mkdir -p output/html
for F in index imputation-plan-a imputation-plan-b imputation-plan-b-redux ; do
  if [ ! -f guide/$F.adoc ]; then
    echo "Error:  guide/$F.adoc not found."
    exit 1
  fi
  ${GUIDES}/run.sh guide/$F.adoc output/html/$F.html +1 "$@"
done

cp -Rp guide/images output/html
