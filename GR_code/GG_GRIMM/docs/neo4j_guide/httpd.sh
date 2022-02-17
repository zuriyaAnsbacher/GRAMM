#!/bin/sh

if [ -z "${GUIDES}" ]; then
  echo "Error:  Location of neo4j-guides is undefined."
  echo "Usage example:"
  echo "  GUIDES=\${HOME}/git/neo4j-guides ./httpd.sh"
  exit 1
fi

if [ ! -f ${GUIDES}/http-server.py ]; then
  echo "Error:  ${GUIDES}/http-server.py not found."
  exit 1
fi

if [ ! -d output/html ]; then
  echo "Error:  output/html subdirectory not found."
  exit 1
fi

cd output/html
exec ${GUIDES}/http-server.py
