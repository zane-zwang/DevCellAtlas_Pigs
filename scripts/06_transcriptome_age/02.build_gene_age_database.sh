#!/usr/bin/env bash
# Purpose: build gene age database.


echo "Start......"
date

gunzip ./nr_db/nr.gz
#gunzip ./accession2taxid/prot.accession2taxid.gz

diamond makedb \
 --in ./nr_db/nr \
 --db ./nr_db/nr \
 --taxonmap ./accession2taxid/prot.accession2taxid \
 --taxonnodes taxdump/nodes.dmp \
 --taxonnames taxdump/names.dmp

date
echo "Done......"
