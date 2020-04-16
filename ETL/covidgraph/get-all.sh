cd $(dirname $0)

[[ -z "$NEO_HOME" ]] && { echo 'Cannot find NEO_HOME, download Neo4j and initialise me with some *-init.sh'; exit 1; }

echo "Getting abstract-cited genes"
./get-abs-genes.sh >/dev/null

echo "Getting body-cited genes"
./get-body-genes.sh >/dev/null


echo "Building PMID list from references in previous results"
./get-pmids.sh >/dev/null
