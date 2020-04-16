cd $(dirname $0)
mkdir -p ../data-dumps/covidgraph
cat ./get-pmids.cypher \
  | $NEO_HOME/cypher-shell -a bolt://covid.petesis.com:7687 -u public -p corona --format=plain \
  | tee ../data-dumps/covidgraph/pmids.txt

