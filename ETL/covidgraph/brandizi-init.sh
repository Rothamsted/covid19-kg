export NEO_HOME=/Applications/local/dev/semantic_web/neo4j-community-3.5.12/bin
export ONDEX_HOME=/Users/brandizi/Documents/Work/RRes/ondex_git/ondex-full/ondex-knet-builder/ondex-mini/target/ondex-mini

export SCRIPTDIR="$(cd $(dirname "${BASH_SOURCE[0]}") && pwd)"

mydir="$(pwd)"
cd "$SCRIPTDIR/../data-dumps/covidgraph"
export DUMPDIR="$(pwd)"
cd "$mydir"

. "$SCRIPTDIR/utils.sh"
