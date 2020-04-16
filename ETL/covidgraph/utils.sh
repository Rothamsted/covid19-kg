export CYCMD="$NEO_HOME/cypher-shell -a bolt://covid.petesis.com:7687 -u public -p corona --format=plain"
mkdir -p "$DUMPDIR"

function cypher
{
	query="$1"
	out="$2"

	cat "$query" | $CYCMD | tee "$DUMPDIR/$out"
}

export -f cypher
