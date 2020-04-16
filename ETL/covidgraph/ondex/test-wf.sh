# This is how you should invoke the main Ondex workflow, once you have data downloaded from Neo4j
#
[[ -z "$ONDEX_HOME" ]] && {
	echo 'Cannot find ONDEX_HOME, download Ondex and initialise me with some *-init.sh';
	exit 1;
}

# the test main WF loads both abstract and body annotations and integrates them with publications imported
# via the PMIDs. Fragments of such WF should be copy-pasted into the WF that creates the whole OXL
# 
./runme.sh $SCRIPTDIR/ondex/test-wf.xml DUMPDIR="$DUMPDIR" MYDIR="$SCRIPTDIR/ondex"
