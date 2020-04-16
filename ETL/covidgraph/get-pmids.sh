# Plain numbers, the format needed by Ondex
sed s/'"'//g "$DUMPDIR"/{abs,body}-genes.txt \
  | awk -F ',' '/pmid\.id/ { next } { print $1 }' \
  | sort | uniq | tee "$DUMPDIR"/pmids.txt
