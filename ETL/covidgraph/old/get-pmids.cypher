MATCH (paper:Paper) - [:PAPER_HAS_BIB_ENTRIES] -> (bib:Bib_entries) 
  - [:BIB_ENTRIES_HAS_BIBREF] -> (bibref:Bibref) 
  - [:BIBREF_HAS_OTHER_IDS] -> (id:Other_ids) 
  - [:OTHER_IDS_HAS_PMID] -> (pmidRef:Pmid) 
  - [:PMID_HAS_PMID] -> (pmid)
RETURN pmid.id;

