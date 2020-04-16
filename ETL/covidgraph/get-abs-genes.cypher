MATCH (paper:Paper) - [:PAPER_HAS_BIB_ENTRIES] -> (bib:Bib_entries) 
  - [:BIB_ENTRIES_HAS_BIBREF] -> (bibref:Bibref) 
  - [:BIBREF_HAS_OTHER_IDS] -> (id:Other_ids) 
  - [:OTHER_IDS_HAS_PMID] -> (pmidRef:Pmid) 
  - [:PMID_HAS_PMID] -> (pmid),
(paper:Paper) 
  - [absl:PAPER_HAS_ABSTRACT] -> (abs:Abstract) 
  - [:ABSTRACT_HAS_ABSTRACT] -> (abs1:Abstract) 
  - [frgl:HAS_FRAGMENT] -> (frag) 
  - [:MENTIONS] -> (gsymb:GeneSymbol) <- [:MAPS*1..2] - (gene:Gene)
WHERE gene.sid STARTS WITH "ENSG"
RETURN pmid.id, gene.sid;

