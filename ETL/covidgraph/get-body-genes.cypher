MATCH (paper:Paper) - [:PAPER_HAS_BIB_ENTRIES] -> (bib:Bib_entries)
  - [:BIB_ENTRIES_HAS_BIBREF] -> (bibref:Bibref)
  - [:BIBREF_HAS_OTHER_IDS] -> (id:Other_ids)
  - [:OTHER_IDS_HAS_PMID] -> (pmidRef:Pmid)
  - [:PMID_HAS_PMID] -> (pmid),
(paper:Paper)
  - [bodyl:PAPER_HAS_BODY_TEXT] -> (body:Body_text)
  - [:BODY_TEXT_HAS_BODY_TEXT] -> (body1:Body_text)
  - [frgl:HAS_FRAGMENT] -> (frag)
  - [:MENTIONS] -> (gsymb:GeneSymbol) <- [:MAPS*1..2] - (gene:Gene)
WHERE gene.sid STARTS WITH "ENSG"
RETURN pmid.id, gene.sid;
