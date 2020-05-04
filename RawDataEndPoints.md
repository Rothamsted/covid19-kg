# CoViD-19 Kentminer, Raw Data Endpoints

The dataset backing our [Knetminer instance about CoViD 19][10] pandemics are also available as 
[SPARQL][20] and [Neo4j][30]. Dumps are also available for these databases, please ask
us if you're interested.

[10]: https://knetminer.rothamsted.ac.uk/COVID-19
[20]: http://knetminer-data.cyverseuk.org/lodestar/sparql
[30]: TODO

We are grateful to [Cyverse UK](https://cyverseuk.org/) for the server resources provided!

## Schema

The data are modelled after our own [BioKNO ontology][40], with some mapping to standard ontologies (we 
plan to do more mappings in future). We use that ontology as a common base for both RDF data and Neo4j 
data (which, by the way, are obtained from RDF using our [rdf2neo][50]).

The example queries below reflect that data model.

[40]: https://github.com/Rothamsted/bioknet-onto
[50]: https://github.com/Rothamsted/rdf2neo


## Example Queries

The following are examples of Cypher queries that can extract data from our [Neo4j database][30]. 
Semantically equivalent queries are available straight on our [SPARQL endpoint][20].


### Drugs occurring in known publications

```haskell
MATCH (drug:Drug) - [:occ_in] -> (pub:Publication) - [:identifier] - (acc)
WHERE acc.identifier IN [ 
  '32226295',
  '111B9A6E91C938696FCDB4CB128B8AE739DBE11C',
  '7852AAFDFB9E59E6AF78A47AF796325434F8922A',
  '6CF87A546884756094DA0D300E85A061C2CC43EA',
  'F863247DC84916A96448088399BADEFD09D54FB8'
]
RETURN drug.prefName AS drugName, acc.identifier as pubAcc, pub.AbstractHeader AS pubTitle, toInteger(pub.YEAR) as pubYear
ORDER BY pubYear DESC, pubAcc, drugName
LIMIT 25
```

### Known drugs, their targets and encoding genes

```haskell
MATCH (prot:Protein) - [:has_target] - (drug:Drug)
OPTIONAL MATCH (gene:Gene) - [:enc] -> (prot)
WHERE toLower(drug.prefName) IN [
  'opril',
  'minoxidil',
  'benzoyl peroxide',
  'isotretinoin',
  'trifluoperazine'
]
RETURN drug.prefName, prot.prefName, gene.identifier
ORDER BY drug.prefName, prot.prefName, gene.identifier
LIMIT 25
```

From this, we learn that there aren't genes encoding drug targets. But there are references:

```haskell
MATCH (gene:Gene) - [:enc] -> (xref:Protein) - [:xref] - (prot:Protein) <- [:has_target] - (drug:Drug),
(gene) - [:identifier] - (geneAcc:Accession)
WHERE toLower ( drug.prefName ) IN [
  'opril',
  'minoxidil',
  'benzoyl peroxide',
  'isotretinoin',
  'trifluoperazine']
RETURN DISTINCT drug.prefName AS drugName, prot.prefName AS protName, geneAcc.identifier AS geneAcc, xref.prefName AS xrefName
ORDER BY drugName, protName, geneAcc
LIMIT 25
```

### Most common pathways related to proteins that are targeted by literature-cited drugs 

This computes a table of the GO Biological processes that have the highest number of associated 
publications, where the association is made considering the literature-cited drugs that target the proteins involved in a process. 

```haskell
MATCH (bp:BioProc) <- [:participates_in] - (prot:Protein) <- [:has_target] - (drug:Drug) - [:occ_in] -> (pub:Publication)
RETURN bp.prefName, COUNT ( DISTINCT pub ) AS pubNo
ORDER BY pubNo DESC
LIMIT 10
```


### Distribution of no of disease citations, using indirect gene/protein mentions

This is another example of statistical figures computing. The first part of the query (before `UNWIND`)
gets the number of publications per associated disease (considering related genes and proteins mentioned 
in the publications). Then the second part computes cumulative frequencies about how many diseases have
fewer than a given publication count.

```haskell
MATCH (bioel) - [:inv_in] -> (dis:Disease),
(bioel) - [:occ_in|pub_in] -> (pub:Publication)
WHERE (bioel:Gene OR bioel:Protein)
WITH COUNT(DISTINCT pub) AS pubNo, dis
UNWIND [ 5, 10, 30, 50, 100, 1E9] AS maxPub
WITH pubNo, dis, maxPub
WHERE pubNo < maxPub
RETURN maxPub, COUNT( DISTINCT dis ) AS disNo
ORDER BY maxPub
```
