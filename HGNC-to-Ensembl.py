import argparse
import pandas as pd

def mapID(gene_file, map_file, output):
    map_df = pd.read_csv(map_file, sep="\t")
    map_df.columns = ['Gene Stable ID', 'Exon Stable ID', 'id']
    del map_df['Exon Stable ID']

    gene_df = pd.read_csv(gene_file, sep="\t")
    merged_pd = gene_df.merge(map_df, on=['id'], how='right')
    merged_pd = merged_pd.dropna().drop_duplicates()
    merged_pd = merged_pd[['Gene Stable ID', 'id', 'name', 'paper_id']]
    merged_pd.to_csv(output, sep="\t", index=None, header=True)

parser = argparse.ArgumentParser(
    description="This script will convert Kaggle data from JSON format to table using metadata only.\n")
parser.add_argument(
    '-gf', '--geneFile', help='Specifies the file for the Gene output form CORD19 data using the ETL I made', required=True)
parser.add_argument(
    '-o', '--output', help='Specifies the output name and directory of the resultant file', required=True)
parser.add_argument(
    '-mf', '--mapFile', help='Specifies the mapping file for HGNC to Ensembl ID from another ETL I made', required=True)

args = vars(parser.parse_args())

mapID(gene_file=args['geneFile'], map_file=args['mapFile'], output=args['output'])
