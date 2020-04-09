import json, glob, os, argparse
import pandas as pd
import numpy as np
from pandas.io.json import json_normalize
from collections import OrderedDict

#base_dir='/home/joseph/covid19/biorxiv_medrxiv'

def sciBiteDF(df, cols, target):

    if f"termite_hits.{target}" in cols:
        list = df[f'termite_hits.{target}'].dropna()
        df_index = df.index[df[f'termite_hits.{target}'].nonzero()].tolist()
        list_index = df[f'termite_hits.{target}']
        list = [x for x in list if x != []] # Remove empty lists
        df_ = json_normalize(pd.DataFrame(list)[0])
        text_list = df.text[df_index].dropna()
        df_['text'] = text_list.reset_index(drop=True)
        return df_.drop_duplicates()
    else:
        return None

def sciMetaDF(target, data):
    if target in data['metadata']['termite_hits']:
        meta_SARSCOV_pd = json_normalize(allData['metadata']['termite_hits'][target])
        if meta_SARSCOV_pd.empty:
            meta_SARSCOV_pd = None
        return meta_SARSCOV_pd
    else:
        return None

def combineDFs(target, data, abstract_cols, all_cols, abstract_df, all_df):
    abstract_df = sciBiteDF(df=abstract_df, cols=abstract_cols, target=target)
    all_df = sciBiteDF(df=all_df, cols=all_cols, target=target)
    meta_df = sciMetaDF(target=target, data=data)
    concat_df = pd.concat([abstract_df, all_df, meta_df], ignore_index=True, sort=False)
    concat_df['paper_id'] = data['paper_id']
    return concat_df


def convertData(output_dir, base_dir, all_output, sarscov, genes, cvprot, drug, hpo):
    print(f"Making directory in {output_dir}/output/")
    if not os.path.exists(f"{output_dir}/output/"):
        os.makedirs(f"{output_dir}/output/")

    json_file_dirs = glob.glob(f"{base_dir}/*.json")
    pd_list, SARSCORS_list, GENE_list, CVPROT_list, DRUG_list, HPO_list = [], [], [], [], [], []
    for file_no, file_dir in enumerate(json_file_dirs):
        with open(file_dir, 'r') as f:
            allData = json.load(f)

        # Obtain all SciBite termite data
        termite_abstract_df, termite_all_df = json_normalize(allData['abstract']), json_normalize(allData['body_text'])
        termite_abstract_cols, termite_all_df_cols = termite_abstract_df.columns, termite_all_df.columns

        text = termite_all_df.text
        # SARS ontologies
        try:
            SARSCORS_list.append(combineDFs(target="SARSCOV", data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df))
            print("Created SARSCORS list")
        except  Exception as e:
            print(f"failed because {e}")
            pass
        # GENE
        try:
            GENE_list.append((combineDFs(target="GENE", data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df)))
            print("Created GENES list")
        except  Exception as e:
            print(f"failed because {e}")
            pass
        # CVPROT
        try:
            CVPROT_list.append((combineDFs(target="CVPROT", data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df)))
            print("Created CVPROT list")
        except  Exception as e:
            print(f"failed because {e}")
            pass
        #DRUG
        try:
            DRUG_list.append((combineDFs(target="DRUG", data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df)))
            print("Created DRUGS list")
        except  Exception as e:
            print(f"failed because {e}")
            pass
        #HPO
        try:
            HPO_list.append((combineDFs(target="HPO", data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df)))
            print("Created HPO list")
        except  Exception as e:
            print(f"failed because {e}")
            pass



        paper_ids, paper_title, paper_authors = [], [], []
        # Extract paper metadata
        paper_ids.append(allData['paper_id'])

        try:
            paper_title.append(allData['metadata']['title'])
        except:
            try:
                paper_title.append(allData['abstract']['text'][0])
            except:
                paper_title.append("Title not given")
                pass
            pass
        try:
            abstract = json_normalize(allData['abstract'])['text']
            abstract_tidy = [x.replace("'", "") for x in list(abstract)]
            abstract_tidy = "\n".join(abstract_tidy)
            # Extract paper names
            if 'first' in json_normalize(allData['metadata']['authors']):
                for i, first_name in enumerate(list(json_normalize(allData['metadata']['authors'])['first'])):
                    name = f"{first_name} {list(json_normalize(allData['metadata']['authors'])['last'])[i]}"
                    paper_authors.append(name)
            else:
                paper_authors.append("Not found")

            paper_data = {
                            'Abstract': abstract_tidy,
                            'AUTHORS': str(paper_authors),
                            'AbstractHeader': paper_title,
                            'paper_id': paper_ids
                         }

            paper_df = pd.DataFrame(paper_data, columns=['paper_id', 'AUTHORS', 'AbstractHeader', 'Abstract'])
            pd_list.append(paper_df)
        except:
            pass

    print(f"Concatenating data\n\n")
    concat_df = pd.concat(pd_list, axis=0, ignore_index=True)
    print(f"Writing data out in {base_dir}/output/covid19_papers.tsv\n\n")
    concat_df.to_csv(f"{output_dir}/output/{all_output}.tsv", sep="\t", index=None, header=True)
    print("Finished outputting main data!\n")

    if sarscov is not None:
        try:
            print(f"Concatenating sars-CoV data\n\n")
            sars_concat = pd.concat(SARSCORS_list, axis=0, ignore_index=True, sort=False)
            sars_concat = sars_concat[['name', 'id', 'hit_count']]
            sars_concat = sars_concat.replace(np.nan, 'N/A', regex=True)
            sars_concat['id'] = sars_concat['id'].str.split("/").str[-1]
            sars_concat['id'] = sars_concat['id'].str.replace("_", ": ")
            sars_concat.to_csv(f"{base_dir}/output/{sarscov}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting SARS_CoV data!\n")
        except:
            pass

    if genes is not None:
        try:
            print(f"Concatenating gene data\n\n")
            genes_concat = pd.concat(GENE_list, axis=0, ignore_index=True, sort=False)
            genes_concat = genes_concat[['name', 'id' , 'paper_id']]
            genes_concat['id'] = genes_concat['id'].str.split("=").str[1]
            genes_concat.to_csv(f"{output_dir}/output/{genes}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting Gene data!\n")
        except:
            pass

    if cvprot is not None:
        try:
            print(f"Concatenating cvprot data\n\n")
            cvprot_list = pd.concat(CVPROT_list, axis=0, ignore_index=True, sort=False)
            cvprot_list = cvprot_list[['name', 'id' , 'paper_id']]
            cvprot_list['id'] = cvprot_list['id'].str.split("/").str[-1]
            cvprot_list.to_csv(f"{output_dir}/output/{cvprot}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting CVPROT data!\n")
        except:
            pass

    if drug is not None:
        try:
            print(f"Concatenating cvprot data\n\n")
            drug_list = pd.concat(DRUG_list, axis=0, ignore_index=True, sort=False)
            drug_list = drug_list[['name', 'id' , 'paper_id']]
            drug_list['id'] = drug_list['id'].str.split("/").str[-1]
            drug_list.to_csv(f"{output_dir}/output/{drug}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting drug data!\n")
        except:
            pass

    if hpo is not None:
        try:
            print(f"Concatenating cvprot data\n\n")
            hpo_list = pd.concat(HPO_list, axis=0, ignore_index=True, sort=False)
            hpo_list['id'] = hpo_list['id'].str.split("/").str[-1]
            hpo_list.to_csv(f"{output_dir}/output/{hpo}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting drug data!\n")
        except:
            pass



parser = argparse.ArgumentParser(
    description="This script will convert Kaggle data from JSON format to table using metadata only.\n")
parser.add_argument(
    '-bdir', '--basedir', help='Specifies the base directory you are using from where you will find the JSON data.', required=True)
parser.add_argument(
    '-all', '--alloutput', help='Specifies the output name for all of the data.', required=True)
parser.add_argument(
    '-o', '--output', help='Specifies the main output directory.', required=True)
parser.add_argument(
    '-scov', '--sarscov', help='Specifies the output name for SARSCOV data.', required=False)
parser.add_argument(
    '-g', '--genes', help='Specifies the output name for GENE data.', required=False)
parser.add_argument(
    '-cv', '--cvprot', help='Specifies the output name for CVPROT data.', required=False)
parser.add_argument(
    '-d', '--drug', help='Specifies the output name for DRUG data.', required=False)
parser.add_argument(
    '-hp', '--hpo', help='Specifies the output name for HPO data.', required=False)
args = vars(parser.parse_args())

convertData(output_dir=args['output'], base_dir=args['basedir'], all_output=args['alloutput'], sarscov=args['sarscov'], genes=args['genes'], cvprot=args['cvprot'], drug=args['drug'], hpo=args['hpo'])
