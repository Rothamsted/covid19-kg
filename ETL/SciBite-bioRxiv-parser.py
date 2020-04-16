import json, glob, os, argparse, re
import pandas as pd
import numpy as np
from pandas.io.json import json_normalize
from collections import OrderedDict


def sciBiteDF(df, cols, target):
    if f"termite_hits.{target}" in cols:
        list_df =  df[[f'termite_hits.{target}', 'text']].dropna()
        # Ensure there is data
        if [x for x in df[f'termite_hits.{target}'].dropna() if x != []] != []:
            # Add the text to the JSON, to then normalize and put into a df
            for x, i in enumerate(list_df[f'termite_hits.{target}']):
                curr_ind = list_df[f'termite_hits.{target}'].index.tolist()[x]
                curr_text = list_df['text'][curr_ind]
                for key in i:
                    try:
                        sentence_loc = key['hit_sentence_locations'][0]
                        key['text'] = curr_text[sentence_loc[0]:sentence_loc[1]] + "<br>"
                    except:
                        key['text'] = curr_text + "<br>"

            # Normalize the df - must iterate through it and concat the normalized data to create the df
            normalized_df = pd.concat([pd.DataFrame(json_normalize(x)) for x in df[f'termite_hits.{target}'].dropna()], ignore_index=False)
            normalized_df['text'] = normalized_df['text'].str.replace("\t", "", regex=True).str.lstrip(".")
            normalized_df_ = normalized_df['id'].groupby([normalized_df.name, normalized_df.text, normalized_df.hit_count]).apply(set).reset_index()
            normalized_df_ = normalized_df_ [['name', 'id', 'text', 'hit_count']]
            df_agg = normalized_df_.groupby(normalized_df_['name']).aggregate({'name': 'first', 'id': 'first', 'hit_count': 'sum', 'text': ''.join})
            df_agg.reset_index(drop=True, inplace=True)
            return df_agg
        else:
            return None
    else:
        return None


def sciMetaDF(target, data):
    if target in data['metadata']['termite_hits']:
        meta_SARSCOV_pd = json_normalize(data['metadata']['termite_hits'][target])
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
    concat_df['id'] = concat_df.id.apply(str).str.replace("{", '').str.replace("}", '').str.replace('\'', '').str.strip()
    normalized_concat_df = concat_df['id'].groupby([concat_df.name, concat_df.text, concat_df.hit_count]).apply(set).reset_index()
    df_agg_concat_df = normalized_concat_df.groupby(normalized_concat_df['name']).aggregate({'name': 'first', 'id': 'first', 'hit_count': 'sum', 'text': ''.join})
    df_agg_concat_df.reset_index(drop=True, inplace=True)
    df_agg_concat_df['paper_id'] = data['paper_id']
    return df_agg_concat_df

def addToList(type, list_, target, data, abstract_cols, all_cols, abstract_df, all_df):
    try:
        if type is not None:
            list_.append(combineDFs(target, data, abstract_cols, all_cols, abstract_df, all_df))
    except:
        pass

def convertData(output_dir, base_dir, all_output, sarscov, genes, cvprot, drug, hpo):
    print(f"Making directory in {output_dir}/output/")
    if not os.path.exists(f"{output_dir}/output/"):
        os.makedirs(f"{output_dir}/output/")

    json_file_dirs = glob.glob(f"{base_dir}/*.json")
    paper_id_list = [x.split("/")[-1].replace(".json", "") for x in json_file_dirs]
    #index = [x for x, s in enumerate(json_file_dirs) if '53442eacc3f233078507fa37b78267399e8c1e3b' in s ]
    pd_list, SARSCORS_list, GENE_list, CVPROT_list, DRUG_list, HPO_list = [], [], [], [], [], []
    for file_no, file_dir in enumerate(json_file_dirs):
        with open(file_dir, 'r') as f:
            allData = json.load(f)

        # Obtain all SciBite termite data
        termite_abstract_df, termite_all_df = json_normalize(allData['abstract']), json_normalize(allData['body_text'])
        termite_abstract_cols, termite_all_df_cols = termite_abstract_df.columns, termite_all_df.columns
        # SARS ontologies

        addToList(type=sarscov, list_=SARSCORS_list, target="SARSCOV",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df)
        # GENE
        addToList(type=genes, list_=GENE_list, target="GENE",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df)
        # CVPROT
        addToList(type=cvprot, list_=CVPROT_list, target="CVPROT",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df)
        #DRUG
        addToList(type=drug, list_=DRUG_list, target="DRUG",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df)
        #HPO
        addToList(type=hpo, list_=HPO_list, target="HPO",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df)

        if all_output is not None:
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
                abstract = json_normalize(allData['abstract'])['text'].tolist()
                abstract_tidy = [x.replace("'", "") for x in abstract]
                abstract_tidy = "\n".join(abstract_tidy)
                # Extract paper names
                if 'first' in json_normalize(allData['metadata']['authors']):
                    for i, first_name in enumerate((json_normalize(allData['metadata']['authors'])['first']).tolist()):
                        name = f"{first_name} {json_normalize(allData['metadata']['authors'])['last'].tolist()[i]}"
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

    if all_output is not None:
        print(f"Concatenating data\n\n")
        concat_df = pd.concat(pd_list, axis=0, ignore_index=True)
        concat_df=concat_df[concat_df.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
        concat_df=concat_df.dropna()
        concat_df=concat_df.drop_duplicates(['AUTHORS','AbstractHeader','paper_id'])
        concat_df.replace('', "Not found", inplace=True)
        concat_df['Abstract'] = concat_df['Abstract'].str.strip()
        concat_df['AbstractHeader'] = concat_df['AbstractHeader'].str.strip()
        concat_df = concat_df[concat_df['paper_id'] != '']
        concat_df = concat_df[concat_df['Abstract'] != 'Not found']
        print(f"Writing data out in {base_dir}/output/covid19_papers.tsv\n\n")
        concat_df.to_csv(f"{output_dir}/output/{all_output}.tsv", sep="\t", index=None, header=True)
        print("Finished outputting main data!\n")

    if sarscov is not None:
        try:
            print(f"Concatenating sars-CoV data\n\n")
            sars_concat = pd.concat(SARSCORS_list, axis=0, ignore_index=True, sort=False)
            sars_concat = sars_concat[['name', 'id', 'hit_count', 'text', 'paper_id']]
            sars_concat=sars_concat.dropna()
            #sars_concat = sars_concat.replace(np.nan, 'N/A', regex=True)
            sars_concat['id'] = sars_concat.id.apply(str)
            sars_concat['id'] = sars_concat['id'].str.replace('{', '').str.replace('\'}', '').str.split("/").str[-1].str.replace("_", ": ")
            sars_concat=sars_concat[sars_concat.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
            sars_concat['text'] = sars_concat['text'].str.strip().replace('', np.nan)
            sars_concat=sars_concat.dropna().drop_duplicates(['name','id','paper_id'])
            sars_concat = sars_concat[sars_concat['paper_id'] != '']
            sars_concat = sars_concat[['name', 'id', 'paper_id', 'hit_count', 'text']]
            sars_concat.to_csv(f"{output_dir}/output/{sarscov}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting SARS_CoV data!\n")
        except Exception as e:
            print(f"SARS error {e}")
            pass

    if genes is not None:
        try:
            print(f"Concatenating gene data\n\n")
            genes_concat = pd.concat(GENE_list, axis=0, ignore_index=True, sort=False)
            genes_concat = genes_concat[['name', 'id', 'hit_count', 'text', 'paper_id']]
            genes_concat['id'] = genes_concat.id.apply(str)
            genes_concat['id'] = genes_concat['id'].str.replace('{', '').str.replace('\'}', '').str.split("=").str[-1]
            genes_concat['text'] = genes_concat['text'].replace('', np.nan)
            genes_concat = genes_concat.dropna().drop_duplicates(['name','id','paper_id'])
            genes_concat = genes_concat[genes_concat.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
            genes_concat = genes_concat[genes_concat['paper_id']  != '']
            genes_concat = genes_concat[['name', 'id', 'paper_id', 'hit_count', 'text']]
            genes_concat.to_csv(f"{output_dir}/output/{genes}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting Gene data!\n")
        except Exception as e:
            print(f"GENES error {e}")
            pass

    if cvprot is not None:
        try:
            print(f"Concatenating cvprot data\n\n")
            cvprot_list = pd.concat(CVPROT_list, axis=0, ignore_index=True, sort=False)
            cvprot_list = cvprot_list[['name', 'id', 'hit_count', 'text', 'paper_id']]
            cvprot_list['id'] = cvprot_list.id.apply(str)
            cvprot_list['id'] = cvprot_list['id'].str.replace('{', '').str.replace('\'}', '').str.split("/").str[-1]
            cvprot_list = cvprot_list[cvprot_list.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
            cvprot_list['text'] = cvprot_list['text'].str.strip().replace('', np.nan)
            cvprot_list=cvprot_list.dropna().drop_duplicates(['name','id','paper_id'])
            cvprot_list = cvprot_list[cvprot_list['paper_id'] != '']
            cvprot_list = cvprot_list[['name', 'id', 'paper_id', 'hit_count', 'text']]
            cvprot_list.to_csv(f"{output_dir}/output/{cvprot}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting CVPROT data!\n")
        except Exception as e:
            print(f"CVPROT error {e}")
            pass

    if drug is not None:
        try:
            print(f"Concatenating drug data\n\n")
            drug_list = pd.concat(DRUG_list, axis=0, ignore_index=True, sort=False)
            drug_list = drug_list[['name', 'id', 'hit_count', 'text', 'paper_id']]
            drug_list['id'] = drug_list.id.apply(str)
            drug_list['id'] = drug_list['id'].str.replace('{', '').str.replace('\'}', '').str.split("/").str[-1]
            drug_list = drug_list[drug_list.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
            drug_list['text'] = drug_list['text'].str.strip().replace('', np.nan)
            drug_list=drug_list.dropna().drop_duplicates(['name','id','paper_id'])
            drug_list = drug_list[drug_list['paper_id'] != '']
            drug_list = drug_list[['name', 'id', 'paper_id', 'hit_count', 'text']]
            drug_list.to_csv(f"{output_dir}/output/{drug}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting drug data!\n")
        except Exception as e:
            print(f"DRUG error {e}")
            pass

    if hpo is not None:
        try:
            print(f"Concatenating hpo data\n\n")
            hpo_list = pd.concat(HPO_list, axis=0, ignore_index=True, sort=False)
            hpo_list = hpo_list[['name', 'id', 'hit_count', 'text', 'paper_id']]
            hpo_list['id'] = hpo_list.id.apply(str)
            hpo_list['id'] = hpo_list['id'].str.replace('{', '').str.replace('\'}', '').str.split("/").str[-1].str.replace("HP", "HP:")
            hpo_list[hpo_list.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
            hpo_list['text'].str.strip().replace('', np.nan, inplace=True)
            hpo_list=hpo_list.dropna().drop_duplicates(['name','id','paper_id'])
            hpo_list = hpo_list[hpo_list['paper_id'] != '']
            hpo_list = hpo_list[['name', 'id', 'paper_id', 'hit_count', 'text']]
            hpo_list.to_csv(f"{output_dir}/output/{hpo}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting hpo data!\n")
        except Exception as e:
            print(f"HPO error {e}")
            pass



parser = argparse.ArgumentParser(
    description="This script will convert Kaggle data from JSON format to table using metadata only.\n")
parser.add_argument(
    '-bdir', '--basedir', help='Specifies the base directory you are using from where you will find the JSON data.', required=True)
parser.add_argument(
    '-all', '--alloutput', help='Specifies the output name for all of the data.', required=False)
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
