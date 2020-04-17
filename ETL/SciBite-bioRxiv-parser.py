import json, glob, os, argparse, re
import pandas as pd
import numpy as np
from pandas.io.json import json_normalize
from collections import OrderedDict
import requests
import math


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

def combineDFs(target, data, abstract_cols, all_cols, abstract_df, all_df, pmc):
    if pmc is None:
        abstract_df = sciBiteDF(df=abstract_df, cols=abstract_cols, target=target)
    all_df = sciBiteDF(df=all_df, cols=all_cols, target=target)
    meta_df = sciMetaDF(target=target, data=data)
    if pmc is None:
        concat_df = pd.concat([abstract_df, all_df, meta_df], ignore_index=True, sort=False)
    else:
        concat_df = pd.concat([all_df, meta_df], ignore_index=True, sort=False)
    concat_df['id'] = concat_df.id.apply(str).str.replace("{", '').str.replace("}", '').str.replace('\'', '').str.strip()
    normalized_concat_df = concat_df['id'].groupby([concat_df.name, concat_df.text, concat_df.hit_count]).apply(set).reset_index()
    df_agg_concat_df = normalized_concat_df.groupby(normalized_concat_df['name']).aggregate({'name': 'first', 'id': 'first', 'hit_count': 'sum', 'text': ''.join})
    df_agg_concat_df.reset_index(drop=True, inplace=True)
    # if pmc is not None:
    #     response = requests.get(f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={data['paper_id']}&format=json&tool=my_tool")
    #     df_agg_concat_df['paper_id'] = response.json()['records'][0]['pmid']
    # else:
    df_agg_concat_df['paper_id'] = data['paper_id']
    return df_agg_concat_df

def addToList(type, list_, target, data, abstract_cols, all_cols, abstract_df, all_df, pmc):
    try:
        if type is not None:
            list_.append(combineDFs(target, data, abstract_cols, all_cols, abstract_df, all_df, pmc))
    except:
        pass

def convertData(output_dir, base_dir, all_output, sarscov, genes, cvprot, drug, hpo, pmc):
    print(f"Making directory in {output_dir}/output/")
    if not os.path.exists(f"{output_dir}/output/"):
        os.makedirs(f"{output_dir}/output/")


    json_file_dirs = glob.glob(f"{base_dir}/*.json")
    if pmc is None:
        paper_id_list = [x.split("/")[-1].replace(".json", "") for x in json_file_dirs]
    else:
        paper_id_list = [x.split("/")[-1].replace(".xml.json", "") for x in json_file_dirs]
    #index = [x for x, s in enumerate(json_file_dirs) if '53442eacc3f233078507fa37b78267399e8c1e3b' in s ]
    pmc_list, pd_list, SARSCORS_list, GENE_list, CVPROT_list, DRUG_list, HPO_list = [], [], [], [], [], [], []
    for file_no, file_dir in enumerate(json_file_dirs):
        with open(file_dir, 'r') as f:
            allData = json.load(f)


        # Obtain all SciBite termite data
        if pmc is None:
            termite_abstract_df, termite_all_df = json_normalize(allData['abstract']), json_normalize(allData['body_text'])
            termite_abstract_cols, termite_all_df_cols = termite_abstract_df.columns, termite_all_df.columns
        else:
            termite_abstract_df, termite_all_df = None, json_normalize(allData['body_text'])
            termite_abstract_cols, termite_all_df_cols = None, termite_all_df.columns
        # SARS ontologies

        addToList(type=sarscov, list_=SARSCORS_list, target="SARSCOV",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df, pmc=pmc)
        # GENE
        addToList(type=genes, list_=GENE_list, target="GENE",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df,pmc=pmc)
        # CVPROT
        addToList(type=cvprot, list_=CVPROT_list, target="CVPROT",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df,pmc=pmc)
        #DRUG
        addToList(type=drug, list_=DRUG_list, target="DRUG",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df,pmc=pmc)
        #HPO
        addToList(type=hpo, list_=HPO_list, target="HPO",data=allData, abstract_cols=termite_abstract_cols, all_cols=termite_all_df_cols, abstract_df=termite_abstract_df, all_df=termite_all_df,pmc=pmc)

        if pmc is None:
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
        else:
            paper_ids= []
            # Extract paper metadata
            paper_ids.append(allData['paper_id'])
            #nlm_ids.append(requests.get(f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={allData['paper_id']}&format=json&tool=my_tool"))
            paper_data = {
                            'PMC': allData['paper_id'],
                        }
            paper_df = pd.DataFrame([paper_data], columns=['PMC'])
            pd_list.append(paper_df)

    if all_output is not None and pmc is not None:
        print(f"Concatenating data\n\n")
        concat_df = pd.concat(pd_list, axis=0, ignore_index=False)
        #concat_df=concat_df[concat_df.PMC.apply(lambda s: np.any([t in s for t in paper_id_list]))]
        concat_df=concat_df.dropna().drop_duplicates(['PMC'])
        pmc_ids = concat_df.PMC.tolist()
        pmc_ids = str(pmc_ids)[1:-1].replace("'", "").replace(", ", ",")
        pmd_list = pmc_ids.split(",")
        if len(pmd_list) > 199:
            response_chunk_list = []
            remove_ids = []
            divisor = len(pmd_list)/199
            chunk = np.array_split(pmd_list, math.ceil(divisor))
            for i in range(len(chunk)):
                id_chunks = str(chunk[i].tolist())[1:-1].replace("'", "").replace(", ", ",")
                response = requests.get(f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={id_chunks}&format=json&tool=my_tool")
                pmid_list = []
                for d in response.json()['records']:
                    try:
                        response_chunk_list.append(d['pmid'])
                    except:
                        remove_ids.append(d['pmcid'])

            concat_df = concat_df[~concat_df['PMC'].isin(remove_ids)] # Remove ID's which weren't matched by the API
            concat_df['PMID'] = response_chunk_list
            pmid_mapping_df = concat_df
            pmid_mapping_df.columns = ["paper_id", "PMID"]
            pmid_mapping_df.to_csv(f"{output_dir}/output/{all_output}.tsv", sep="\t", index=None, header=True)
        else:
            concat_df['PMID'] = pmid_list
            pmid_mapping_df = concat_df
            pmid_mapping_df.columns = ["paper_id", "PMID"]
            pmid_mapping_df.to_csv(f"{output_dir}/output/{all_output}.tsv", sep="\t", index=None, header=True)
        print(f"Writing data out in {base_dir}/output/covid19_papers.tsv\n\n")

    if all_output is not None and pmc is None:
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
            sars_concat['id'] = sars_concat['id'].str.replace('{', '').str.replace('\'}', '').str.split("/").str[-1].str.replace("_", ": ").str.replace("'", "")
            sars_concat=sars_concat[sars_concat.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
            sars_concat['text'] = sars_concat['text'].str.strip().replace('', np.nan)
            sars_concat=sars_concat.dropna().drop_duplicates(['name','id','paper_id'])
            sars_concat = sars_concat[sars_concat['paper_id'] != '']
            if pmc is None:
                sars_concat = sars_concat[['name', 'id', 'paper_id', 'hit_count', 'text']]
                sars_concat.to_csv(f"{output_dir}/output/{sarscov}.tsv", sep="\t", index=None, header=True)
            else:
                sars_concat_merged = pmid_mapping_df.merge(sars_concat, on=['paper_id', 'paper_id'], how = 'inner')
                sars_concat_merged = sars_concat_merged.dropna().drop_duplicates()
                sars_concat_merged.drop('paper_id', axis=1, inplace=True)
                sars_concat_merged = sars_concat_merged[['name', 'id', 'PMID', 'hit_count', 'text']]
                sars_concat_merged.to_csv(f"{output_dir}/output/{sarscov}.tsv", sep="\t", index=None, header=True)
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
            genes_concat['id'] = genes_concat['id'].str.replace('{', '').str.replace('\'}', '').str.split("=").str[-1].str.replace("'", "")
            genes_concat['text'] = genes_concat['text'].replace('', np.nan)
            genes_concat = genes_concat.dropna().drop_duplicates(['name','id','paper_id'])
            genes_concat = genes_concat[genes_concat.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
            genes_concat = genes_concat[genes_concat['paper_id']  != '']
            if pmc is None:
                genes_concat = genes_concat[['name', 'id', 'paper_id', 'hit_count', 'text']]
                genes_concat.to_csv(f"{output_dir}/output/{genes}.tsv", sep="\t", index=None, header=True)
            else:
                genes_concat_merged = pmid_mapping_df.merge(genes_concat, on=['paper_id', 'paper_id'], how = 'inner')
                genes_concat_merged = genes_concat_merged.dropna().drop_duplicates()
                genes_concat_merged.drop('paper_id', axis=1, inplace=True)
                genes_concat_merged = genes_concat_merged[['name', 'id', 'PMID', 'hit_count', 'text']]
                genes_concat_merged.to_csv(f"{output_dir}/output/{genes}.tsv", sep="\t", index=None, header=True)
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
            if pmc is None:
                cvprot_list = cvprot_list[['name', 'id', 'paper_id', 'hit_count', 'text']]
                cvprot_list.to_csv(f"{output_dir}/output/{cvprot}.tsv", sep="\t", index=None, header=True)
            else:
                cvprot_list_merged = pmid_mapping_df.merge(cvprot_list, on=['paper_id', 'paper_id'], how = 'inner')
                cvprot_list_merged = cvprot_list_merged.dropna().drop_duplicates()
                cvprot_list_merged.drop('paper_id', axis=1, inplace=True)
                cvprot_list_merged = cvprot_list_merged[['name', 'id', 'PMID', 'hit_count', 'text']]
                cvprot_list_merged.to_csv(f"{output_dir}/output/{cvprot}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting CVPROT data!\n")
        except Exception as e:
            print(f"CVPROT error {e}")
            pass

    if drug is not None:
        try:
            print(f"Concatenating drug data\n\n")
            drug_list = pd.concat(DRUG_list, axis=0, ignore_index=True, sort=False)
            drug_list['id'] = drug_list.id.apply(str)
            drug_list['id'] = drug_list['id'].str.replace('{', '').str.replace('\'}', '').str.split("/").str[-1]
            drug_list = drug_list[drug_list.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
            drug_list['text'] = drug_list['text'].str.strip().replace('', np.nan)
            drug_list=drug_list.dropna().drop_duplicates(['name','id','paper_id'])
            drug_list = drug_list[drug_list['paper_id'] != '']
            if pmc is None:
                drug_list = drug_list[['name', 'id', 'paper_id', 'hit_count', 'text']]
                drug_list.to_csv(f"{output_dir}/output/{drug}.tsv", sep="\t", index=None, header=True)
            else:
                drug_list_merged = pmid_mapping_df.merge(drug_list, on=['paper_id', 'paper_id'], how = 'inner')
                drug_list_merged = drug_list_merged.dropna().drop_duplicates()
                drug_list_merged.drop('paper_id', axis=1, inplace=True)
                drug_list_merged = drug_list_merged[['name', 'id', 'PMID', 'hit_count', 'text']]
                drug_list_merged.to_csv(f"{output_dir}/output/{drug}.tsv", sep="\t", index=None, header=True)
            print("Finished outputting drug data!\n")
        except Exception as e:
            print(f"DRUG error {e}")
            pass

    if hpo is not None:
        try:
            print(f"Concatenating hpo data\n\n")
            hpo_list = pd.concat(HPO_list, axis=0, ignore_index=True, sort=False)
            hpo_list['id'] = hpo_list.id.apply(str)
            hpo_list['id'] = hpo_list['id'].str.replace('{', '').str.replace('\'}', '').str.split("/").str[-1].str.replace("HP", "HP:").str.replace("_", "")
            hpo_list[hpo_list.paper_id.apply(lambda s: np.any([t in s for t in paper_id_list]))]
            hpo_list['text'].str.strip().replace('', np.nan, inplace=True)
            hpo_list=hpo_list.dropna().drop_duplicates(['name','id','paper_id'])
            hpo_list = hpo_list[hpo_list['paper_id'] != '']
            if pmc is None:
                hpo_list = hpo_list[['name', 'id', 'paper_id', 'hit_count', 'text']]
                hpo_list.to_csv(f"{output_dir}/output/{hpo}.tsv", sep="\t", index=None, header=True)
            else:
                hpo_list_merged = pmid_mapping_df.merge(hpo_list, on=['paper_id', 'paper_id'], how = 'inner')
                hpo_list_merged = hpo_list_merged.dropna().drop_duplicates()
                hpo_list_merged.drop('paper_id', axis=1, inplace=True)
                hpo_list_merged = hpo_list_merged[['name', 'id', 'PMID', 'hit_count', 'text']]
                hpo_list_merged.to_csv(f"{output_dir}/output/{hpo}.tsv", sep="\t", index=None, header=True)
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
parser.add_argument(
    '-p', '--pmc', help='Arg for PMC', required=False)
args = vars(parser.parse_args())

convertData(output_dir=args['output'], base_dir=args['basedir'], all_output=args['alloutput'], sarscov=args['sarscov'], genes=args['genes'], cvprot=args['cvprot'], drug=args['drug'], hpo=args['hpo'], pmc=args['pmc'])
