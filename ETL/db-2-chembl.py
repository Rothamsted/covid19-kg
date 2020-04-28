from chembl_webresource_client.new_client import new_client
import pandas as pd
import progressbar as pb
import argparse, time, os

def timeFormat(secs):
    """ Formats Seconds to relevant time interval"""

    secs = float("{0:.2f}".format(secs))
    days = secs//86400
    hours = (secs - days*86400)//3600
    minutes = (secs - days*86400 - hours*3600)//60
    seconds = secs - days*86400 - hours*3600 - minutes*60
    time = ("{0} day{1}, ".format(days, "s" if days!=1 else "") if days else "") + \
    ("{0} hour{1} ".format(hours, "s" if hours!=1 else "") if hours else "") + \
    ("{0} minute{1} ".format(minutes, "s" if minutes!=1 else "") if minutes else "") + \
    ("{0} second{1} ".format(seconds, "s" if seconds!=1 else "") if seconds else "")
    return time


def searchChEMBL(chunks, drugname,chembl_id_list, drugbank_list_updated):

    widgets = ['Test: ', pb.Percentage(), ' ', pb.Bar(marker='0',left='[',right=']'), ' ', pb.ETA(), ' ', pb.FileTransferSpeed()]
    pbar = pb.ProgressBar(widgets=widgets, maxval=chunks)
    pbar.start()
    j = None
    count = 0
    for i in range(chunks):
        batch = drug_names[i*50:(i+1)*50]
        for drug in batch:
            try:
                if drug != j:
                    chembl_id_list.append(new_client.molecule.search(drug)[0]['molecule_hierarchy']['parent_chembl_id'])
                    drugbank_list_updated.append(drug)
                    j = drug
                    count+=1
            except:
                try:
                    if new_client.molecule.search(drug)[0]['molecule_chembl_id'] != '[]':
                        chembl_id_list.append(new_client.molecule.search(drug)[0]['molecule_chembl_id'])
                        drugbank_list_updated.append(drug)
                except:
                    pass
        time.sleep(0.1)
        pbar.update(i)


def mapData(db_dir, base_dir, release_ver):
    db_dir = f"{db_dir}/DrugBank.csv"
    drugbank_df = pd.read_csv(db_dir)
    drugbank_df = drugbank_df.drop_duplicates()
    # Functions
    loc_names = lambda x,y: drugbank_df.iloc[x:y, :]['Name'].tolist()
    chunk_it = lambda x:  (len(x) - 1) // 50 + 1

    drug_names_1, drug_names_2 = loc_names(0,3391), loc_names(3391,6782)
    drug_names_3, drug_names_4 = loc_names(6782,10172), loc_names(10172, 13563)

    chembl_id_list_1,drugbank_list_updated_1  = [], []
    chembl_id_list_2,drugbank_list_updated_2  = [], []
    chembl_id_list_3,drugbank_list_updated_3  = [], []
    chembl_id_list_4,drugbank_list_updated_4  = [], []

    chunks_1, chunks_2 = chunk_it(drug_names_1), chunk_it(drug_names_2)
    chunks_3, chunks_4 = chunk_it(drug_names_3), chunk_it(drug_names_4)
    cv_start_time = time.time()

    searchChEMBL(chunks_1, drug_names_1,chembl_id_list_1, drugbank_list_updated_1)
    print("Finished batch 1 of 4\n")
    searchChEMBL(chunks_2, drug_names_2,chembl_id_list_2, drugbank_list_updated_2)
    print("Finished batch 2 of 4\n")
    searchChEMBL(chunks_3, drug_names_3,chembl_id_list_3, drugbank_list_updated_3)
    print("Finished batch 3 of 4\n")
    searchChEMBL(chunks_4, drug_names_4,chembl_id_list_4, drugbank_list_updated_4)
    print("Finished batch 4 of 4\n")

    chembl_id_list_ = chembl_id_list_1 + chembl_id_list_2 + chembl_id_list_3 + chembl_id_list_4
    drugbank_list_updated_ = drugbank_list_updated_1 + drugbank_list_updated_2 + drugbank_list_updated_3 + drugbank_list_updated_4

    chembl_to_db_json =  {
                            'Name': drugbank_list_updated_,
                            'ChEMBL ID': chembl_id_list_
                         }

    chembl_to_db_pd = pd.DataFrame(chembl_to_db_json, columns=['Name', 'ChEMBL ID'])
    chembl_to_db_pd_updated = chembl_to_db_pd.merge(drugbank_df, on=['Name', 'Name'], how='inner')
    chembl_to_db_pd_updated = chembl_to_db_pd_updated[['DrugBank ID', 'ChEMBL ID', 'Name']]
    chembl_to_db_pd_updated.columns = ['DrugBank ID', 'ChEMBL ID', 'Drug Name']
    chembl_to_db_pd_updated = chembl_to_db_pd_updated.drop_duplicates()

    chembl_to_db_pd_updated.to_csv(f"{base_dir}{release_ver}/humanKnet/organisms/homo_sapiens/drug/db_chembl_mapping.tsv", sep="\t", index=None, header=True)
    print(f"\nData written out to {base_dir}{release_ver}/humanKnet/organisms/homo_sapiens/drug/db_chembl_mapping.tsv\nTask completed in {timeFormat(time.time()-cv_start_time)}\n")

parser = argparse.ArgumentParser(
    description="This script will convert DrugBank ID's to ChEMBL ID's\n")
parser.add_argument(
    '-db', '--drugbank', help='Specifies the drugbank directory you are using from where you will find the DrugBank csv data containg all drug data', required=True)
parser.add_argument(
    '-bdir', '--basedir', help='Specifies the base directory for the data', required=True)
parser.add_argument(
    '-r', '--releasever', help='Specifies the Ensembl release version (Coherent with RRes ETL file directories nomenclature)', required=True)
args = vars(parser.parse_args())
mapData(db_dir=args['drugbank'], base_dir=args['basedir'], release_ver=args['releasever'])
