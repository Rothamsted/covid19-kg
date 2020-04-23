from chembl_webresource_client.new_client import new_client
import pandas as pd
import progressbar as pb
import argparse, time

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

def mapData(db_dir, base_dir, release_ver):
    db_dir = f"{db_dir}/DrugBank.csv"
    drugbank_df = pd.read_csv(db_dir)
    drugbank_df = drugbank_df.drop_duplicates()
    drug_names = drugbank_df['Name'].tolist()

    chembl_id_list,drugbank_list_updated  = [], []
    drug_names = list(dict.fromkeys(drug_names))
    chunks = (len(drug_names) - 1) // 50 + 1 # Chunk into 50 to be handled by EBI as fast as possible.

    widgets = ['Test: ', pb.Percentage(), ' ', pb.Bar(marker='0',left='[',right=']'), ' ', pb.ETA(), ' ', pb.FileTransferSpeed()]
    pbar = pb.ProgressBar(widgets=widgets, maxval=chunks)
    cv_start_time = time()
    pbar.start()

    for i in range(chunks):
        batch = drug_names[i*50:(i+1)*50]
        for drug in batch:
            try:
                chembl_id_list.append(new_client.molecule.search(drug)[0]['molecule_chembl_id'])
                drugbank_list_updated.append(drug)
                print(drug)
            except:
                pass
        time.sleep(0.5) # limit throttling on EBI network, if present.
        pbar.update(i)

    chembl_to_db_json =  {
                            'Name': drugbank_list_updated,
                            'ChEMBL ID': chembl_id_list
                         }

    chembl_to_db_pd = pd.DataFrame(chembl_to_db_json, columns=['Name', 'ChEMBL ID'])
    chembl_to_db_pd_updated = chembl_to_db_pd.merge(drugbank_df, on=['Name', 'Name'], how='inner')
    chembl_to_db_pd_updated = chembl_to_db_pd_updated[['DrugBank', 'ChEMBL ID', 'Name']]
    chembl_to_db_pd_updated.to_csv(f"{base_dir}{release_ver}/humanKnet/organisms/homo_sapiens/drug/db_chembl_mapping.tsv", sep="\t", index=None, header=True)
    print(f"\nData written out to {base_dir}{release_ver}/humanKnet/organisms/homo_sapiens/drug/db_chembl_mapping.tsv\nTask completed in {timeFormat(time()-cv_start_time)}\n")

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
