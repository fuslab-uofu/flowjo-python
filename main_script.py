# This is a sample Python script.
import import_tools
import pickle
import pandas as pd
from get_remote_data import get_remote_data


def define_paths():

    global server_json       # Needed for ssh file transfer only
    global remotepath        # Location of remote data (CSVs)
    global remotepath_mrti   # Loation of 'tempData.csv' for MRTI data
    global datapath          # Location to store local data (CSVs)
    global fcspath           # Location to store local FCS data (CSVs)
    global exportpath        # Location to store the generated flowtree Python Objects

    server_json = '/Users/sarajohnson/Library/CloudStorage/Box-Box/Python Projects/flowjo_export_analysis/CSV Data/server_login.json'
    remotepath = '/v/raid10/users/sjohnson/Data/Flow Cytometry Data/IACUC 21-11013/Data Export/'
    remotepath_mrti = '/v/raid10/users/sjohnson/Experiment Analysis/IACUC21-11013/'
    datapath = '/Users/sarajohnson/Library/CloudStorage/Box-Box/Python Projects/flowjo_export_analysis/CSV Data/'
    fcspath = datapath + 'FCS Data/'
    exportpath = '/Users/sarajohnson/Library/CloudStorage/Box-Box/Python Projects/flowjo_export_analysis/ExportedData/'


def main(get_remote, toss_files):

    define_paths()

    # ----- LOAD CSV FILES -------

    # Transfer data to local machine over SFTP connection
    if get_remote:
        get_remote_data(remotepath, fcspath, server_json)

    # Load CSVs and concatenate into dataframe
    print('Loading CSVs for Tumor data: ')
    df_tumor, df_tumor_nans, mdh_tumor = import_tools.preprocess_csvs(fcspath,
                                                                      search_string='Tumor',
                                                                      exclude_files=toss_files)
    print('Loading CSVs for Spleen data: ')
    df_spleen, df_spleen_nans, mdh_spleen = import_tools.preprocess_csvs(fcspath,
                                                                         search_string='Spleen',
                                                                         exclude_files=toss_files)

    # ----- CREATE FLOWTREES -------

    # Create population flowtree for TUMOR tissue dataframe
    print('Creating population tree for Tumor data: ')
    csv_popnames_master = datapath + 'Tumor Population Names.csv'
    tumor_root = import_tools.create_pop_tree(csv_popnames_master,
                                              df_tumor,
                                              tissue_type='tumor',
                                              show_tree=False)

    # Calculate Counts for TUMOR tree
    tumor_root.calculate_counts_tree()
    tumor_root.show_pop_name_tree()  # display tree with population nicknames

    # Create population flowtree for SPLEEN tissue dataframe
    print('Creating population tree for Spleen data: ')
    csv_popnames_master = datapath + 'Spleen Population Names.csv'
    spleen_root = import_tools.create_pop_tree(csv_popnames_master,
                                               df_spleen,
                                               tissue_type='spleen',
                                               show_tree=False)

    # Calculate Counts for SPLEEN tree
    spleen_root.calculate_counts_tree()
    spleen_root.show_pop_name_tree()  # display tree with population nicknames

    # ----- ADD MRTI DATA -------

    # Transfer data to local machine over SFTP connection
    if get_remote:
        get_remote_data(remotepath_mrti, datapath, '.csv')

    mrtidata = pd.read_csv(datapath + 'tempData.csv')
    mrtidata.rename(columns={'MouseID': 'SampleID'}, inplace=True)
    mrtidata.drop(['Group', 'AblationType'], axis=1, inplace=True)

    # Create new column for voxels above 50-degrees MTP
    mrtidata.insert(3, 'heated50', mrtidata['heated5060'] + mrtidata['heated60'])

    # Merge MRTI data with FCS data and add to flowtree root
    _ = tumor_root.append_mrti_data(mrtidata)
    _ = spleen_root.append_mrti_data(mrtidata)

    # ----- SAVE FLOWTREES -------

    # Save flowtrees to local path
    with open(exportpath + 'tumor_tree_master.pkl', 'wb') as f:
        pickle.dump(tumor_root, f)

    with open(exportpath + 'spleen_tree_master.pkl', 'wb') as f:
        pickle.dump(spleen_root, f)


if __name__ == "__main__":

    toss_files = ['Spleen Table - C3 -Column Names.csv',
                  'Spleen Population Names.csv',
                  'Tumor Population Names.csv']

    get_remote = True

    main(get_remote, toss_files)
