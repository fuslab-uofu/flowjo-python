import flowtree
import pickle



if __name__ == "__main__":

    global exportpath

    exportpath = '/Users/sarajohnson/Library/CloudStorage/Box-Box/Python Projects/flowjo_export_analysis/ExportedData/'
    tumor_tree_file = 'tumor_tree_master.pkl'
    spleen_tree_file = 'spleen_tree_master.pkl'

    # --- CHANGE THESE FOR EXPORT ----

    filename_suffix = 'subset'

    tumor_sub_pops = ['Tumor Cells', 'CD45+ Subset', 'T Cells', 'CD4+ T Cells', 'CD4+ Memory T Cells 2',
                      'CD8+ T Cells', 'CD8+ Memory T Cells 2', 'B Cells', 'MHCII+ B Cells', 'NK Cells',
                      'Macrophages', 'MHCII+ Macrophages', 'DCs', 'CD11b+ cDCs Migratory', 'Neutrophils', 'MDSCs']

    tumor_out_of_populations = ['CD45+ Subset', 'Live']

    spleen_sub_pops = ['CD45+ Subset', 'T Cells', 'CD4+ T Cells', 'CD4+ Memory T Cells 2',
                'CD8+ T Cells', 'CD8+ Memory T Cells 2', 'B Cells', 'Memory B Cells', 'NK Cells',
                'Macrophages', 'MHCII+ Macrophages', 'DCs', 'CD8+ cDCs', 'CD11b+ cDCs', 'Neutrophils', 'MDSCs']

    spleen_out_of_populations = ['CD45+ Subset', 'Live']


    # --- TUMOR FLOWTREE SUBSET ---

    with open(exportpath + tumor_tree_file, 'rb') as f:
        tumor_root = pickle.load(f)

    tumor_root.export_freqs_as_dataframe(tumor_sub_pops, tumor_out_of_populations,
                                         csv_filename=exportpath + 'Tumor_' + filename_suffix + '.csv',
                                         to_csv=True,
                                         merge_mrti=True,
                                         freq_of_parent=True)

    # --- SPLEEN FLOWTREE SUBSET ---

    with open(exportpath + spleen_tree_file, 'rb') as f:
        spleen_root = pickle.load(f)

    sub_pops = ['CD45+ Subset', 'T Cells', 'CD4+ T Cells', 'CD4+ Memory T Cells 2',
                'CD8+ T Cells', 'CD8+ Memory T Cells 2', 'B Cells', 'Memory B Cells', 'NK Cells',
                'Macrophages', 'MHCII+ Macrophages', 'DCs', 'CD8+ cDCs', 'CD11b+ cDCs', 'Neutrophils', 'MDSCs']

    spleen_root.export_freqs_as_dataframe(spleen_sub_pops, spleen_out_of_populations,
                                          csv_filename=exportpath + 'Spleen_' + filename_suffix + '.csv',
                                          to_csv=True,
                                          merge_mrti=True,
                                          freq_of_parent=True)

