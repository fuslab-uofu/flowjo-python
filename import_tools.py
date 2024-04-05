# This is a sample Python script.
import logging
import bigtree
import flowtree
import os
import csv
import pandas as pd


def load_csvs_to_dataframe(local_dir, search_string=None, exclude_files=None):
    """
    Loads CSV files in directory and concatenates into a DataFrame.

    Parameters
    ----------
    local_dir : string
        absolute path to CSV data on local machine.
    search_string : string
        includes all CSV filenames with this string.
    exclude_files : list[string]
        excludes filenames in list.

    Returns
    -------
    DataFrame
        the DataFrame of all concatenated CSV data

    """

    # Get list of all files in local directory
    files = os.listdir(local_dir)

    # Remove excluded_files from list
    files = [x for x in files if x not in exclude_files]

    # Keep files that contain the search_string
    if search_string:
        files = [x for x in files if search_string in x]

    print(files)

    # Loop through files, compare Headers
    df_list = []
    for file in files:
        # load CSV as DataFrame
        df_data = pd.read_csv(local_dir + file)

        # drop Mean and SD rows
        df_data = df_data.drop(df_data[df_data['Unnamed: 0'] == 'Mean'].index)
        df_data = df_data.drop(df_data[df_data['Unnamed: 0'] == 'SD'].index)

        # get column names from current data
        # cur_column_names = list(df_data.columns)

        # add dataframe to list to concatenate
        df_list.append(df_data)

    df_all = pd.concat(df_list, axis=0)

    return df_all


def check_for_nans(df):
    """
    Report columns in a concatenated DataFrame that were not merged for all samples due to
    mismatched column header names. Mismatched column headers will result in missing row data (NaNs).
    Columns with missing data for all rows will not be reported.

    Parameters
    ----------
    df:
        DataFrame, concatenated from multiple CSV files

    Returns
    -------
    df_nan:
        DataFrame, subset of columns with missing rows (NaNs).

    """
    # Get list of columns with NaN values
    nan_columns = df.columns[df.isna().any()].tolist()
    df_nan = df[nan_columns]

    # Calculate number of NaNs in each column
    nan_counts = df.isna().sum()

    # If NaN counts for any column is ~equal to column length, report warning
    if not (nan_counts == len(df_nan)).any():
        logging.warning('Column names of CSV files are inconsistent. See returned DataFrame for details.')

    # Output DataFrame with NaN columns for inspection
    new_col_names = []
    for col in nan_columns:
        # Make column name the final pop name for easy comparison
        hier, stat = parse_col_name(col)
        new_col_names.append(hier[-1])

    # Concatenate NaN columns with metadata columns
    df_nan = df_nan.rename(columns=dict(zip(nan_columns, new_col_names)))
    df_nan = pd.concat([df[['Unnamed: 0', 'Treatment Batch']], df_nan], axis=1)

    return df_nan


def parse_col_name(col_name):
    """
    Split the full gate path into individual populations and the summary statistic (ie: Counts, Freq of Parent)

    Parameters
    ----------
    col_name:
        The name of the column to parse

    Returns
    -------
    hier:
        List of str, individual gate names returned in path order
    stat:
        String, name of the population summary statistic

    """
    hier = col_name.split('/')
    hier_end = hier[-1].split(' | ')
    hier = hier[:-1] + [hier_end[0]]
    stat = hier_end[-1]

    return hier, stat


def preprocess_csvs(local_path, search_string=None, exclude_files=None):
    """
    Preprocesses data from several CSVs. Loads CSV files from local path, removes NaN columns,
    creates list of metadata columns (mdh_col), and converts several columns to Categorical Type.

    Parameters
    ----------
    local_path : string
       Pathname where CSVs are located
    search_string : string
        Optional. Searches for CSV files with this sub-string in file name.
    exclude_files : list, string
        Optional. List of file names to exclude from loading and concatenating.

    Returns
    -------
    df_data : DataFrame
        Preprocessed, concatenated data from all CSV files
    df_data_nans : DataFrame
        Concatenated data from all CSV files, including only NaN columns (for inspection)
    mdh_col : list
        A list of 'df_data' column names that contain Metadata about each row of 'df_data'

    """
    # Load all data CSVs in local_path and merge into single dataframe
    df_data = load_csvs_to_dataframe(local_path, search_string, exclude_files)

    # Check for NaN columns
    df_data_nans = check_for_nans(df_data)

    # Drop columns that are entirely NaN
    df_data = df_data.dropna(axis=1, how='all')

    # A few other changes
    df_data = df_data.rename(columns={'Count': 'Event Count', 'Unnamed: 0': 'FCS'})

    # Reset row index of DataFrame
    df_data = df_data.sort_values('SampleID', axis=0)
    df_data = df_data.reset_index().drop('index', axis=1)

    # Get list of metadata columns in df
    mdh_col = [col for col in df_data.columns.to_list() if ' | ' not in col]

    # Set mdh_cols to type strings
    df_data[mdh_col[:-1]] = df_data[mdh_col[:-1]].astype('str')

    # Find categorical columns and define the category sort order
    df_data = make_categorical_column(df_data, ['Pilot03', 'Pilot04', 'C1G1', 'C1G2', 'C2G1', 'C2G2', 'C3'])
    df_data = make_categorical_column(df_data, ['Control', 'Ablation', 'Hyperthermia'])
    df_data = make_categorical_column(df_data, ['Contra', 'Ipsi'], ignore_nans=True)

    return df_data, df_data_nans, mdh_col


def make_categorical_column(df, cat_order, ignore_nans=False):
    """
    Searches for a column in 'df' that has a row or rows with a value in 'cat_order.'
    Defines the column as Categorical Type, with sort-order equal to 'cat_order'

    Parameters
    ----------
    df: DataFrame
        DataFrame to search columns of
    cat_order : list
        A list containing all entries of the category, in the specified Categorical order.
    ignore_nans : bool
        Default False. If False, check if any values in categorical column are not matched to a
        value in 'cat_order.' If True, error will not be raised.

    Returns
    -------
    df : DataFrame
        Dataframe 'df' with the converted Categorical Type column

    """
    for col in df.columns:
        for cat in cat_order:
            if cat in df[col].to_list():
                df[col] = pd.Categorical(df[col], categories=cat_order, ordered=True)
                print(col + ' was set to categorical type')

                if not ignore_nans:
                    if df[col].isna().any():
                        raise Exception('Some rows of ' + col + ' were not assigned to Category. Check CSV files for typos.')

                return df

    return df



def create_pop_tree(csv_popnames, df, tissue_type=None, show_tree=True):
    # Load MASTER population names and paired path names into DataFrame
    names_xref = pd.read_csv(csv_popnames, names=['pop_name', 'PATH_NAME_XREF'])

    # Load the path names from the CSV data into DataFrame
    names_data = pd.DataFrame([x for x in df.columns if ' | ' in x], columns=['PATH_NAME_XREF'])
    names_data['PATH_NAME_CLEAN'] = [x.split(' | ')[0] for x in names_data['PATH_NAME_XREF'] if ' | ' in x]

    # Cross-reference the loaded DataFrame column names with the MASTER population names
    try:
        names_data = pd.merge(names_xref, names_data, on='PATH_NAME_XREF', validate='one_to_one')
    except pd.errors.MergeError as e:
        logging.warning('Mismatched population naming key with loaded DataFrame')
        print("MergeError:", e)

    # Create tree from the validated DataFrame path headers
    root = flowtree.PopNode(name="Cells")
    root = bigtree.add_dataframe_to_tree_by_path(root, names_data, path_col='PATH_NAME_CLEAN',
                                                 attribute_cols=['pop_name'])

    # Specify the metadata (MDH) columns for the DataFrame
    mdh_col = [col for col in df.columns.to_list() if ' | ' not in col]

    # Loop through tree descendents, adding data attribute to each tree node
    all_nodes = [root] + [x for x in root.descendants]
    for node in all_nodes:

        # Get full path name from the node
        full_path = node.path_name

        # Search for the full path name in the validated list of populations
        df_col_name = names_data[names_data['PATH_NAME_CLEAN'] == full_path[1:]]
        df_col_name = df_col_name['PATH_NAME_XREF'].to_list()

        # Collect the data and MDH columns from the loaded DataFrame
        data = df[mdh_col + df_col_name]
        data = data.rename(columns={df_col_name[0]: 'Data'})

        # Assign the data to the correct Node attribute
        _, stat = parse_col_name(df_col_name[0])
        if 'count' in stat.lower():

            node.count = data
            node.avg_count = data['Data'].mean()

        elif 'freq' in stat.lower():

            node.freq_of_parent = data
            node.avg_freq_of_parent = data['Data'].mean()

    # Add tissue-type attribute to the tree nodes
    if tissue_type:
        root.tissue = tissue_type

        for n in root.descendants:
            n.tissue = tissue_type

    # Show Tree
    if show_tree:
        root.show(attr_list=['pop_name', 'avg_freq_of_parent'])

    return root


def export_df_to_csv(df, csv_filename, headers=None):

    if headers:
        if not len(df.columns.to_list()) == len(headers):
            logging.warning('Export to CSV error: Header column list does not match DataFrame header columns.')

    with open(csv_filename, 'w') as f:
        # f.write('\n'.join([','.join(line) for line in header]) + '\n')
        if headers:
            writer = csv.writer(f)
            writer.writerow(headers)

        df.to_csv(f, index=False, header=True)


#def check_pop_names(df, csv):
#     pop_gates = [x.split(' | ')[0] for x in df.columns if ' | ' in x]
#     pop_names = []
#     for p in pop_gates:
#         hier, _ = parse_col_name(p)
#         pop_names.append(hier[-1])
#
#     pop_names_lower = [x.lower() for x in pop_names]
#
#     for i, p in enumerate(pop_names_lower):
#         count = pop_names_lower.count(p)
#         if count > 1:
#             hier, _ = parse_col_name(pop_gates[i])
#             print(hier[-2] + '_' + hier[-1])
#             pop_names[i] = hier[-2] + '_' + hier[-1]