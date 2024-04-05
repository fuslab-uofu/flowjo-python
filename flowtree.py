import bigtree
import logging
import pandas as pd
import csv

# def my_decorator(func):
#     def wrapper(*args, **kwargs):
#         print("Before calling", func.__name__)
#         result = func(*args, **kwargs)
#         print("After calling", func.__name__)
#         return result
#     return wrapper
#
# # Apply the decorator to every function in the module
#
# this_module = bigtree.__name__
# for name, obj in vars(this_module).items():
#     if callable(obj):
#         setattr(this_module, name, my_decorator(obj))


def filter_node(node, var_column, keep_values=None, drop_values=None):
    """
    Filter rows of all DataFrame attributes of the node by values of a specific column.
    Typically used to filter by metadata columns (ie: Treatment Batch, etc).

    Parameters
    ----------
    node : object
        PopNode object, node to filter
    var_column : str
        Column to filter data on
    keep_values : list
        Keep rows where var_column is in list of keep_values
    drop_values : list
        Drop rows where var_column is in list of drop_values
    """

    df_fields = {k: v for k, v in vars(node).items() if isinstance(v, pd.DataFrame)}

    for attr_name, df in df_fields.items():
        df_temp = df.copy()

        if keep_values:
            df_temp = df_temp.loc[df_temp[var_column].isin(keep_values)]

        if drop_values:
            df_temp = df_temp.loc[~df_temp[var_column].isin(keep_values)]

        setattr(node, attr_name, df_temp)
        del df, df_temp


class PopNode(bigtree.Node):

    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)
        self.name = name
        self.pop_name = None

    def get_list_pop_names(self):
        """Returns list of 'pop_name' attribute for all nodes"""

        pop_names = [self.root.pop_name]
        for node in self.root.descendants:
            pop_names.append(node.pop_name)
            print(node.pop_name)

        return pop_names

    def show_pop_name_tree(self):
        """Prints flowtree structure with 'pop_name' attribute for all nodes"""

        new_self = self.copy()
        new_self.name = new_self.pop_name
        for node in new_self.descendants:
            node.name = node.pop_name

        new_self.show()

    def find_popname(self, value):
        """Finds node in flowtree by 'pop_name' attribute."""

        try:
            node = bigtree.find_attr(self.root, attr_name='pop_name', attr_value=value)
            return node
        except Exception as err:
            print("The error is: ", err)
            raise RuntimeError('Population name was not found in PopTree.')

    def calculate_counts_tree(self):
        """Calculates event Counts on all descendants of self.
        Creates the 'counts' attribute for each descendent node of self."""

        for node in self.descendants:
            node.calculate_counts()

    def calculate_counts(self):
        """Calculates event Counts on self by propagating freq_of_parent from flowtree 'Cells' node.
        Creates the 'counts' attribute for self."""

        # print('Calculating counts for ' + self.pop_name + ':')
        backgate_nodes = self.go_to(bigtree.find_name(self.root, 'Cells'))

        # Starting at the highest ancestor of the population of interest:
        for node in reversed(backgate_nodes):

            # If counts have not been calculated
            if not hasattr(node, 'counts'):

                # If there is no parent - use flowtree root
                if (node.parent is None) or node.is_root:

                    # Calculate counts from the Event Count column of mdh
                    # print('Calculating counts for backgate ' + node.path_name)
                    counts = (node.freq_of_parent['Data'] * node.freq_of_parent['Event Count'].astype("float64")) / 100
                    data = node.freq_of_parent.copy().drop('Data', axis=1)
                    data['Data'] = counts
                    node.counts = data

                    # Check dataframe is the same size as input
                    if not (node.freq_of_parent.shape[0] == node.counts.shape[0]):
                        logging.warning('Number of samples for Freq of Parent and Counts are not consistent.')

                # Else if there is a parent
                else:

                    # Get counts from parent
                    # print('Calculating counts for backgate ' + node.path_name)
                    merged_data = node.freq_of_parent.merge(node.parent.counts, on='FCS',
                                                            suffixes=('_freq', '_count'))
                    counts = (merged_data['Data_freq'] * merged_data['Data_count']) / 100
                    data = node.freq_of_parent.copy().drop('Data', axis=1)
                    data['Data'] = counts
                    node.counts = data

                    # Check dataframe is the same size as input
                    if not (node.freq_of_parent.shape[0] == node.counts.shape[0]):
                        logging.warning('Number of samples for Freq of Parent and Counts are not consistent.')

    def get_freq_of_ancestor(self, ancestor, child=None):
        """Calculates population frequency of a specified ancestor.
        Specify alternate population node with optional 'child' input.

        Parameters
        ----------
        ancestor : str
            'node_path' or 'popname' attribute of ancestor population node
        child : str
            'node_path' or 'popname' attribute of child population node

        Returns
        ----------
        freq : object
            DataFrame with calculated frequency stat
        stat_name : str
            Stat description using node 'name' attributes
        stat_name_pop : str
            Stat description using node 'popname' attributes
        """

        # search for child node, if not self (if input by user)
        if child:
            if '/' in child:
                child_node = bigtree.find_full_path(self.root, child)

            else:
                child_node = self.find_popname(child)

        else:
            child_node = self

        # search for ancestor node based on the input string (type is full_path or popname)
        if '/' in ancestor:

            ancestor_node = child_node.go_to(bigtree.find_full_path(child_node.root, ancestor))

            if self.pop_name == ancestor:
                ancestor_node = ancestor_node[0]
            else:
                ancestor_node = ancestor_node[-1]

        else:

            ancestor_node = child_node.go_to(child_node.find_popname(ancestor))

            if self.pop_name == ancestor:
                ancestor_node = ancestor_node[0]
            else:
                ancestor_node = ancestor_node[-1]

        # If counts have not been calculated for child node, calculate counts
        if not hasattr(child_node, 'counts'):
            # raise Exception('Node missing "Count" attribute: ' + node.name)
            child_node.calculate_counts()

        # Initialize dataframe to hold counts of ancestor populations
        pop_counts = ancestor_node.counts
        merge_on = pop_counts.drop(columns='Data').columns.to_list()

        # Merge Counts data from both child and ancestor nodes
        pop_counts = pop_counts.merge(child_node.counts, on=merge_on, suffixes=['_Anc', '_Des'], validate='one_to_one')

        # Calculate the child population frequency of ancestor
        pop_counts['Data'] = 100*(pop_counts['Data_Des']/pop_counts['Data_Anc'])

        # Drop the count columns
        freq = pop_counts.drop(columns=['Data_Des', 'Data_Anc'])

        # Export the new stat name
        stat_name = child_node.path_name + ' | Freq of ' + ancestor_node.path_name + ' (%)'
        stat_name_pop = child_node.pop_name + ' | Freq of ' + ancestor_node.pop_name + ' (%)'

        return freq, stat_name, stat_name_pop

    def get_freq(self):
        """Returns node 'freq_of_parent' attribute, or creates it if it does not exist.

        Returns
        ----------
        freq : object
            DataFrame with calculated freq_of_parent stat
        stat_name : str
            Stat description using node 'name' attributes
        stat_name_pop : str
            Stat description using node 'popname' attributes
        """

        if hasattr(self, 'freq_of_parent'):
            freq = self.freq_of_parent
            stat_name = self.name + ' | Freq of ' + self.parent.name + ' (%)'
            stat_name_pop = self.pop_name + ' | Freq of ' + self.parent.pop_name + ' (%)'

        else:
            freq, stat_name, stat_name_pop = self.get_freq_of_ancestor(self, self.parent.pop_name)

        return freq, stat_name, stat_name_pop

    def get_count(self):
        """Returns node 'counts' attribute, or creates it if it does not exist.

        Returns
        ----------
        count : object
            DataFrame with calculated Counts stat
        stat_name : str
            Stat description using node 'name' attributes
        stat_name_pop : str
            Stat description using node 'popname' attributes
        """

        if hasattr(self, 'counts'):
            count = self.counts
        else:
            self.calculate_counts()
            count = self.counts

        stat_name = self.name + ' | Count'
        stat_name_pop = self.pop_name + ' | Count'

        return count, stat_name, stat_name_pop

    def append_mrti_data(self, df_mrti, on_column="SampleID"):
        """Merges MRTI dataframe with flowtree root dataframes on 'SampleID' column.
        Creates flowtree root node 'mrti' attribute.

        Parameters
        ----------
        df_mrti : object
            DataFrame containing MRTI statistics
        on_column : str
            Column to perform DataFrame outer merge on. Default is "SampleID" to merge on
            the individual samples in the flowtree DataFrames

        Returns
        ----------
        df_root_new : object
            DataFrame of MRTI statistics merged with flowtree DataFrames
        """

        df_root = self.root.freq_of_parent.copy()
        df_root.drop('Data', axis=1, inplace=True)
        df_root_new = pd.merge(df_root, df_mrti, on=on_column, how='outer')
        df_root_new = df_root_new[~df_root_new['FCS'].isna()]

        # check that df_root is the same size as before merge
        if not len(df_root) == len(df_root_new):
            logging.warning('MRTI Merge error: number of FCS files in dataset has changed.')

        self.root.mrti = df_root_new

        return df_root_new  # was df_root. changed 4/5/2024. Check if correct.

    def export_tree_as_dataframe(self, to_csv=True, csv_filename='tree_data', merge_mrti=True):
        """
        Export flowtree 'counts' and 'freq_of_parent' attributes to DataFrame and/or CSV.
        Option to include 'mrti' attribute as well.

        Parameters
        ----------
        to_csv : bool
            if True, export DataFrame to csv.
        csv_filename : str
            Pre-fix for exported CSV filename
        merge_mrti : bool
            if True, include 'mrti' data from the flowtree in exported DataFrame/CSV file

        Returns
        ----------
        df_counts : object
            DataFrame with merged 'counts' attribute for all nodes in the flowtree
        df_freq : object
            DataFrame with merged 'freq_of_parent' attribute for all nodes in the flowtree.
        """

        # Export Freq of Parent
        df_out = self.freq_of_parent.copy()
        merge_on = df_out.drop(columns='Data').columns.to_list()

        if merge_mrti:
            df_out = pd.merge(self.root.mrti, df_out, on=merge_on, validate='one_to_one')

        # root node exported data
        df_out.rename(columns={'Data': self.pop_name + ' | Freq of Parent (%)'}, inplace=True)

        # root descendant nodes exported data
        header_fullpath = df_out.columns.to_list()
        for node in self.descendants:
            df = node.freq_of_parent.copy()
            header_fullpath.append(node.path_name + ' | Freq of Parent (%)')
            df.rename(columns={'Data': node.pop_name + ' | Freq of Parent (%)'}, inplace=True)
            df_out = pd.merge(df_out, df, on=merge_on, validate='one_to_one')

        # export to csv
        if to_csv:
            if not len(df_out.columns.to_list()) == len(header_fullpath):
                logging.warning('Export to CSV error: Header column list does not match DataFrame header columns.')

            with open(csv_filename + '_Freq_of_Parent.csv', 'w') as f:
                # f.write('\n'.join([','.join(line) for line in header]) + '\n')
                writer = csv.writer(f)
                writer.writerow(header_fullpath)
                df_out.to_csv(f, index=False, header=True)

        df_freq = df_out.copy()
        del df_out

        # Export Counts
        df_out = self.counts.copy()
        merge_on = df_out.drop(columns='Data').columns.to_list()

        if merge_mrti:
            df_out = pd.merge(self.root.mrti, df_out, on=merge_on, validate='one_to_one')

        # root node exported data
        df_out.rename(columns={'Data': self.pop_name + ' | Freq of Parent (%)'}, inplace=True)

        # root descendant nodes exported data
        header_fullpath = df_out.columns.to_list()
        for node in self.descendants:
            df = node.counts.copy()
            header_fullpath.append(node.path_name + ' | Count')
            df.rename(columns={'Data': node.pop_name + ' | Count'}, inplace=True)
            df_out = pd.merge(df_out, df, on=merge_on, validate='one_to_one')

        # export to csv
        if to_csv:
            if not len(df_out.columns.to_list()) == len(header_fullpath):
                logging.warning('Export to CSV error: Header column list does not match DataFrame header columns.')

            with open(csv_filename + '_Counts.csv', 'w') as f:
                # f.write('\n'.join([','.join(line) for line in header]) + '\n')
                writer = csv.writer(f)
                writer.writerow(header_fullpath)
                df_out.to_csv(f, index=False, header=True)

        df_counts = df_out.copy()

        return df_counts, df_freq

    def export_freqs_as_dataframe(self, sub_populations, populations,
                                  to_csv=True, csv_filename='tree_data',
                                  merge_mrti=True, freq_of_parent=True):
        """
        Export custom population frequency data from the flowtree to CSV.

        Parameters
        ----------
        sub_populations : list[str]
            a list of strings, which are the population alias names (pop_name) for each "child" population to calculate
            the frequencies of
        populations : list[str]
            a list of stings, which are the population alias names (pop_name) for each ancestor population that each
            child frequency will be calculated from
        to_csv : bool
            if True, export DataFrame to csv.
        csv_filename : str
            Pre-fix for exported CSV filename
        merge_mrti : bool
            if True, include 'mrti' data from the flowtree in exported DataFrame/CSV file
        freq_of_parent : bool
            if True, include the 'freq_of_parent' attribute for every node in 'sub_populations'

        Returns
        ----------
        df_out : object
            DataFrame of all merged data
        header_fullpath : list
            Column names in df_out, using the full pathnames of each frequency statistic.
        """

        df_out = self.counts.copy()
        df_out.drop(columns='Data', inplace=True)
        merge_on = df_out.columns.to_list()  # drop(columns='Data').columns.to_list()

        if merge_mrti:
            df_out = pd.merge(self.root.mrti, df_out, on=merge_on, validate='one_to_one')

        header_fullpath = df_out.columns.to_list()

        # if freq_of_parent:
        #     for sub_pop in sub_populations:
        #
        #         sub_pop_node = bigtree.find_attr(self.root, 'pop_name', sub_pop)
        #         df = sub_pop_node.freq_of_parent.copy()
        #         header_fullpath.append(sub_pop_node.path_name + ' | Freq of Parent (%)')
        #         df.rename(columns={'Data': sub_pop_node.pop_name + ' | Freq of Parent (%)'}, inplace=True)
        #         df_out = pd.merge(df_out, df, on=merge_on, validate='one_to_one')

        for sub_pop in sub_populations:

            # pop_node = bigtree.find_attr(self.root, 'pop_name', pop)
            sub_pop_node = bigtree.find_attr(self.root, 'pop_name', sub_pop)

            if freq_of_parent:

                sub_pop_node = bigtree.find_attr(self.root, 'pop_name', sub_pop)
                df = sub_pop_node.freq_of_parent.copy()
                header_fullpath.append(sub_pop_node.path_name + ' | Freq of Parent (%)')
                df.rename(columns={'Data': sub_pop_node.pop_name + ' | Freq of Parent (%)'}, inplace=True)
                df_out = pd.merge(df_out, df, on=merge_on, validate='one_to_one')

            for pop in populations:

                # check that ancestor is an ancestor of child_node:
                if pop not in [x.pop_name for x in sub_pop_node.ancestors]:
                    logging.warning(pop + ' is not an ancestor of ' + sub_pop)
                    continue

                df, stat_name, stat_name_pop = sub_pop_node.get_freq_of_ancestor(pop, child=None)
                header_fullpath.append(stat_name)
                df.rename(columns={'Data': stat_name_pop}, inplace=True)
                df_out = pd.merge(df_out, df, on=merge_on, validate='one_to_one')

        # export to csv
        if to_csv:
            if not len(df_out.columns.to_list()) == len(header_fullpath):
                logging.warning('Export to CSV error: Header column list does not match DataFrame header columns.')

            with open(csv_filename, 'w') as f:
                # f.write('\n'.join([','.join(line) for line in header]) + '\n')
                writer = csv.writer(f)
                writer.writerow(header_fullpath)
                df_out.to_csv(f, index=False, header=True)

        return df_out, header_fullpath

    def filter_by_cat(self, var_column, keep_values=None, drop_values=None):
        """
        For all nodes in flowtree, filter rows of 'counts' and 'freq_of_parent' by values of var_column.
        Typically used to filter by metadata columns (ie: Treatment Batch, etc)

        Parameters
        ----------
        var_column : str
            Column to filter data on
        keep_values : list
            Keep rows where var_column is in list of keep_values
        drop_values : list
            Drop rows where var_column is in list of drop_values
        """

        filter_node(self.root, var_column, keep_values, drop_values)

        for node in self.root.descendants:
            # df_fields = {k: v for k, v in vars(node).items() if isinstance(v, pd.DataFrame)}

            filter_node(node, var_column, keep_values, drop_values)

        return self

    # if '/' in ancestor:
    #     backgate_nodes = self.go_to(bigtree.find_full_path(self.root, ancestor))
    #
    # else:
    #     backgate_nodes = self.go_to(bigtree.find_attr(self.root, attr_name='pop_name', attr_value=ancestor))
    #
    # # Initialize dataframe to hold counts of backgated populations
    # df_backgate_counts = self.freq_of_parent.drop(columns='Data')
    # merge_on = df_backgate_counts.columns.to_list()
    #
    # # Starting at the highest ancestor of the population of interest:
    # for i, node in enumerate(reversed(backgate_nodes)):
    #
    #     # If counts have not been calculated
    #     if not hasattr(node, 'counts'):
    #         #raise Exception('Node missing "Count" attribute: ' + node.name)
    #         node.calculate_counts()
    #
    #     df_backgate_counts = df_backgate_counts.merge(node.counts,on=merge_on,validate='one_to_one')
    #     df_backgate_counts.rename(columns={'Data': i}, inplace=True)
    #
    # for i in range(1, len(backgate_nodes)):
    #
    #     freq = df_backgate_counts[i]/df_backgate_counts[i-1]
    #
    #
    #
