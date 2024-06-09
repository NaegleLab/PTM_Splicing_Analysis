import pandas as pd
import numpy as np
import scipy.stats as stats

#utility packages
from tqdm import tqdm
import os

#PTM-POSE
from PTM_POSE import flanking_sequences as fs
from PTM_POSE import database_interfacing as di

def calculateMW_EffectSize(group1, group2):
    """
    Given two lists of values, calculate the effect size and p-value of the Mann-Whitney U test

    Parameters
    ----------
    group1: list or array
        first group of values
    group2: list or array
        second group of values
    
    Returns
    -------
    p: float
        p-value of Mann-Whitney U test
    r: float
        effect size of Mann-Whitney U test
    """
    stat, p = stats.mannwhitneyu(group1, group2)
    n1 = len(group1)
    n2 = len(group2)
    u1 = n1*n2/2
    u2 = n1*n2*(n1+n2+1)/12
    z = (stat - u1)/np.sqrt(u2)
    r = abs(z)/np.sqrt(n1+n2)
    return p, r

def adjust_single_p(p, num_tests, rank, prev_p):
    """
    Given a single p-value, adjust for multiple testing using the Benjamini-Hochberg method

    Parameters
    ----------
    p: float
        p-value to adjust
    num_tests: int
        total number of tests
    rank: int
        rank of the p-value in the sorted list
    prev_p: float
        previous p-value in the list
    """
    adj_p = p*(num_tests/rank)
    if adj_p < prev_p:
        adj_p = prev_p
    elif adj_p > 1:
        adj_p = 1
    else:
        prev_p = adj_p
    return adj_p, prev_p
    

def adjustP(sorted_p, method = 'BH'):
    """
    Given a sorted list of p-values, adjust for multiple testing using the Benjamini-Hochberg method or Bonferroni correction

    Parameters
    ----------
    sorted_p: list
        sorted list of p-values in ascending order
    method: str
        method to use for multiple testing correction. Options include 'BH' (Benjamini-Hochberg) or 'Bonf' (Bonferroni)

    Returns
    -------
    adj_p_list: list
        list of adjusted p-values
    """
    adj_p_list = []
    prev_p = 0
    for i in range(len(sorted_p)):
        if method == 'BH':
            adj_p, prev_p = adjust_single_p(sorted_p[i], len(sorted_p), i+1, prev_p)
            adj_p_list.append(adj_p)
        elif method == 'Bonf':
            adj_p = sorted_p[i]*len(sorted_p)
            if adj_p > 1:
                adj_p = 1
            adj_p_list.append(adj_p)
    return adj_p_list

#create function for calculating parameters for hypergeometric test, then running
def hypergeom(M, n, N, k):
    """
    Given parameters for a hypergeometric test, calculate the p-value of the test

    Parameters
    ----------
    M: int
        total number of genes or PTMs in the dataset
    n: int
        number of genes or PTMs associated with annotation in the dataset
    N: int
        total number of genes or PTMs in the subset
    k: int
        number of genes or PTMs associated with annotation in the subset

    Returns
    -------
    p: float
        p-value of hypergeometric test
    """
    p = stats.hypergeom(M=M, 
                n = n, 
                N=N).sf(k-1)
    return p

def convertToFishers(M, n, N, x):
    """
    Given parameters for a hypergeometric test, convert to a 2x2 table for fishers exact test

        Parameters
    ----------
    M: int
        total number of genes or PTMs in the dataset
    n: int
        number of genes or PTMs associated with annotation in the dataset
    N: int
        total number of genes or PTMs in the subset
    k: int
        number of genes or PTMs associated with annotation in the subset

    Returns
    -------
    table: 2x2 list
        2x2 table for fishers exact test
    """
    table = [[x, n-x],
            [N-x, M- (n+N) + x]]
    return table

def getEnrichment(function_class, all_data, subset_list = None, fishers = True):
    """
    Given a pivot table containing annotations for a set of genes across the entire dataset and a subset of genes in a group of interest, calculate the enrichment of given annotation

    Parameters
    ----------
    function_class: str
        column name of the annotation to test for enrichment
    all_data: pd.DataFrame
        dataframe containing all data to test for enrichment
    subset_list: list
        list of genes or PTMs to test for enrichment
    fishers: bool
        whether to use fishers exact test for enrichment, if False will perform one-tailed hypergeometric test instead

    Returns
    -------
    n: int
        number of genes or PTMs associated with annotation in the dataset
    k: int
        number of genes or PTMs associated with annotation in the subset
    p: float
        p-value of enrichment test
    M: int
        total number of genes or PTMs in the dataset
    N: int
        total number of genes or PTMs in the subset
    odds: float
        odds ratio of enrichment
    """
    M = all_data.shape[0]
    n = all_data.dropna(subset = [function_class]).shape[0]
    mask = all_data.index.isin(subset_list)
    subset_data = all_data.loc[mask]
    N = subset_data.shape[0]
    if not function_class in subset_data.columns:
        k = 0
    else:
        k = subset_data.dropna(subset = [function_class]).shape[0]
    
    if fishers:
        table = convertToFishers(M, n, N, k)
        odds, p = stats.fisher_exact(table)
    else:
        p = hypergeom(M, n, N, k)
    return n, k, p, M, N, odds

def reverse_ESRP1_label_MXE(row):
    """
    Given a row of data from splice seq, determine the ESRP1 label for an MXE event
    """
    if row['Individual exon'] in row['First Exon']:
        return row['ESRP1_MW']
    elif row['ESRP1_MW'] == 'High':
        return 'Low'
    else:
        return 'High'
    
#custom functions
def expand_exon_data(exons, remove_duplicates = False):
    """
    Given exon data from splice seq, separate events with multiple exons into separate groups

    Parameters
    ----------
    exons: pd.DataFrame
        dataframe containing PSI exon data from splice seq
    remove_duplicates: bool
        whether to remove exons that appear across multiple events with conflicting direction (up and downregulated)

    Returns
    -------
    exons: pd.DataFrame
        Updated dataframe with exons separated into individual rows when multiple appearing in the same event
    """
    #split exons separately
    exons['Individual exon'] = exons['exons'].apply(lambda x: x.split(':'))
    exons = exons.explode('Individual exon')

    #if desired remove conflicting exons
    if remove_duplicates:
        exons = exons.drop_duplicates(subset = ['Symbol','Individual exon', 'ESRP1_MW'])
        exons = exons.drop_duplicates(subset = ['Symbol','Individual exon'], keep = False)
    return exons



class ESRP1_analysis:
    """
    Class for analyzing ESRP1 expression and splicing data in TCGA. Will automatically load data from TCGA data directory and extract significant splicing events between ESRP1 high and low groups

    Parameters
    ----------
    tissue: str
        tissue type to analyze, used for extracting data files associated with tissue. Default is 'PRAD' (prostate adenocarcinoma)
    include_ME: bool
        whether to include mutually exclusive exons in the analysis. Default is True

    Attributes
    ----------
    tissue: str
        tissue type analyzed
    figshare_dir: str
        directory containing figshare data
    include_ME: bool
        whether to include mutually exclusive exons in the analysis
    alpha: float
        maximum p-value for significance
    effect_size: float
        minimum effect size (r) for significance
    min_psi_range: float
        minimum required variation in psi values across patients for significance
    cutoff : float
        cutoff for separating ESRP1 high and low groups
    

    Methods
    -------
    get_ESRP1_groups(cutoff = 1)
        Given ESRP1 expression data, separate patients into high and low groups based on a given cutoff
    compare_PSI_for_events_MW(cutoff = 1, min_patients_in_group = 3)
        Given spliceseq PSI data, compare ESRP1-high and ESRP1-low groups using a Mann-Whitney U test
    compare_PSI_for_events_Corr(min_patients = 10)
        Given spliceseq PSI data, compare ESRP1-high and ESRP1-low groups using spearman correlation
    add_exon_annotation()
        Annotate PSI data with additional exon information
    add_ptms_to_data()
        Given annotated spliceseq data, add PTMs to the dataset and create PSI_ptms attribute
    extract_significant_ptms(alpha = 0.05, effect_size = 0.25, min_psi_range = 0.25, by = "MW", duplicate_handling = 'conflicting_sig')
        Given the PSI_ptms data (generate if not yet created), identify the PTMs that are differentially spliced between ESRP1 high and low groups
    """
    def __init__(self, figshare_dir, tissue = "PRAD", include_ME = True, alpha = 0.01, effect_size = 0.3, psi_range = 0.25, cutoff = 1):
        analysis_dir = f'{figshare_dir}/Analysis_For_Paper/'
        database_dir = f'{figshare_dir}/External_Data/'
        
        #load ESRP1 expression data
        mRNA = pd.read_csv(database_dir + f"/TCGA/{tissue}/cBio/mRNA/{tissue}_mRNA expression z-scores relative to all samples (log RNA Seq V2 RSEM).txt", sep = '\t') # Z-score of ESPR1 
        mRNA = mRNA.dropna(subset = 'ESRP1')

        #load processed data if it exists
        self.load_data(tissue = tissue, odir = analysis_dir + '/TCGA/', spliceseq_dir=database_dir + '/TCGA/')
        if not include_ME:
            self.PSI_events = self.PSI_events[self.PSI_events['splice_type'] != 'ME']


        #remove patients with no measured ESRP1 mRNA from splice data, or vice versa
        mRNA['Edited Patient ID'] = mRNA['SAMPLE_ID'].apply(lambda x: '_'.join(x.split('-')[0:3])).to_list()
        patients_to_drop = [col for col in self.PSI_events.columns if col not in mRNA['Edited Patient ID'].values and 'TCGA' in col]
        self.PSI_events = self.PSI_events.drop(patients_to_drop, axis = 1)
        mRNA.set_index('Edited Patient ID', inplace = True)

        #remove patients with no measured splice data from ESRP1 data
        patients_to_drop = [col for col in mRNA.index if col not in self.PSI_events.columns]
        mRNA = mRNA[~mRNA.index.isin(patients_to_drop)].copy()

        self.mRNA = mRNA
        self.mapped_ptms = pd.read_csv(analysis_dir + '/TCGA/splicegraph_ptms.csv', dtype = {'Chromosome': str})
        

        self.splicegraph = pd.read_csv(analysis_dir + "/TCGA/splicegraph_exons.csv") #SpliceSeq Data

        self.get_ESRP1_groups(cutoff = cutoff)

        #store properties
        self.tissue = tissue
        self.odir = analysis_dir + '/TCGA/'
        self.figshare_dir = figshare_dir
        self.analysis_dir = analysis_dir
        self.database_dir = database_dir
        self.alpha = alpha
        self.effect_size = effect_size
        self.min_psi_range = psi_range
        
        #get significant values if exon and ptm info is present
        if self.PSI_exons is not None:
            self.extract_significant_exons() 

        if self.PSI_ptms is not None:
            self.extract_significant_ptms()

        #store properties
        self.tissue = tissue
        self.odir = odir


    def get_ESRP1_groups(self, cutoff = 1):
        """
        Given a specific cutoff, identify patients with high and low expression of ESRP1
        """
        ## Pull low and high 
        ESPR_low = self.mRNA[self.mRNA["ESRP1"] < -1]
        ESRP1_low_id = ESPR_low["SAMPLE_ID"].str.split("-").apply(lambda x: x[0:3]).apply(lambda x: '_'.join(x)).to_list()
        ESPR_high = self.mRNA[self.mRNA["ESRP1"] > 1]
        ESRP1_high_id = ESPR_high["SAMPLE_ID"].str.split("-").apply(lambda x: x[0:3]).apply(lambda x: '_'.join(x)).to_list()
        print(f'Number of ESRP1 low patients: {len(ESRP1_low_id)} ({round(len(ESRP1_low_id)/self.mRNA["SAMPLE_ID"].nunique()*100, 2)}%))')
        print(f'Number of ESRP1 high patients: {len(ESRP1_high_id)} ({round(len(ESRP1_high_id)/self.mRNA["SAMPLE_ID"].nunique()*100, 2)}%)')

        self.ESRP1_low = ESRP1_low_id
        self.ESRP1_high = ESRP1_high_id
        self.ESRP1_cutoff = cutoff

    def compare_PSI_for_events_MW(self, min_patients_in_group = 3):
        """
        Given spliceseq PSI data, compare ESRP1-high and ESRP1-low groups for individual splice events (not exons) using a Mann-Whitney U test

        Parameters
        ----------
        min_patients_in_group: int
            minimum number of patients required in each group to perform statistical test
        """
        if not hasattr(self, 'ESRP1_low'):
            self.get_ESRP1_groups(cutoff = self.ESRP1_cutoff)

        # Hold indexes of PSI data where values are statistically significant + if mean is higher than other group place it in that group. 
        direction = [] 
        p_list = []
        effect_list =[]
        delta_PSI_list = []
        for index, row in tqdm(self.PSI_events.iterrows(), total = self.PSI_events.shape[0], desc = 'Comparing PSI for ESRP1-high and low groups'):
            #grab PSI data for high and low groups
            high_sample = self.PSI_data.loc[index, self.ESRP1_high].values
            high_sample = list(high_sample[~pd.isnull(high_sample)])
            low_sample = self.PSI_data.loc[index, self.ESRP1_low].values
            low_sample = list(low_sample[~pd.isnull(low_sample)])
            #if each sample has at least 3 data points, compare the groups
            if len(low_sample) >= min_patients_in_group and len(high_sample) >= min_patients_in_group:
                p_value, effect_size = calculateMW_EffectSize(high_sample, low_sample)
            else:
                p_value = np.nan
                effect_size = np.nan
            p_list.append(p_value)
            effect_list.append(effect_size)
            
            #if statistical test occurred and p-value was obtained, extract the change in PSI and direction
            delta_PSI = np.mean(high_sample) - np.mean(low_sample)
            if p_value != p_value:
                direction.append(np.nan)
            else:
                delta_PSI = np.mean(high_sample) - np.mean(low_sample)
                if delta_PSI > 0:
                    direction.append('High')
                else: 
                    direction.append('Low')
                delta_PSI_list.append(delta_PSI)
                
        self.PSI_events['ESRP1_MW'] = direction
        self.PSI_events['p_MW'] = p_list
        self.PSI_events['Effect Size_MW'] = effect_list
        self.PSI_events['deltaPSI_MW'] = delta_PSI_list

        #get Benjamini-Hochberg adjusted p-value
        self.PSI_events = self.PSI_data.sort_values(by = 'p_MW', ascending = True)
        self.PSI_events['p-adj_MW'] = adjustP(self.PSI_data['p_MW'].values)

        self.min_patients_in_group = min_patients_in_group


    def compare_PSI_for_exons_MW(self, agg_type = 'mean', min_patients = 3):
        """
        Given spliceseq PSI data, compare average PSI for individual exons across ESRP1-high and ESRP1-low groups using Mann Whitney U test. For exons found across multiple events, use agg_type to determine how to aggregate the data

        Parameters
        ----------
        agg_type: str
            method to use for aggregating exons found in multiple events. Options include 'mean' (average PSI across events) or 'max' (maximum PSI across events)
        min_patients: int
            minimum number of patients required in each group to perform statistical test
        """
        #separate events into individual exons and aggregate exons found in multiple events
        psi_data = expand_exon_data(self.PSI_events, remove_duplicates = False)
        cols = [col for col in psi_data.columns if 'TCGA' in col]
        psi_data = psi_data[['symbol', 'Individual exon']+cols]
        psi_data = psi_data.groupby(['symbol', 'Individual exon'], as_index = False).agg(agg_type)

        direction = [] 
        p_list = []
        effect_list =[]
        delta_PSI_list = []
        for index, row in tqdm(psi_data.iterrows(), total = psi_data.shape[0], desc = 'Comparing PSI for ESRP1-high and low groups'):
            #grab PSI data for high and low groups
            high_sample = psi_data.loc[index, self.ESRP1_high].values
            high_sample = list(high_sample[~pd.isnull(high_sample)])
            low_sample = psi_data.loc[index, self.ESRP1_low].values
            low_sample = list(low_sample[~pd.isnull(low_sample)])
            #if each sample has at least 3 data points, compare the groups
            if len(low_sample) >= min_patients and len(high_sample) >= min_patients:
                p_value, effect_size = calculateMW_EffectSize(high_sample, low_sample)
            else:
                p_value = np.nan
                effect_size = np.nan
            p_list.append(p_value)
            effect_list.append(effect_size)
            
            #if statistical test occurred and p-value was obtained, extract the change in PSI and direction
            delta_PSI = np.mean(high_sample) - np.mean(low_sample)
            if p_value != p_value:
                direction.append(np.nan)
                delta_PSI_list.append(np.nan)
            else:
                delta_PSI = np.mean(high_sample) - np.mean(low_sample)
                if delta_PSI > 0:
                    direction.append('High')
                else: 
                    direction.append('Low')
                delta_PSI_list.append(delta_PSI)
                        
                
        psi_data['ESRP1_MW'] = direction
        psi_data['p_MW'] = p_list
        psi_data['Effect Size_MW'] = effect_list
        psi_data['deltaPSI_MW'] = delta_PSI_list

        #get Benjamini-Hochberg adjusted p-value
        psi_data = psi_data.sort_values(by = 'p_MW', ascending = True)
        psi_data['p-adj_MW'] = adjustP(psi_data['p_MW'].values)

        psi_data = psi_data.reset_index()
        data_cols = [col for col in psi_data if 'TCGA' in col]
        psi_data['psi_range'] = psi_data[data_cols].max(axis = 1) - psi_data[data_cols].min(axis = 1)
        self.PSI_exons = psi_data


    def add_ptms_to_data(self):
        """
        Given annotated spliceseq exon data, add PTMs to the dataset and create PSI_ptms attribute
        """
        PSI_ptms = self.PSI_exons.copy()


        PSI_ptms['Individual exon'] = PSI_ptms['Individual exon'].astype(float)


        #add to PSI data
        PSI_ptms = PSI_ptms.merge(self.mapped_ptms, left_on = ['symbol','Individual exon'], right_on = ['Gene', 'Exon'], how = 'left')
        PSI_ptms = PSI_ptms.dropna(subset = 'PTM')
        PSI_ptms = PSI_ptms.drop(['symbol', 'Exon'], axis = 1)
        self.PSI_ptms = PSI_ptms

    def extract_significant_exons(self, by = "MW", duplicate_handling = 'conflicting_sig'):
        """
        Given the PSI_exons data (generate if not yet created), identify the PTMs that are differentially spliced between ESRP1 high and low groups

        Parameters
        ----------
        by: str
            method to use for significance testing. Either "MW" (mann whitney U test) or "Corr" (spearman rank correlation)
        duplicate_handling: str
            how to handle PTMs involved in multiple splice events. Options include:
                1. 'strict': remove any duplicate entries of PTMs
                2. 'loose': don't remove any duplicates
                3. 'conflicting_sig': remove any entries that appear as significant in both high and low groups
                4. 'conflicting_any': remove any entries that appear with multiple classifications (insignificant, high, or low). This is different from conflicting_sig in that it will also remove any PTM that is significant for one event but not another.
        """
        if not hasattr(self, 'PSI_exons'):
            self.compare_PSI_for_exons_MW()

        
        #identify events that are significant
        sig_exons = self.PSI_exons.copy()
        sig_exons['Significant Event'] = True
        sig_exons.loc[(sig_exons[f'p-adj_{by}'] >= self.alpha), 'Significant Event'] = False
        sig_exons.loc[(sig_exons[f'Effect Size_{by}'] <= self.effect_size), 'Significant Event'] = False
        sig_exons.loc[(sig_exons['psi_range'] <= self.min_psi_range), 'Significant Event'] = False

        #remove duplicate PTM entries according to duplicate_handling
        if duplicate_handling == 'strict': #remove any duplicate entries of PTMs
            sig_exons = sig_exons[sig_exons['Significant Event']]
            sig_exons = sig_exons.drop_duplicates(subset = ['symbol', 'Individual exon'])
        elif duplicate_handling == 'loose': #only remove insignificant events, keep duplicate PTM entries
            sig_exons = sig_exons[sig_exons['Significant Event']]
        elif duplicate_handling == 'conflicting_sig':   #remove events that suggest both up and down regulated
            sig_exons = sig_exons[sig_exons['Significant Event']]
            sig_exons = sig_exons.drop_duplicates(subset = ['symbol','Individual exon', f'ESRP1_{by}'])
            sig_exons = sig_exons.drop_duplicates(subset = ['symbol', 'Individual exon'], keep = False)
        elif duplicate_handling == 'conflicting_any':   #remove events that suggest multiple of up, down, or unregulated
            sig_exons = sig_exons.drop_duplicates(subset = ['symbol','Individual exon', f'ESRP1_{by}', 'Significant Event'])
            sig_exons = sig_exons.drop_duplicates(subset = ['symbol', 'Individual exon'], keep = False)
            sig_exons = sig_exons[sig_exons['Significant Event']]

        #record significant ptms and parameters used to get them
        self.significance_by = by
        self.sig_exons = sig_exons
        self.duplicate_handling = duplicate_handling


    def extract_significant_ptms(self, by = "MW", duplicate_handling = 'conflicting_sig'):
        """
        Given the PSI_ptms data (generate if not yet created), identify the PTMs that are differentially spliced between ESRP1 high and low groups

        Parameters
        ----------
        by: str
            method to use for significance testing. Either "MW" (mann whitney U test) or "Corr" (spearman rank correlation)
        duplicate_handling: str
            how to handle PTMs involved in multiple splice events. Options include:
                1. 'strict': remove any duplicate entries of PTMs
                2. 'loose': don't remove any duplicates
                3. 'conflicting_sig': remove any entries that appear as significant in both high and low groups
                4. 'conflicting_any': remove any entries that appear with multiple classifications (insignificant, high, or low). This is different from conflicting_sig in that it will also remove any PTM that is significant for one event but not another.
        """
        if not hasattr(self, 'PSI_ptms'):
            self.add_ptms_to_data()

        
        #identify events that are significant
        sig_ptms = self.PSI_ptms.copy()
        sig_ptms['Significant Event'] = True
        sig_ptms.loc[(sig_ptms[f'p-adj_{by}'] >= self.alpha), 'Significant Event'] = False
        sig_ptms.loc[(sig_ptms[f'Effect Size_{by}'] <= self.effect_size), 'Significant Event'] = False
        sig_ptms.loc[(sig_ptms['psi_range'] <= self.min_psi_range), 'Significant Event'] = False

        #remove duplicate PTM entries according to duplicate_handling
        if duplicate_handling == 'strict': #remove any duplicate entries of PTMs
            sig_ptms = sig_ptms[sig_ptms['Significant Event']]
            sig_ptms = sig_ptms.drop_duplicates(subset = 'PTM')
        elif duplicate_handling == 'loose': #only remove insignificant events, keep duplicate PTM entries
            sig_ptms = sig_ptms[sig_ptms['Significant Event']]
        elif duplicate_handling == 'conflicting_sig':   #remove events that suggest both up and down regulated
            sig_ptms = sig_ptms[sig_ptms['Significant Event']]
            sig_ptms = sig_ptms.drop_duplicates(subset = ['PTM', f'ESRP1_{by}'])
            sig_ptms = sig_ptms.drop_duplicates(subset = 'PTM', keep = False)
        elif duplicate_handling == 'conflicting_any':   #remove events that suggest multiple of up, down, or unregulated
            sig_ptms = sig_ptms.drop_duplicates(subset = ['PTM', f'ESRP1_{by}', 'Significant Event'])
            sig_ptms = sig_ptms.drop_duplicates(subset = 'PTM', keep = False)
            sig_ptms = sig_ptms[sig_ptms['Significant Event']]

        #record significant ptms and parameters used to get them
        self.significance_by = by
        self.sig_ptms = sig_ptms
        self.duplicate_handling = duplicate_handling

    def get_event_sequences(self, events_to_assess):
        """
        Given a set of splice events, extract the DNA sequences of the flanking and spliced regions for each event

        Parameters
        ----------
        events_to_assess: pd.DataFrame
            dataframe containing splice events to assess

        Returns
        -------
        events_to_assess: pd.DataFrame
            updated dataframe containing sequences for each event
        """
        from_region = []
        from_seq = []
        to_region = []
        to_seq = []
        spliced_seq = []
        for _, event in events_to_assess.iterrows():
            first_exon_region, spliced_regions, second_exon_region = fs.get_spliceseq_event_regions(event, self.splicegraph)
            
            region_list = [first_exon_region] + spliced_regions + [second_exon_region]
            seqs = di.get_region_sequences_from_list(region_list, coordinate_type = 'hg19')
            from_seq.append(seqs[0])
            to_seq.append(seqs[-1])
            spliced_seq.append(''.join(seqs[1:-1]))

            #add region info as well
            from_region.append(','.join([str(first_exon_region[2]), str(first_exon_region[3])]))
            to_region.append(','.join([str(second_exon_region[2]), str(second_exon_region[3])]))


        #add information to dataframe
        events_to_assess['From Region'] = from_region
        events_to_assess['To Region'] = to_region
        events_to_assess['From Sequence'] = from_seq
        events_to_assess['Spliced Sequence'] = spliced_seq
        events_to_assess['To Sequence'] = to_seq
        return events_to_assess

    def get_changed_flanking_sequences(self, ptm_coordinates, flank_size = 5):
        """
        Given a set of significant splice events, identify cases in which PTMs are located near the splice boundary and have potentially altered flanking sequences

        Parameters
        ----------
        ptm_coordinates: pd.DataFrame
            dataframe containing PTM coordinates
        flank_size: int
            size of flanking region to extract for comparison
        
        """
        #load spliceseq
        self.splicegraph['Region ID'] = self.splicegraph['Symbol'] + '_' + self.splicegraph['Exon'].astype(str)
        self.splicegraph.index = self.splicegraph['Region ID'].values

        #get significant events
        if 'p-adj_MW' not in self.PSI_events.columns:
            self.compare_PSI_for_events_MW()
            self.save_data(tissue = self.tissue, odir = self.odir)

        data_for_flanks = self.PSI_events.copy() 
        data_for_flanks = data_for_flanks[(data_for_flanks['p-adj_MW'] <= self.alpha) & (data_for_flanks['Effect Size_MW'] >= self.effect_size)]
        data_for_flanks = data_for_flanks[data_for_flanks['psi_range'] >= self.min_psi_range]
        data_for_flanks = data_for_flanks[['as_id', 'splice_type','symbol', 'from_exon', 'exons', 'to_exon', 'novel_splice', 'p-adj_MW', 'Effect Size_MW', 'deltaPSI_MW']].drop_duplicates()

        #load splicegraph ptms and add strand info
        self.mapped_ptms = self.mapped_ptms.merge(self.splicegraph[['Region ID', 'Strand']], how = 'left', on = 'Region ID')

        #grab events with PTMs at splice boundary
        first_flank = data_for_flanks.copy() #first flank information (from_exon)
        first_flank['Region ID'] = first_flank['symbol']+'_'+first_flank['from_exon'].astype(str)
        regions_with_ptms= self.mapped_ptms.loc[(self.mapped_ptms['Proximity to Region End (bp)'] < flank_size * 3), 'Region ID'].unique() #extract ptms that are nearby the splice boundary
        first_flank = first_flank[first_flank['Region ID'].isin(regions_with_ptms)] #restrict to events with PTMs at splice boundary

        second_flank = data_for_flanks.copy() #second flank information (to_exon)
        second_flank['Region ID'] = second_flank['symbol']+'_'+second_flank['to_exon'].astype(str)
        regions_with_ptms = self.mapped_ptms.loc[(self.mapped_ptms['Proximity to Region Start (bp)'] < flank_size * 3), 'Region ID'].unique() #extract ptms that are nearby the splice boundary
        second_flank = second_flank[second_flank['Region ID'].isin(regions_with_ptms)]


        #get sequences for each event, or load if already generated
        if os.path.exists(self.odir + f'event_sequences.csv'):
            events_to_assess = pd.read_csv(self.odir + f'event_sequences.csv')
        else:
            #combine information from both flanks
            events_to_assess = pd.concat([first_flank, second_flank]).drop_duplicates()
            events_to_assess = self.get_event_sequences(events_to_assess)
            #save to file
            events_to_assess.to_csv(self.odir+f'event_sequences.csv', index = False)

        trim_splicegraph = self.mapped_ptms[self.mapped_ptms['Region ID'].isin(events_to_assess['Region ID'])].copy()
        #add gene location to splicegraph
        trim_splicegraph = trim_splicegraph.merge(ptm_coordinates[['Source of PTM', 'Gene Location (hg19)', 'Expected Flanking Sequence']], on = 'Source of PTM', how = 'left')
                
        #add ptms to flanking information and combine into final dataframe
        from_ptms = first_flank.merge(trim_splicegraph.loc[trim_splicegraph['Proximity to Region End (bp)'] < flank_size * 3, ['Region ID', 'Gene', 'Strand','Source of PTM', 'PTM', 'Residue','Modification Class', 'Gene Location (hg19)', 'Expected Flanking Sequence']], on = 'Region ID', how = 'inner')
        from_ptms['Which Flank'] = 'First'
        to_ptms = second_flank.merge(trim_splicegraph.loc[trim_splicegraph['Proximity to Region Start (bp)'] < flank_size * 3, ['Region ID', 'Gene','Strand', 'Source of PTM', 'PTM', 'Residue', 'Modification Class', 'Gene Location (hg19)', 'Expected Flanking Sequence']], on = 'Region ID', how = 'inner')
        to_ptms['Which Flank'] = 'Second'
        flank_ptms = pd.concat([from_ptms, to_ptms]).copy()
        flank_ptms = flank_ptms.reset_index(drop = True)

        def get_flank_loc(ptm):
            if ptm['Strand'] == '+' and ptm['Which Flank'] == 'First':
                return ptm['Gene Location (hg19)']  - int(event['From Region'].split(',')[0])
            elif ptm['Strand'] == '+' and ptm['Which Flank'] == 'Second':
                return ptm['Gene Location (hg19)']  - int(event['To Region'].split(',')[0])
            elif ptm['Strand'] == '-' and ptm['Which Flank'] == 'First':
                return int(event['From Region'].split(',')[1]) - ptm['Gene Location (hg19)']
            else:
                return int(event['To Region'].split(',')[1]) - ptm['Gene Location (hg19)']

        for _,event in events_to_assess.iterrows():    
            #grab inclusion and exclusion event sequences (exclusion does not contain spliced sequences)
            inclusion_seq = event['From Sequence'] + event['Spliced Sequence'] + event['To Sequence']
            exclusion_seq = event['From Sequence'] + event['To Sequence']

            #grab ptms associated with splice event
            ptms_in_region = flank_ptms[flank_ptms['as_id'] == event['as_id']].copy()

            #calculate the location of the ptm in the relevant flank (one where ptm is located)
            ptms_in_region['PTM Loc in Flank'] = ptms_in_region.apply(get_flank_loc, axis = 1)

            for i, ptm in ptms_in_region.iterrows():
                #grab where ptm is located in both the inclusion and exclusion event
                inclusion_ptm_loc, exclusion_ptm_loc = fs.get_ptm_locs_in_spliced_sequences(ptm['PTM Loc in Flank'], event['From Sequence'], event['Spliced Sequence'], event['To Sequence'], which_flank = ptm['Which Flank'], order_by = 'Translation')

                #extract expected flanking sequence based on location in sequence
                inclusion_flank = fs.get_flanking_sequence(inclusion_ptm_loc, inclusion_seq, ptm_residue = ptm['Residue'], flank_size = flank_size, full_flanking_seq = False)
                exclusion_flank = fs.get_flanking_sequence(exclusion_ptm_loc, exclusion_seq, ptm_residue = ptm['Residue'], flank_size = flank_size, full_flanking_seq = False)

                #add to dataframe
                flank_ptms.loc[i, 'Inclusion Flanking Sequence'] = inclusion_flank
                flank_ptms.loc[i, 'Exclusion Flanking Sequence'] = exclusion_flank
        
        flank_ptms['Expected Flanking Sequence'] = flank_ptms['Expected Flanking Sequence'].apply(lambda x: x[int((len(x)-1)/2-flank_size):int((len(x)-1)/2+flank_size+1)] if x == x else np.nan)
        flank_ptms['Matched'] = flank_ptms['Inclusion Flanking Sequence'] == flank_ptms['Exclusion Flanking Sequence']
        self.flanking_sequences = flank_ptms.copy()

        #save to file
        flank_ptms.to_csv(self.odir + f'changed_flank_sequences_{self.tissue}.csv', index = False)

    def load_data(self, odir, tissue = 'PRAD', spliceseq_dir = None):
        #load processed data if it exists
        if os.path.exists(odir + f"PSI_events_{tissue}.csv"):
            self.PSI_events = pd.read_csv(odir + f"PSI_events_{tissue}.csv")
        else:
            if spliceseq_dir is not None:
                self.PSI_events = pd.read_csv(spliceseq_dir + f"PSI_download_{tissue}.txt", sep = '\t')
            else:
                raise ValueError('No event data from spliceseq found, please add to folder or indicate where PSI_download data is with spliceseq_dir parameter.')

        #load processed data if it exists
        if os.path.exists(odir + f'PSI_ptms_{tissue}.csv'):
            self.PSI_ptms = pd.read_csv(odir + f'PSI_ptms_{tissue}.csv')
        else:
            self.PSI_ptms = None

        if os.path.exists(odir + f'PSI_exons_{tissue}.csv'):
            self.PSI_exons = pd.read_csv(odir + f'PSI_exons_{tissue}.csv')
        else:
            self.PSI_exons = None

        if os.path.exists(odir + f'{tissue}_annotated_ptms.csv'):
            self.annotated_ptms = pd.read_csv(odir + f'{tissue}_annotated_ptms.csv')
        else:
            self.annotated_ptms = None

        if os.path.exists(odir + f'changed_flank_sequences_{tissue}.csv'):
            self.flanking_sequences = pd.read_csv(odir + f'changed_flank_sequences_{tissue}.csv')
        else:
            self.flanking_sequences = None
        

    def save_data(self, tissue, odir):
        #load processed data if it exists
        if self.PSI_events is not None:
            self.PSI_events.to_csv(f'{odir}/PSI_events_{tissue}.csv', index = False)
        
        if self.PSI_ptms is not None:
            self.PSI_ptms.to_csv(f'{odir}/PSI_ptms_{tissue}.csv', index = False)
        
        if self.PSI_exons is not None:
            self.PSI_exons.to_csv(f'{odir}/PSI_exons_{tissue}.csv', index = False)

        if self.annotated_ptms is not None:
            self.annotated_ptms.to_csv(f'{odir}/{tissue}_annotated_ptms.csv', index = False)
        
        if self.flanking_sequences is not None:
            self.flanking_sequences.to_csv(f'{odir}/changed_flank_sequences_{tissue}.csv', index = False)
        
        


    




