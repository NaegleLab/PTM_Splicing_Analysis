import pandas as pd
from ExonPTMapper import config, mapping, get_splice_events
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np
import seaborn as sns


def getOverallDensities_Canonical(mapper, exploded_mod):
	"""
	Compute overall PTM density in the proteome, as well as the density broken down by modification type
	
	Parameters
	----------
	mapper: ExonPTMapper object
    
	exploded_mod: Pandas dataframe
        Version of mapper.ptm_info with unique modification types separated into different rows
	"""
	ptm_info = mapper.ptm_info[mapper.ptm_info['Isoform Type'] == 'Canonical'].copy()
	exploded_mod = exploded_mod[exploded_mod['Isoform Type'] == 'Canonical'].copy()
	#grab first transcript associated with canonical protein, and use this information (avoid double counting the same protein)
	ptm_info['First Transcript'] = ptm_info['Transcripts'].apply(lambda x: x.split(';')[0])
	ptm_number = ptm_info.groupby('First Transcript').size()

	#get protein length for	 transcripts
	tmp_transcripts = mapper.transcripts.dropna(subset = 'Amino Acid Sequence').copy()
	tmp_transcripts['Protein Length'] = tmp_transcripts['Amino Acid Sequence'].apply(len)

	#get density of ptms relative to total protein length
	overall_ptm_density = (ptm_number.sum()/tmp_transcripts.loc[ptm_number.index, 'Protein Length'].sum())*1000
	mod_densities = (exploded_mod.groupby('Modification Class').size()/tmp_transcripts.loc[ptm_number.index, 'Protein Length'].sum())*1000


	#get overall densities
	overall_density = mod_densities.reset_index()
	overall_density = overall_density.rename({0:'Density in Proteome (per 1000aa)'}, axis = 1)
	overall_proteome = pd.DataFrame({'Modification Class': 'All', 'Density in Proteome (per 1000aa)':overall_ptm_density}, index = [300])
	overall_density = pd.concat([overall_proteome,overall_density])
	overall_density.index = overall_density['Modification Class']
	overall_density = overall_density.drop(['Modification Class'], axis = 1)
	return overall_density

def getOverallDensities_Isoforms(mapper):
	iso_ptms = mapper.isoform_ptms[mapper.isoform_ptms['Mapping Result'] == 'Success'].copy()
	num_ptms = iso_ptms.groupby('Isoform ID').size()
	num_residues = iso_ptms.groupby('Isoform ID')['Isoform Length'].mean()
	mod_classes = ['All']
	ptm_density = [(num_ptms/num_residues.loc[num_ptms.index]).mean()*1000]

	exploded_iso = iso_ptms.copy()
	exploded_iso['Modification Class'] = exploded_iso['Modification Class'].apply(lambda x: x.split(';') if x == x else np.nan)
	exploded_iso = exploded_iso.explode('Modification Class').drop_duplicates()
	for mod in exploded_iso['Modification Class'].unique():
		mod_classes.append(mod)
		num_ptms = exploded_iso[exploded_iso['Modification Class'] == mod].groupby('Isoform ID').size()
		num_residues = exploded_iso[exploded_iso['Modification Class'] ==  mod].groupby('Isoform ID')['Isoform Length'].mean()
		ptm_density.append((num_ptms/num_residues.loc[num_ptms.index]).mean()*1000)

	overall_density = pd.DataFrame({'Density in Proteome (per 1000aa)':ptm_density}, index = mod_classes)
	return overall_density

def addPTMsToExons(exons, mapper, exploded_ptms, mod_groups, exon_column = 'Exon'):
	"""
	Given dataframe of tissue specific exons, add PTM information to each exon
	"""
	#check canonical exons for ptm info first
	#ts_ptms = exons.merge(exploded_ptms, left_on = exon_column, right_on = 'Exon stable ID', how = 'left')
	#ts_ptms = ts_ptms[[exon_column, 'PTM', 'Modification', 'Exon AA Seq (Full Codon)']].drop_duplicates()
	#ts_ptms = ts_ptms.rename({'PTM':'Source of PTM'}, axis = 1)
	#check alternative exons for ptm info next
	ts_ptms2 = exons.merge(mapper.alternative_ptms[mapper.alternative_ptms['Mapping Result'] == 'Success'], left_on = exon_column, right_on = 'Exon ID (Alternative)', how = 'left')
	ts_ptms2 = ts_ptms2[[exon_column, 'Source of PTM', 'Modification Class', 'Exon AA Seq (Full Codon)']].drop_duplicates()
	#combine information
	#ts_ptms = pd.concat([ts_ptms, ts_ptms2])
	ts_ptms = ts_ptms2
	ts_ptms = ts_ptms.rename({exon_column:'Exon stable ID'}, axis = 1)
	#separate unique Modification
	ts_ptms['Modification Class'] = ts_ptms['Modification Class'].apply(lambda x: x.split(';') if x == x else np.nan)
	ts_ptms = ts_ptms.explode('Modification Class')
	#ts_ptms = ts_ptms.merge(mod_groups, on = 'Modification', how = 'left')
	#ts_ptms = ts_ptms.drop(['Modification'], axis = 1)
	ts_ptms = ts_ptms.drop_duplicates()
	return ts_ptms
	
def getBuljanPTMs(trim_exons, mapper, exploded_ptms, mod_groups, figshare_dir = './PTM_Splicing_FigShare/'):
	"""
	Extract which exons were identified as tissue specific by Buljan et al. (2012), and which PTMs are found in these exons.
	"""
	#load tissue-specific data
	tse_exons = pd.read_csv(figshare_dir + 'External_Data/TissueSpecificExons/Buljan2012.txt', sep = '\t')

	#merge tissue-specific information with our exon data
	tse_exons = trim_exons.merge(tse_exons, right_on = ['ensembl_gene', 'ensembl_transcript', 'exon'], 
								   left_on = ['Gene stable ID', 'Transcript stable ID', 'Exon rank in transcript'], how = 'inner')

	#identify ptms present in canonical tse_exons
	tse_ptms = addPTMsToExons(tse_exons, mapper, exploded_ptms, mod_groups, exon_column = 'Exon stable ID')
	return tse_ptms

def getGP_PTMs(mapper, trim_exons, exploded_ptms, mod_groups, figshare_dir = '/TissueSpecificData/gonzalez-porta_5fold.txt'):
	#load data
	gp_transcripts = pd.read_csv(figshare_dir + '/External_Data/TissueSpecificExons/Gonzalez-porta2013.txt', sep = ' ')
	#gp_transcripts['Sample A Tissue'] = gp_transcripts['sample_A:major_transcript_sample_A:fpkm'].apply(lambda x: x.split(':')[0])
	gp_transcripts['Sample A Transcript'] = gp_transcripts['sample_A:major_transcript_sample_A:fpkm'].apply(lambda x: x.split(':')[1])
	#gp_transcripts['Sample A FPKM'] = gp_transcripts['sample_A:major_transcript_sample_A:fpkm'].apply(lambda x: x.split(':')[2])
	#gp_transcripts['Sample B Tissue'] = gp_transcripts['sample_B:major_transcript_sample_B:fpkm'].apply(lambda x: x.split(':')[0])
	gp_transcripts['Sample B Transcript'] = gp_transcripts['sample_B:major_transcript_sample_B:fpkm'].apply(lambda x: x.split(':')[1])
	#gp_transcripts['Sample B FPKM'] = gp_transcripts['sample_B:major_transcript_sample_B:fpkm'].apply(lambda x: x.split(':')[2])

	#identify tissue-specific exons (skipped exons in one of the tissue-specific transcripts)
	sevents = []
	for e in gp_transcripts.index:
		tmp_data = gp_transcripts.loc[e]
		main = tmp_data['Sample A Transcript']
		alt = tmp_data['Sample B Transcript']
		main_exons = trim_exons[trim_exons['Transcript stable ID'] == main]
		alternative_exons = trim_exons[trim_exons['Transcript stable ID'] == alt]
		if main_exons.shape[0] > 0 and alternative_exons.shape[0] > 0:
			strand = mapper.genes.loc[main_exons.iloc[0,0], 'Strand']
			for i in main_exons.index:
				main_exon = main_exons.loc[i]
				splice_event = get_splice_events.identifySpliceEvent(main_exon, alternative_exons, strand, mapper.transcripts)
				sevents.append(splice_event)

			for i in alternative_exons.index:
				alt_exon = alternative_exons.loc[i]
				splice_event = get_splice_events.identifySpliceEvent(main_exon, alternative_exons, strand, mapper.transcripts)
				sevents.append(splice_event)
	#grab skipped exon events
	data = pd.DataFrame(sevents, columns = ['Exon', 'Event', 'Alternative Exon','Impacted Region', 'Atype', 'protein_region', 'Warnings'])
	skipped_exons = data[data['Event'] == 'skipped']
	
	#add exon sequence info to skipped exon
	skipped_exons = skipped_exons.merge(trim_exons[['Exon stable ID', 'Exon AA Seq (Full Codon)']].drop_duplicates(), left_on = 'Exon', right_on = 'Exon stable ID')
	
	#check canonical exons for ptm info first
	gp_ts_ptms = addPTMsToExons(skipped_exons, mapper, exploded_ptms, mod_groups, exon_column = 'Exon')
	return gp_ts_ptms

def getApprisDensity(mapper, trim_exons, exploded_ptms, mod_groups, figshare_dir = 'PTM_Splicing_FigShare/', fname = '../../../../FunctionalTableLong.tsv'):
	tissue_data = pd.read_csv(figshare_dir + '/External_Data/TissueSpecificExons/Rodriguez2020.tsv', sep = '\t')
	
	#identify tissue-specific exons (skipped exons in one of the tissue-specific transcripts)
	sevents = []
	for e in tissue_data.index:
		tmp_data = tissue_data.loc[e]
		main = tmp_data['Main Transcript']
		alt = tmp_data['Alternative Transcript']
		main_exons = trim_exons[trim_exons['Transcript stable ID'] == main]
		alternative_exons = trim_exons[trim_exons['Transcript stable ID'] == alt]
		if main_exons.shape[0] > 0 and alternative_exons.shape[0] > 0:
			strand = mapper.genes.loc[main_exons.iloc[0,0], 'Strand']
			for i in main_exons.index:
				main_exon = main_exons.loc[i]
				splice_event = get_splice_events.identifySpliceEvent(main_exon, alternative_exons, strand, mapper.transcripts)
				sevents.append(splice_event)

			for i in alternative_exons.index:
				alt_exon = alternative_exons.loc[i]
				splice_event = get_splice_events.identifySpliceEvent(main_exon, alternative_exons, strand, mapper.transcripts)
				sevents.append(splice_event)
	#grab skipped exon events
	data = pd.DataFrame(sevents, columns = ['Exon', 'Event', 'Alternative Exon','Impacted Region', 'Atype', 'protein_region', 'Warnings'])
	skipped_exons = data[data['Event'] == 'skipped']
	
	#add exon sequence info to skipped exon
	skipped_exons = skipped_exons.merge(trim_exons[['Exon stable ID', 'Exon AA Seq (Full Codon)']].drop_duplicates(), left_on = 'Exon', right_on = 'Exon stable ID')
	
	appris_ts_ptms = addPTMsToExons(skipped_exons, mapper, exploded_ptms, mod_groups, exon_column = 'Exon')
	return appris_ts_ptms

	
def getExonCharacteristics(exon_data, mod_type = 'All'):
    if mod_type == 'All':
        num_ptms_in_exon = exon_data.groupby('Exon stable ID')['Source of PTM'].nunique()
        num_ptms_in_exon.name = 'Number of PTMs'
    else:
        num_ptms_in_exon = exon_data.groupby(['Exon stable ID', 'Modification Class'])['Source of PTM'].nunique()
        num_ptms_in_exon.name = 'Number of PTMs'
        num_ptms_in_exon = num_ptms_in_exon.reset_index()
        num_ptms_in_exon = num_ptms_in_exon[num_ptms_in_exon['Modification Class'] == mod_type]
        num_ptms_in_exon = num_ptms_in_exon.set_index('Exon stable ID')['Number of PTMs']
    exon_density = exon_data.merge(num_ptms_in_exon, left_on = 'Exon stable ID', right_index = True, how = 'left')
    exon_density = exon_density[['Exon stable ID', 'Number of Residues', 'Number of PTMs']].drop_duplicates()
    exon_density['PTM Density'] = exon_density['Number of PTMs']/exon_density['Number of Residues']

    #get overall ptm density
    mod_density = exon_density['Number of PTMs'].sum()/exon_density['Number of Residues'].sum() * 1000
    return exon_density, mod_density

def getTSData(mapper, exploded_ptms, exploded_mod, mod_groups, figshare_dir = './PTM_Splicing_FigShare/', size_cutoff = 150):
	"""
	Extract tissue-specific exons from three different studies and calculate PTM density in these exons.

	Parameters
	----------
	mapper: PTM_mapper object
		contains information about PTMs and their genomic locations
	exploded_ptms: Pandas dataframe
		Version of mapper.ptm_info with unique PTM/transcript pairs separated into different rows
	exploded_mod: Pandas dataframe
		Version of mapper.ptm_info with unique modification types separated into different rows
	mod_groups: Pandas dataframe
		Table indicating how to convert between modification class and subtype
	figshare_dir: str
		Directory where figshare data is stored
	size_cutoff: int
		Minimum number of instances for a PTM class to be included in the analysis

	Returns
	-------
	all_data: Pandas dataframe
		Table containing information about PTMs in tissue-specific exons
	densities: Pandas dataframe
		Table containing information about PTM density in tissue-specific exons and the overall proteome

	"""
	#process exons to remove exons without coding information
	trim_exons = mapper.exons.copy()
	trim_exons['Exon Start (Protein)'] = pd.to_numeric(trim_exons['Exon Start (Protein)'], errors = 'coerce')
	trim_exons = trim_exons.dropna(subset = 'Exon Start (Protein)')
	trim_exons = trim_exons.dropna(subset = 'Exon AA Seq (Full Codon)')
	
	print('Getting PTM Density across the entire proteome')
	overall_ptm_density = getOverallDensities_Canonical(mapper, exploded_mod)
	print('Getting PTMs in tissue specific exons from Buljan et al. (2012)')
	tse_ptms = getBuljanPTMs(trim_exons, mapper, exploded_ptms, mod_groups, figshare_dir = figshare_dir)
	tse_ptms['Data Source'] = 'Buljan2012'
	print('Getting PTMs in tissue specific exons from Rodriguez et al. (2020)')
	appris_ts_ptms = getApprisDensity(mapper, trim_exons, exploded_ptms, mod_groups, figshare_dir = figshare_dir)
	appris_ts_ptms['Data Source'] = 'Rodriguez2020'
	print('Getting PTMs in tissue specific exons from Gonzalez-Porta et al. (2013)')
	gp_ts_ptms = getGP_PTMs(mapper, trim_exons, exploded_ptms, mod_groups, figshare_dir = figshare_dir)
	gp_ts_ptms['Data Source'] = 'Gonzalez-Porta2013'
	all_data = pd.concat([tse_ptms , appris_ts_ptms, gp_ts_ptms])
	all_data = all_data.drop_duplicates(subset = ['Exon stable ID', 'Source of PTM', 'Modification Class'])
	
	#overall_ptm_density = getOverallDensities(mapper, exploded_mod)
	
	#get overall densities
	sizes = exploded_mod.groupby('Modification Class').size()
	sizes = sizes[sizes > size_cutoff].sort_values(ascending = False)
	mods_to_keep = sizes.index
	#plt_data1 = mod_densities[mods_to_keep].reset_index()
	#plt_data1['Type'] = 'Entire Proteome'
	#overall_proteome = pd.DataFrame({'Modification Class': 'All', 0:overall_ptm_density,'Type':'Entire Proteome'}, index = [300])
	#plt_data1 = pd.concat([overall_proteome, plt_data1])
	#plt_data1.index = plt_data1['Modification Class']
	
	print('Calculating PTM density in tissue-specific exons')
	#get number of aa in each exon
	all_data['Number of Residues'] = all_data['Exon AA Seq (Full Codon)'].apply(len)
	all_data2 = all_data[['Exon stable ID', 'Source of PTM', 'Number of Residues', 'Modification Class']].copy()
	mod_density = {}
	for mod in mods_to_keep:
		exon_density, mod_density[mod] = getExonCharacteristics(all_data2, mod_type = mod)
	ts_density = pd.Series(mod_density)
	#get density for all ptms
	all_data['Exon Length'] = all_data['Exon AA Seq (Full Codon)'].apply(len)
	num_aa_tse = all_data.drop_duplicates('Exon stable ID')['Exon Length'].sum()
	ts_density['All'] = (all_data['Source of PTM'].nunique()/num_aa_tse)*1000
	ts_density.name = 'Tissue-Specific Density (per 1000aa)'

	densities = pd.concat([overall_ptm_density, ts_density], axis = 1).loc[['All']+list(mods_to_keep)]
	densities['Normalized PTM Density'] = densities['Tissue-Specific Density (per 1000aa)']/densities['Density in Proteome (per 1000aa)']

	#get ts densities, normalize to overall
	#all_data['Exon Length'] = all_data['Exon AA Seq (Full Codon)'].apply(len)
	#num_aa_tse = all_data.drop_duplicates('Exon stable ID')['Exon Length'].sum()
	#tse_mods = [mod for mod in mods_to_keep if mod in all_data['Modification Class'].values]
	#plt_data2 = ((all_data.groupby('Modification Class').size()/num_aa_tse)*1000)[tse_mods].reset_index()
	#overall_tse = pd.DataFrame({'Modification Class': 'All', 0:(all_data['Source of PTM'].nunique()/num_aa_tse)*1000}, index = [100])
	#plt_data2 = pd.concat([overall_tse,plt_data2])
	#plt_data2.index = plt_data2['Modification Class']
	#plt_data1 = plt_data1.loc[['All'] + tse_mods]
	#plt_data2[0] = plt_data2.loc[plt_data1.index, 0]/ plt_data1[0]
	
	return all_data, densities