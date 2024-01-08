import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
import matplotlib.colors as mcl
import matplotlib.patches as patches
import networkx as nx



def stackedBar(xlabels, ylist, colormap = 'inferno', width = 0.5, x_label = '',y_label = '', title = '',figsize = (3,3),
            ax = None, legend = None, leg_loc = [1.05,0.8], title_fontsize = 10, legend_fontsize = 9, legend_title = 'Event Type'):
    """
    Provided a 2D list of y-values, plots a stacked bar plot.

    Parameters
    ----------
    xlabels : list
        labels for each of the bars in the plot (x variables)
    ylist : Nested list
        Contains the needed data for each xvariable. In the format [[values for bar1],[values for bar2],...]
    colormap : string, optional
        Idea is to select color gradient to choose from. Currently does not change anything however. 
        The default is 'BuGn'.
    width : float, optional
        width of the bars to plot. The default is 0.5.
    y_label : string, optional
        Y axis label. The default is ''.
    title : string, optional
        title for the graph. The default is ''.
    legend : list of strings, optional
        Labels to use for the legend, which will indicate  what each bar color in a stack indicates. 
        The default is None.

    Returns
    -------
    Plots a stacked bar plot with each of the categories in a different color. Right now, color scheme
    will need to be manually changed in function.

    """
    #choose colors
    if isinstance(colormap, str):
        if colormap == 'terrain':
            colors = cm.terrain(np.linspace(0,1,len(ylist)+3))
        elif colormap == 'inferno':
            colors = cm.inferno(np.linspace(0,1,len(ylist)))
        elif colormap == 'colorblind':
            if len(ylist) == 5:
                colors = sns.color_palette('colorblind')
                colors = [colors[0],colors[1],colors[2],colors[4],colors[5]]
            else:
                colors = sns.color_palette('colorblind', n_colors = len(ylist))
    elif isinstance(colormap, list) or isinstance(colormap, np.ndarray):
        #check to make sure length of colormap is correct
        if len(colormap) == len(ylist):
            colors = colormap
        else:
            raise ValueError('Length of colormap list does not match number of bars to plot')
    else:
        colors = cm.nipy_spectral(np.linspace(0,1,len(ylist)+2))
    
    #if no axes provided, create figure and axes
    if ax == None:
        fig = plt.figure(figsize = figsize)
        ax = fig.add_axes([0,0,1,1])
    #plot first set of bars
    ax.bar(xlabels, ylist[0], color=colors[0])
    #plot second set of bars on top first
    ax.bar(xlabels, ylist[1], bottom = ylist[0], color=colors[1])
    #continue with remaining values in ylist
    if len(ylist) > 2:
        bottom = [ylist[0][j]+ylist[1][j] for j in range(len(ylist[2]))]
        ax.bar(xlabels, ylist[2],bottom=bottom, color=colors[2])
    if len(ylist) > 3:
        for i in range(3,len(ylist)):
            bottom = [bottom[j] + ylist[i-1][j] for j in range(0,len(ylist[i]))]
            ax.bar(xlabels, ylist[i], bottom = bottom, color=colors[i])
    
    #add labels, ticks, and title
    ax.set_ylabel(y_label, fontsize = 12)
    ax.set_xlabel(x_label, fontsize = 12)
    ax.set_title(title)
    ax.set_yticks(np.arange(0, 1.1, 0.2))
    #
    if isinstance(legend, list):
        ax.legend(labels=legend, loc = leg_loc, title_fontsize = title_fontsize, fontsize = legend_fontsize, ncols = 1, title = legend_title)
    plt.xticks(rotation = 35, ha = 'right')
    
    
def stackedBarH(xlabels, ylist, colormap = 'inferno', width = 0.5, x_label = '',
                ax = None, y_label = '', title = '',figsize = (3,3), legend = None, leg_loc = [1.05,0.8], legend_type = 'vertical', legend_fontsize = 9, legend_title = '', legend_bbox_to_anchor = None):
    """
    Given a 2D list of y-values, plots a horizontal stacked bar plot.

    Parameters
    ----------
    xlabels : list
        labels for each of the bars in the plot (x variables)
    ylist : Nested list
        Contains the needed data for each xvariable. In the format [[values for bar1],[values for bar2],...]
    colormap : string, optional
        Idea is to select color gradient to choose from. Currently does not change anything however. 
        The default is 'BuGn'.
    width : float, optional
        width of the bars to plot. The default is 0.5.
    y_label : string, optional
        Y axis label. The default is ''.
    title : string, optional
        title for the graph. The default is ''.
    legend : list of strings, optional
        Labels to use for the legend, which will indicate  what each bar color in a stack indicates. 
        The default is None.
    legend_type : string, optional
        Direction to include legend labels in. Options are 'vertical' or 'horizontal'. The default is 'vertical'.
    legend_fontsize : int, optional
        Fontsize for legend labels. The default is 9.
    legend_title : string, optional
        Title for legend. The default is ''.

    Returns
    -------
    Plots a stacked bar plot with each of the categories in a different color. Right now, color scheme
    will need to be manually changed in function.

    """
    #pick the color scheme
    if isinstance(colormap, str):
        if colormap == 'terrain':
            colors = cm.terrain(np.linspace(0,1,len(ylist)+3))
        elif colormap == 'inferno':
            colors = cm.inferno(np.linspace(0,1,len(ylist)))
        elif colormap == 'colorblind':
            colors = sns.color_palette('colorblind')
    elif isinstance(colormap, list) or isinstance(colormap, np.ndarray):
        colors = colormap
    else:
        colors = cm.nipy_spectral(np.linspace(0,1,len(ylist)+2))

    #if axes is not provided, create figure
    if ax == None:
        fig = plt.figure(figsize = figsize)
        ax = fig.add_axes([0,0,1,1])
    
    #plot the first set of bars
    ax.barh(xlabels, ylist[0], width, color=colors[0])
    #plot the second set of bars to the right of the first
    ax.barh(xlabels, ylist[1], width, left = ylist[0], color=colors[1])
    #continue with this pattern, adding to the previous bar until no more values in ylist
    if len(ylist) > 2:
        bottom = [ylist[0][j]+ylist[1][j] for j in range(len(ylist[2]))]
        ax.barh(xlabels, ylist[2], width,left=bottom, color=colors[2])
    if len(ylist) > 3:
        for i in range(3,len(ylist)):
            bottom = [bottom[j] + ylist[i-1][j] for j in range(0,len(ylist[i]))]
            ax.barh(xlabels, ylist[i], width, left = bottom, color=colors[i])
    
    #set the labels, ticks,  and title
    ax.set_xlabel(y_label, fontsize = 10)
    ax.set_ylabel(x_label, fontsize = 10)
    ax.set_title(title)
    ax.set_xticks(np.arange(0, 1.1, 0.2))
    ax.set_xlim([0,1])
    #add legend if provided
    if isinstance(legend, list):
        if legend_type == 'vertical':
            if legend_bbox_to_anchor is None:
                ax.legend(labels=legend, loc = leg_loc, fontsize = legend_fontsize, title = legend_title)
            else:
                ax.legend(labels=legend, loc = leg_loc, fontsize = legend_fontsize, title = legend_title, bbox_to_anchor = legend_bbox_to_anchor)
        elif legend_type == 'horizontal':
            if legend_bbox_to_anchor is None:
                ax.legend(labels=legend, loc = leg_loc, fontsize = legend_fontsize, title = legend_title, ncol = len(legend))
            else:
                ax.legend(labels=legend, loc = leg_loc, fontsize = legend_fontsize, title = legend_title, ncol = len(legend), bbox_to_anchor = legend_bbox_to_anchor)
        else:
            raise ValueError('legend_type must be either "vertical" or "horizontal"')
    plt.yticks(rotation = 0)
    return ax
    
def addStatAnnot(p, x1, x2, ax, h = 2, col = 'k', equal_var = False, start_height = None):
    """
    Adds statistical annotation to a plot.

    Parameters
    ----------
    p : float
        p-value to be annotated.
    x1 : float
        x-coordinate of the first bar.
    x2 : float
        x-coordinate of the second bar.
    ax : matplotlib axes object
        axes object to add annotation to.
    h : float, optional
        height of the annotation bar. The default is 2.
    col : string, optional
        color of the annotation bar. The default is 'k'.
    equal_var : bool, optional
        whether or not to assume equal variance. The default is False.
    start_height : float, optional
        height to start annotation bar at. If none, will autocalculate to be the maximum y value in data. The default is None. 
    
    Returns
    -------
    p : float
        p-value to be annotated.
    
    """
    # statistical annotation # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
    if start_height is None:
        y = data[y].max() + h
    else:
        y = start_height
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1, c=col)
    #if p < 1e-10:
    #    annot = '***'
    #elif p < 1e-5:
    #    annot = '**'
    #elif p < 0.05:
    #    annot = '*'

    if p < 0.05:
        annot = '*'
    else:
        annot = 'n.s.'
    ax.text((x1+x2)*.5, y+0.01, annot, ha='center', va='bottom', color=col, fontsize = 11)
    return p

def plotSIRT3_Network(canonical_scores, alternative_scores):
    """
    """
    #extract information for SIRT3, and combine canonical and alternative information
    sirt3_scores = alternative_scores[alternative_scores['PTM'] == 'Q9NTG7_S159']
    sirt3_scores = sirt3_scores.drop('PTM', axis = 1)
    sirt3_scores = sirt3_scores.rename({'site_percentile':'Alternative SLiM Percentile'}, axis = 1)
    sirt3_scores2 = canonical_scores[canonical_scores['PTM'] == 'Q9NTG7_S159']
    sirt3_scores2 = sirt3_scores2.drop(['PTM', 'Known Interaction'], axis = 1)
    sirt3_scores2 = sirt3_scores2.rename({'site_percentile':'Canonical SLiM Percentile'}, axis = 1)
    sirt3_scores = sirt3_scores.merge(sirt3_scores2, on = 'kinase')
    sirt3_scores = sirt3_scores[['kinase', 'Known Interaction', 'Canonical SLiM Percentile', 'Alternative SLiM Percentile', 'Change']]

    #grab scores for sirt3 S159, restrict to interactions with CDKs
    sirt3_scores = alternative_scores[alternative_scores['PTM'] == 'Q9NTG7_S159']
    sirt3_scores = sirt3_scores[(sirt3_scores['kinase'].str.contains('CDK')) & 
                ~(sirt3_scores['kinase'].str.contains('CDKL')) &
                (sirt3_scores['Change'].apply(abs) > 0)]
    sirt3_scores['substrate'] = 'SIRT3\nS159'

    #setup subplot
    fig, ax = plt.subplots(figsize = (3,2.5), nrows = 2, height_ratios = [1,0.05])
    fig.subplots_adjust(wspace = 0.02)

    #create directed graph and add edges
    G = nx.DiGraph()
    for i, row in sirt3_scores.iterrows():
        G.add_edge(row['kinase'],row['substrate'],weight=row['Change'])

    #color edges based on sign of difference (red for increase, blue for decrease)
    #normalizer = mcl.Normalize(vmin = -3, vmax = 3)
    normalizer = mcl.Normalize(vmin = -20, vmax = 20)
    cmap = cm.ScalarMappable(normalizer,cmap = 'coolwarm')
    # Define the edge colors based on the sign of the weight
    edge_colors = [cmap.to_rgba(w) for u, v, w in G.edges.data('weight')]

    # Define the edge widths based on the absolute value of the weight
    edge_widths = [2 for u, v, w in G.edges.data('weight')]

    # Draw the graph with matplotlib
    pos = nx.circular_layout(G)
    pos['SIRT3\nS159'] = np.array([0,0])
    nodes = nx.draw_networkx_nodes(G, pos, node_size = 500, node_color = 'lightgray', ax = ax[0])
    nx.draw_networkx_labels(G, pos, font_size = 7, ax = ax[0])
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, ax = ax[0])
    fig.colorbar(cmap,cax=ax[1], orientation='horizontal', label='Predicted Change in Kinase Percentile')
    ax[0].axis('off')

def plotTranscript(mapper, transcript_id, ptm_loc = None, ax = None, trans_type = 'Canonical',coding_color = 'red', noncoding_color = 'white', fontsize = 18):
    gene_id = mapper.transcripts.loc[transcript_id, 'Gene stable ID']
    transcript_exons = mapper.exons[mapper.exons['Transcript stable ID'] == transcript_id]

    #get relevant transcript info
    transcript = mapper.transcripts.loc[transcript_id]

    #establish plot range   
    if ax is None:
        fig, ax = plt.subplots(figsize = (3,0.5))
    
    #transcript box
    gene_start = mapper.genes.loc[gene_id, 'Gene start (bp)']
    gene_end = mapper.genes.loc[gene_id, 'Gene end (bp)']
    strand = mapper.genes.loc[gene_id, 'Strand']
    ax.set_xlim([gene_start-2500,gene_end+2500])

    #get location of cds start and cds end
    gene_cds_start, gene_cds_end = getGenomicCodingRegion(mapper.transcripts.loc[transcript_id], transcript_exons, mapper.genes.loc[gene_id,'Strand']) 
    if strand == 1:
        ax.annotate(transcript_id, (gene_start-1500, 0.6), ha = 'right', va = 'center', fontsize =fontsize)
        ax.annotate(f'({trans_type})', (gene_start-1500, 0.4), ha = 'right', va = 'center', fontsize = fontsize)
    else:
        ax.annotate(transcript_id, (gene_end+1500, 0.6), ha = 'right', va = 'center', fontsize = fontsize)
        ax.annotate(f'({trans_type})', (gene_end+1500, 0.4), ha = 'right', va = 'center', fontsize = fontsize)
    ax.plot([gene_start, gene_end], [0.5, 0.5], c = 'k')
    for i,row in transcript_exons.iterrows():
        fiveprime = int(row[ "Exon Start (Gene)"])
        threeprime = int(row["Exon End (Gene)"])
        #check if exon is fully noncoding, fully coding, or if cds start/stop exists in exon. Plot exon accordingly
        if threeprime < gene_cds_start or fiveprime > gene_cds_end:
            rect = patches.Rectangle((fiveprime,0.2), threeprime - fiveprime, 0.6, facecolor = noncoding_color, edgecolor = 'black', zorder = 2)
            ax.add_patch(rect)
        elif fiveprime >= gene_cds_start and threeprime <= gene_cds_end:
            rect = patches.Rectangle((fiveprime,0.2), threeprime - fiveprime, 0.6, facecolor = coding_color, edgecolor = 'black', zorder = 2)
            ax.add_patch(rect)
        else:
            noncoding_rect = patches.Rectangle((fiveprime,0.2), threeprime - fiveprime, 0.6, facecolor = noncoding_color, edgecolor = 'black', zorder = 2)
            ax.add_patch(noncoding_rect)

            if fiveprime < gene_cds_start:
                fiveprime = gene_cds_start
            if threeprime > gene_cds_end:
                threeprime = gene_cds_end

            noncoding_rect = patches.Rectangle((fiveprime,0.2), threeprime - fiveprime, 0.6, facecolor = coding_color, edgecolor = 'black', zorder = 3)
            ax.add_patch(noncoding_rect)
            
    if ptm_loc is not None:
        ax.plot([ptm_loc, ptm_loc], [0.8,0.9], c = 'black')
        circle = patches.Ellipse((ptm_loc, 0.95), height = 0.2, width = 500, color = 'gold', 
                                 ec = 'black', zorder = 2)
        ax.add_artist(circle)
        ax.annotate('P',(ptm_loc, 0.95), va = 'center', ha='center', fontsize = 9) 
        ax.set_ylim([0,1.1])
    else:
        ax.set_ylim([0,1])

    if strand == -1:
        ax.invert_xaxis()

    ax.axis('off')
    
def compareTranscripts(mapper, canonical_transcript, alternative_transcript,ptm = None, fontsize = 11, gs = None, fig = None):
    """
    Plot canonical and alternative transcripts together on same axis

    Parameters
    ----------
    mapper : PTM_mapper
        PTM_mapper object
    canonical_transcript : str
        Ensembl Transcript ID of canonical transcript
    alternative_transcript : str
        Ensembl Transcript ID of alternative transcript
    ptm : str, optional
        PTM to highlight. The default is None.

    Returns
    -------
    None. Plots figure.
    """
    if gs is None:
        fig, ax = plt.subplots(figsize = (10,2.5), nrows = 2)
        fig.subplots_adjust(hspace = 0)
    else:
        #construct subplots within gridspec
        gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs, hspace = 0)
        ax = [fig.add_subplot(gs2[0]), fig.add_subplot(gs2[1])]
    
    if ptm is not None:
        ptm = mapper.ptm_coordinates[mapper.ptm_coordinates['Source of PTM'] == ptm]
        ptm_loc = ptm['Gene Location (NC)']
    else:
        ptm_loc = None
        
    
    plotTranscript(mapper,canonical_transcript,ax = ax[0], ptm_loc = ptm_loc, trans_type = 'Canonical', fontsize = fontsize)
    plotTranscript(mapper,alternative_transcript,ax=ax[1], ptm_loc = ptm_loc, trans_type = 'Alternative', fontsize = fontsize)
    ax[0].set_title('SIRT3', fontsize = fontsize)

def getGenomicCodingRegion(transcript, exons, strand):
    """
    Get genomic coordinates of coding region of transcript

    Parameters
    ----------
    transcript : pd.Series
        Series containing transcript information (row from mapper.transcripts dataframe)
    exons : pd.DataFrame
        Dataframe containing exon information (mapper.exons)
    strand : int
        Strand of gene: 1 for forward, -1 for reverse

    Returns
    -------
    gene_cds_start : int
        Genomic coordinate of start of coding region
    gene_cds_end : int
        Genomic coordinate of end of coding region
    """
    #get CDS Start/Stop and map to location in gene if within exon range
    cds_start = int(transcript['Relative CDS Start (bp)'])
    cds_stop = int(transcript['Relative CDS Stop (bp)'])
    coding_start_exon = exons[(exons['Exon Start (Transcript)'] <= cds_start) & (exons['Exon End (Transcript)'] >= cds_start)].squeeze()
    coding_stop_exon = exons[(exons['Exon Start (Transcript)'] <= cds_stop) & (exons['Exon End (Transcript)'] >= cds_stop)].squeeze()
    if strand == 1:
        gene_cds_start = coding_start_exon['Exon Start (Gene)'] + (cds_start - coding_start_exon['Exon Start (Transcript)'])
        gene_cds_end = coding_stop_exon['Exon Start (Gene)'] + (cds_stop - coding_stop_exon['Exon Start (Transcript)'])
    else:
        gene_cds_end = coding_start_exon['Exon End (Gene)'] - (cds_start - coding_start_exon['Exon Start (Transcript)'])
        gene_cds_start = coding_stop_exon['Exon End (Gene)'] - (cds_stop - coding_stop_exon['Exon Start (Transcript)'])
        
    return gene_cds_start, gene_cds_end

def range_normalize(x, min_val, max_val, desired_range = [0,1]):
    """
    Normalize a value to a desired range

    Parameters
    ----------
    x : float
        value to normalize
    min_val : float
        minimum value of range to normalize to
    max_val : float
        maximum value of range to normalize to
    desired_range : list, optional
        desired range to normalize to. The default is [0,1].

    Returns
    -------
    float
        normalized value

    """
    return (x-min_val)/(max_val-min_val)*(desired_range[1]-desired_range[0])+desired_range[0]

def plot_exon(mapper, splice_seq, PRAD_ptms, gene_name, exon_num = 4.1, ax = None, title = True, normed_loc = None, include_coordinate_axis = True):
    """
    Given a gene name and exon number, plot the exon and any PTMs found in the exon

    Parameters
    ----------
    mapper : PTM mapper object
        PTM mapper object containing PTM coordinate data
    splice_seq : pandas dataframe
        dataframe containing exon coordinate data
    PRAD_ptms : pandas dataframe
        dataframe containing PTM data mapped to splice seq exons in splice_seq dataframe
    gene_name : str
        gene name to plot
    exon_num : float, optional
        exon number to plot. The default is 4.1.

    Returns
    -------
    None.

    Outputs
    -------
    plot of exon with PTMs
    """
    plt_data = splice_seq[splice_seq['Symbol'] == gene_name]
    plt_data = plt_data[plt_data['Exon'] == exon_num].squeeze()
    min_loc = plt_data['Chr_Start'].min() - 10
    max_loc = plt_data['Chr_Stop'].max() + 10

    if ax is None:
        fig, ax = plt.subplots(figsize = (6,1))
    if normed_loc is None:
        ax.set_xlim([min_loc, max_loc])
    #else:
    #    plt.xlim(normed_loc)
    #create rectangle object that starts at start and ends at stop
    start = plt_data['Chr_Start'] if normed_loc is None else range_normalize(plt_data['Chr_Start'], min_loc, max_loc, desired_range=normed_loc)
    stop = plt_data['Chr_Stop'] if normed_loc is None else range_normalize(plt_data['Chr_Stop'], min_loc, max_loc, desired_range=normed_loc)

    rect = patches.Rectangle((start, 0.4),stop-start,0.2,facecolor='red', linewidth = 1, edgecolor = 'k')
    # Add the patch to the Axes
    ax.add_patch(rect)

    #add ptms from exon to plot
    PRAD_ptms = PRAD_ptms[PRAD_ptms['Adj p'] <= 0.05]
    PRAD_ptms = PRAD_ptms[PRAD_ptms['symbol'] == gene_name]
    PRAD_ptms = PRAD_ptms[PRAD_ptms['Exon'] == exon_num]
    PRAD_ptms = PRAD_ptms.drop_duplicates(subset = 'PTM')
    #sort to make sure first residues are plotted first
    PRAD_ptms['Position'] = PRAD_ptms['PTM'].apply(lambda x: int(x.split('_')[1][1:]))
    PRAD_ptms = PRAD_ptms.sort_values(by = 'Position', ascending = True)
    for index, i in zip(PRAD_ptms.index, range(PRAD_ptms.shape[0])):
        ptm = mapper.ptm_coordinates[mapper.ptm_coordinates['Source of PTM'] == PRAD_ptms.loc[index, 'PTM']].squeeze()
        ptm_loc = ptm['HG19 Location'] if normed_loc is None else range_normalize(ptm['HG19 Location'], min_loc, max_loc, desired_range=normed_loc)
        #extract mod type
        mod_type = ptm['Modifications']
        if 'Phospho' in mod_type:
            color = 'gold'
        elif 'Glyco' in mod_type:
            color = 'lightpink'
        elif 'Methyl' in mod_type or 'methyl' in mod_type:
            color = 'lightblue'
        elif 'Ubiquitination' in mod_type:
            color = 'orange'
        elif 'Acetyl' in mod_type or 'acetyl' in mod_type:
            color = 'lightgreen'
        elif 'Sumo' in mod_type:
            color = 'brown'
        else:
            color = 'lightgrey'

        if i%2 == 0:
            #ax.plot([ptm_loc, ptm_loc], [0.5, 0.7], color = 'k', linewidth = 1)
            ax.plot(ptm_loc, 0.7, 'o', color = color, markersize = 10, mec='k', linewidth = 0.5)
                #annotate ptm circle with residue
            ax.annotate(ptm['Residue'], xy = (ptm_loc, 0.65), ha = 'center', va = 'center', fontsize = 9)
            ax.annotate(ptm['Source of PTM'].split('_')[1][1:], xy = (ptm_loc, 1), ha = 'center', va = 'center', fontsize = 7)
        else:
            #ax.plot([ptm_loc, ptm_loc], [0.2, 0.4], color = 'k', linewidth = 1)
            ax.plot(ptm_loc, 0.3, 'o', color = color, markersize = 10, mec='k', linewidth = 0.5)
            ax.annotate(ptm['Residue'], xy = (ptm_loc, 0.25), ha = 'center', va = 'center', fontsize = 9)
            ax.annotate(ptm['Source of PTM'].split('_')[1][1:], xy = (ptm_loc, 0), ha = 'center', va = 'center', fontsize = 7)

    ax.set_ylim([0,1])
    ax.set_yticks([])
    if include_coordinate_axis:
        ax.set_xlabel('Genomic Location')
        #keep only bottom axis line
        plt.gca().spines['left'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        #add line to signal gene
        ax.axhline(y=0.5, color='k', zorder = 0, linewidth = 1)
    else:
        ax.set_xticks([])
        ax.set_xlabel('')
        ax.axis('off')

    #add title
    if title:
        ax.set_title(f'{gene_name}: Exon {exon_num:.1f}', fontsize = 11)

def PSI_boxplots_SpecificGenes(mapper, splice_seq, cancer_ptms, ESRP1_low_id, ESRP1_high_id, gene_and_exon = [('CD44',3.2),('CTNND1', 4.3), ('CTNND1',21)],
                 max_p = 0.05, min_effect_size = 0.25, figsize = (3,3.5), gs = None, fig = None):
    """

    """
    cancer_ptms2 = cancer_ptms[cancer_ptms['Adj p'] <= max_p]
    cancer_ptms2 = cancer_ptms2[cancer_ptms2['Effect Size'] >= min_effect_size]

    #restrict ids to ones found in ptm dataframe
    ESRP1_low_id = [i for i in ESRP1_low_id if i in cancer_ptms2.columns]
    ESRP1_high_id = [i for i in ESRP1_high_id if i in cancer_ptms2.columns]

    plt_data = []
    for ge in gene_and_exon:
        #extract information for specific gene and exon combo
        tmp = cancer_ptms2[cancer_ptms2['symbol'] == ge[0]]
        tmp = tmp[tmp['exons'] == ge[1]]
        tmp = tmp.drop_duplicates()
        #grab a single ptm entry
        tmp = tmp.iloc[0]

        #convert to dataframe with necessary context
        row_data = pd.DataFrame({'PSI': list(tmp[ESRP1_low_id])+list(tmp[ESRP1_high_id]), 
                    'Type': list(np.repeat('Low', len(ESRP1_low_id)))+list(np.repeat('High', len(ESRP1_high_id))),
                    'Protein': np.repeat(tmp['symbol'], len(ESRP1_low_id)+len(ESRP1_high_id)),
                    'PTM': np.repeat(tmp['PTM'], len(ESRP1_low_id)+len(ESRP1_high_id)),
                        'p': np.repeat(tmp['Adj p'], len(ESRP1_low_id)+len(ESRP1_high_id)),
                        'Exon': np.repeat(str((tmp['exons'])), len(ESRP1_low_id)+len(ESRP1_high_id))})
        plt_data.append(row_data)
    #combine data into one datframe and add Label
    plt_data = pd.concat(plt_data)
    plt_data['Label'] = plt_data['Protein'] + '\n' + 'Exon ' + plt_data['Exon']

    #plot boxplots
    if gs is None:
        fig, ax = plt.subplots(figsize = figsize, nrows = 2, height_ratios=[1,0.25])
        fig.subplots_adjust(hspace = 0.5)
    else:
        #construct subplots within gridspec
        gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs, hspace = 0.5, height_ratios = [1,0.25])
        ax = [fig.add_subplot(gs2[0]), fig.add_subplot(gs2[1])]

    box = sns.boxplot(x = 'Label', y = 'PSI', hue = 'Type', data = plt_data, ax = ax[0])
    box.legend_.remove()
    plt.xticks(fontsize = 9)
    ax[0].set_xlabel(' ')
    ax[0].set_ylabel('Percent Spliced In (PSI)')

    location = np.linspace(0,1,len(gene_and_exon)+1)
    #add exons to plot
    for i, ge in enumerate(gene_and_exon):
        plot_exon(mapper, splice_seq, cancer_ptms, ge[0], exon_num = ge[1], normed_loc=[location[i],location[i+1]], ax = ax[1], title = False, include_coordinate_axis = False)

def PSI_boxplots_BiologicalProcess(mapper, splice_seq, cancer_ptms, ESRP1_low_id, ESRP1_high_id, process_table, process_name,
                 max_p = 0.05, min_effect_size = 0.25, figsize = (3,3.5), gs = None, fig = None):
    """
    """
    #restrict process table to ptms involved in process of interest
    process_table = process_table.dropna(subset = process_table.columns[process_table.columns.str.contains(process_name)], how = 'all')

    cancer_ptms2 = cancer_ptms[cancer_ptms['Adj p'] <= max_p]
    cancer_ptms2 = cancer_ptms2[cancer_ptms2['Effect Size'] >= min_effect_size]
    cancer_ptms2 = cancer_ptms2.drop_duplicates(subset = 'PTM')

    induces = []
    inhibits = []
    altered = []
    for ptm, row in process_table.iterrows():
        if row[f'{process_name}, induced'] == 1:
            induces.append(cancer_ptms2[cancer_ptms2['PTM'] == ptm].squeeze())
        elif row[f'{process_name}, inhibited'] == 1:
            inhibits.append(cancer_ptms2[cancer_ptms2['PTM'] == ptm].squeeze())
        elif row[f'{process_name}, altered'] == 1:
            altered.append(cancer_ptms2[cancer_ptms2['PTM'] == ptm].squeeze())

    plt_data = {}

    for type,name in zip([induces, inhibits, altered], ['Induces', 'Inhibits', 'Altered']):
        type_data = []
        if len(type) > 0:
            for row in type:
                row_data = pd.DataFrame({'Percent Spliced In (PSI)': list(row[ESRP1_low_id])+list(row[ESRP1_high_id]), 
                            'Type': list(np.repeat('Low', len(ESRP1_low_id)))+list(np.repeat('High', len(ESRP1_high_id))),
                            'Protein': np.repeat(row['symbol'], len(ESRP1_low_id)+len(ESRP1_high_id)),
                            'PTM': np.repeat(row['PTM'], len(ESRP1_low_id)+len(ESRP1_high_id)),
                            'p': np.repeat(row['Adj p'], len(ESRP1_low_id)+len(ESRP1_high_id)),
                            'Exon': np.repeat(str(float(row['exons'])), len(ESRP1_low_id)+len(ESRP1_high_id))})
                type_data.append(row_data)
            type_data = pd.concat(type_data)
            type_data['Label'] = type_data['Protein'] + '\n' + 'Exon ' + type_data['Exon']
            plt_data[name] = type_data 

    
    #setup subplots
    #calculate width ratios
    total_num_exons = len(inhibits) + len(induces) + len(altered)
    width_ratios = [len(type)/total_num_exons for type in [induces, inhibits, altered] if len(type) > 0]

    if gs is None:
        fig, ax = plt.subplots(figsize = figsize, nrows = 2, ncols = len(plt_data), width_ratios = width_ratios, height_ratios=[1,0.25], sharey = 'row')
        fig.subplots_adjust(hspace = 0.5, wspace = 0.1)
    else:
        #construct subplots within gridspec
        gs2 = gridspec.GridSpecFromSubplotSpec(2, len(plt_data), subplot_spec=gs, hspace = 0.5, wspace = 0.1, width_ratios = width_ratios, height_ratios = [1,0.25])
        ax = np.array([[fig.add_subplot(gs2[0,i]) for i in range(len(plt_data))],[fig.add_subplot(gs2[1,i]) for i in range(len(plt_data))]])
        ax[0,0].set_ylim([0,1])
        ax[0,1].set_ylim([0,1])
        ax[0,1].set_yticks([])

    titles = ['Induced', 'Inhibited', 'Altered']
    for i, type in zip(range(len(plt_data)), plt_data.keys()):
        if len(plt_data) > 1:
            box = sns.boxplot(x = 'Label', y = 'Percent Spliced In (PSI)', hue = 'Type', data = plt_data[type], ax = ax[0,i])
            ax[0,i].set_xlabel('')
            if i != 0:
                ax[0,i].set_ylabel('')
            ax[0,i].set_title(type, fontsize = 11)
            box.legend_.remove()
            #plt.title('Genes Regulated by ESRP1', fontsize = 11)
            plt.xticks(fontsize = 9)
            ax[0, 0].set_xlabel(' ')

            #iterate through each and plot exon
            exon_size = 1/plt_data[type]['Exon'].nunique()
            prev_start = -exon_size
            for ptm in plt_data[type]['PTM'].unique():
                new_start = prev_start + exon_size
                row = plt_data[type][plt_data[type]['PTM'] == ptm].iloc[0]
                plot_exon(mapper, splice_seq, cancer_ptms[cancer_ptms['PTM'] == ptm], row['Protein'], float(row['Exon']), normed_loc = [new_start, new_start+exon_size], include_coordinate_axis = False, title = False, ax = ax[1,i])
                prev_start = new_start
        else:
            box = sns.boxplot(x = 'Label', y = 'Percent Spliced In (PSI)', hue = 'Type', data = plt_data[type], ax = ax[0])
            ax[0].set_xlabel('')
            ax[0].set_title(type, fontsize = 11)
            box.legend_.remove()
            #plt.title('Genes Regulated by ESRP1', fontsize = 11)
            plt.xticks(fontsize = 9)

            #iterate through each and plot exon
            exon_size = 1/plt_data[type]['Exon'].nunique()
            prev_start = -exon_size
            for ptm in plt_data[type]['PTM'].unique():
                new_start = prev_start + exon_size
                row = plt_data[type][plt_data[type]['PTM'] == ptm].iloc[0]
                plot_exon(mapper, splice_seq, cancer_ptms[cancer_ptms['PTM'] == ptm], row['Protein'], float(row['Exon']), normed_loc = [new_start, new_start+exon_size], include_coordinate_axis = False, title = False, ax = ax[1])
                prev_start = new_start

    fig.suptitle(f'PTMs related to {process_name}', y = 1.01)

def visualize_exons(mapper, trans_id):
    """
    Given transcript ID, plot the full transcript broken up by each exon (with an exon being a distinct rectangle) and annotated with exon number
    """

    exons = mapper.exons[mapper.exons['Transcript stable ID'] == trans_id]
    length = exons['Exon End (Transcript)'].max()

    #create figure
    fig, ax = plt.subplots(figsize = (15,1))
    ax.axis('off')
    ax.set_xlim([0,length])
    ax.set_ylim([0,1])

    #plot exons
    for i, row in exons.iterrows():
        exon_start = row['Exon Start (Transcript)']
        exon_end = row['Exon End (Transcript)']
        exon_num = row['Exon rank in transcript']
        color = 'lightblue' if exon_num != 26 else 'red'
        ax.add_patch(patches.Rectangle((exon_start,0.2), exon_end - exon_start, 0.6, facecolor = color, edgecolor = 'black', linewidth = 2))
        #annotate with exon number, but if exon length is too small, annotate with exon number above exon
        if exon_end - exon_start < 0.02*length:
            ax.annotate(exon_num, (exon_start + (exon_end - exon_start)/2 ,0.9), ha = 'center', va = 'center', fontsize = 10)
        else:
            ax.annotate(exon_num, (exon_start + (exon_end - exon_start)/2,0.5), ha = 'center', va = 'center', fontsize = 10)


def high_low_colorbar(expression, gene_name = None, cutoff = 1, fontsize = 10, xlabel = 'TCGA Prostate Adenocarcinoma Patient', ax = None, figsize = (3, 0.5)):
    #sort expression data from low to high
    expression = expression.sort_values(ascending = True)

    #make colorbar from scratch with coolwarm palette, ranging from -5 to 5
    cmap = plt.cm.coolwarm
    norm = mcl.Normalize(vmin=-4, vmax=4)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

    if ax is None:
        fig,ax = plt.subplots(figsize = figsize)

    ax.tick_params(axis='x', labelsize=fontsize -1)
    ax.plot(0,0)
    ax.set_ylim([0,1])
    ax.set_xlim([0,expression.shape[0]])
    low_line = False
    high_line = False
    for i, val in enumerate(expression):
        rect = patches.Rectangle((i,0.05),1,1, color = sm.to_rgba(val))
        ax.add_patch(rect)

        if val > -1*cutoff and not low_line:
            ax.axvline(i, color = 'black')
            low_line = True
            #record where first patient not hitting cutoff is
            low_patient_cutoff = i

        if val > cutoff and not high_line:
            ax.axvline(i, color = 'black')
            high_line = True
            #record where first patient not hitting cutoff is
            high_patient_cutoff = i

    #annotate high and low regions
    ax.annotate('Low', (low_patient_cutoff/2, 0.5), ha = 'center', va = 'center', fontsize = fontsize)
    ax.annotate('High', (expression.shape[0]-(expression.shape[0]-high_patient_cutoff)/2, 0.5), ha = 'center', va = 'center', fontsize = fontsize)
    
    #add arrow indicating direction of expression
    if gene_name is not None:
        ax.arrow(low_patient_cutoff-5, 1.2, high_patient_cutoff-low_patient_cutoff + 10, 0, head_width = 0.2, color = 'black', head_length = 20, clip_on = False)
        ax.annotate(f'{gene_name} Expression', (expression.shape[0]/2, 1.4), ha = 'center', va = 'center', fontsize = fontsize, annotation_clip = False)

    
    ax.set_xlabel(xlabel, fontsize = fontsize)
    ax.set_yticks([])

