import os
import sys
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

import pandas as pd
import numpy as np
import anndata as ad
import scipy.stats as stats

import tacco as tc
import scanpy as sc


def _get_colors(colors, types):
    if colors is not None:
        colors = pd.Series(colors)
        if types.isin(colors.index).all():
            #colors = colors[colors.index.intersection(types)]
            types = colors.index
        else:
            raise Exception('Not all types %s are given colors with %s!' % (types, colors))
    else:
        colors = pd.Series(get_default_colors(len(types)), index=types)
    return colors, types


def subplots(
    n_x=1,
    n_y=1,
    axsize=(5,5),
    hspace=0.15,
    wspace=0.15,
    x_padding=None,
    y_padding=None,
    title=None,
    sharex='none',
    sharey='none',
    width_ratios=None,
    height_ratios=None,
    x_shifts=None,
    y_shifts=None,
):
    if x_padding is not None:
        wspace = x_padding / axsize[0]
    if y_padding is not None:
        hspace = y_padding / axsize[1]
    if width_ratios is None:
        effective_n_x = n_x
        effective_wspace = wspace
    else:
        effective_n_x = sum(width_ratios) / max(width_ratios)
        effective_wspace = wspace * n_x / effective_n_x
    if height_ratios is None:
        effective_n_y = n_y
        effective_hspace = hspace
    else:
        effective_n_y = sum(height_ratios) / max(height_ratios)
        effective_hspace = hspace * n_y / effective_n_y
    fig_height = axsize[1] * (effective_n_y + hspace * (n_y - 1))
    fig_width = axsize[0] * (effective_n_x + wspace * (n_x - 1))
    top = 1.0
    if title != None:
        title_space = 0.75
        fig_height += title_space
        top = 1 - title_space / fig_height
    fig, axs = plt.subplots(n_y,n_x,figsize=(fig_width,fig_height), squeeze=False, sharex=sharex, sharey=sharey, gridspec_kw={'wspace':effective_wspace,'hspace':effective_hspace,'left':0,'right':1,'top':top,'bottom':0,'width_ratios': width_ratios,'height_ratios': height_ratios})
    if title is not None:
        fig.suptitle(title, fontsize=16, y=1)
    if x_shifts is not None or y_shifts is not None:
        if x_shifts is None:
            x_shifts = [0.0] * n_x
        else:
            y_shifts = [0.0] * n_y
        for i_x in range(n_x):
            for i_y in range(n_y):
                [left,bottom,width,height] = axs[i_y,i_x].get_position().bounds
                axs[i_y,i_x].set_position([left+x_shifts[i_x],bottom+y_shifts[i_y],width,height])
    return fig, axs


def _anno_hist(mist, intervals, val, anno, norm=False):
    Nobs,Nanno = anno.shape
    assert(Nobs==len(mist))
    assert(Nobs==len(val))
    Nd = len(intervals)-1
    hist = np.zeros((Nd,))
    sums = np.zeros((Nd,Nanno))
    for i in range(Nobs):
        di = mist[i]
        _di = np.argmax(di <= intervals) - 1
        if di > 0:
            hist[_di] += val[i]
            sums[_di] += val[i] * anno[i]
    for d in range(Nd):
        sums[d] /= hist[d]
    if norm:
	sums1 = sums[~np.isnan(sums).any(axis=1)]
	sums1=sums1/(sums1.max(axis=0) + np.spacing(0))
	sums[:sums1.shape[0], :sums1.shape[1]] = sums1
    return sums


def plot_annotation_coordinate(
    adata,
    annotation_key,
    coordinate_key,
    group_key=None,
    reference_key=None,
    max_coordinate=None,
    delta_coordinate=None,
    axsize=None,
    colors=None,
    stacked=True,
    norm=False,
    custom_stick=True,
    max_dist=2000,
    layer_dists=[],
    layer_labels=[],
    lw=1,
    verbose=1,
):
    if group_key is None:
        group_adatas = {'':adata}
    else:
        group_adatas = {group:adata[df.index] for group,df in adata.obs.groupby(group_key) if len(df)>0 }
    if annotation_key in adata.obs:
        annotation = adata.obs[annotation_key]
        if hasattr(annotation, 'cat'):
            annotation = pd.get_dummies(annotation)
        else:
            annotation = pd.DataFrame(annotation)
    elif annotation_key in adata.obsm:
        annotation = adata.obsm[annotation_key].copy()
    else:
        raise ValueError(f'The `annotation_key` {annotation_key!r} is neither in `adata.obs` nor `adata.obsm`!')
    if pd.api.types.is_list_like(coordinate_key):
        if len(coordinate_key) == 2 and coordinate_key[0] in adata.obsm and coordinate_key[1] in adata.obsm[coordinate_key[0]]:
            coordinates = adata.obsm[coordinate_key[0]][coordinate_key[1]]
        else:
            raise ValueError(f'The `coordinate_key` {coordinate_key!r} is list/like, but not something of length 2 containing a `adata.obsm` key and a column name therein!')
        #coordinate_key = f'{coordinate_key[0]}:{coordinate_key[1]}'
        coordinate_key = f'{coordinate_key[1]}'
    elif coordinate_key in adata.obs:
        coordinates = adata.obs[coordinate_key]
    else:
        raise ValueError(f'The `coordinate_key` {coordinate_key!r} is not in `adata.obs` and not a valid specification for something in `adata.obsm`!')
    if max_coordinate is None or max_coordinate == np.inf:
        max_coordinate = coordinates.to_numpy().max()
    max_coordinate = float(max_coordinate)
    if delta_coordinate is None:
        delta_coordinate = max_coordinate / 100
    n_intervals = int(max_coordinate / delta_coordinate)
    max_coordinate = n_intervals * delta_coordinate
    intervals = np.arange(0,max_coordinate+delta_coordinate*0.5,delta_coordinate)
    midpoints = (intervals[1:] + intervals[:-1]) * 0.5
    if reference_key is None:
        reference_val = pd.Series(np.ones(len(adata.obs.index),dtype=float),index=adata.obs.index)
    elif reference_key in adata.obs:
        reference_val = adata.obs[reference_key]
    else:
        raise ValueError(f'The `reference_key` {reference_key!r} is not in `adata.obs`!')
    annotation_categories = annotation.columns
    colors,annotation_categories = _get_colors(colors, annotation_categories)
    if stacked:
        if axsize is None:
            axsize = (6,0.7)
        fig,axs = subplots(len(group_adatas),len(annotation.columns),axsize=axsize,sharex=True,sharey='row',y_padding=0,x_padding=0)
    else:
        if axsize is None:
            axsize = (6,4)
        fig,axs = subplots(len(group_adatas),1,axsize=axsize,sharex=True,sharey='row',y_padding=0,x_padding=0)
    for gi,(group,group_adata) in enumerate(group_adatas.items()):
        group_annotation = annotation.loc[group_adata.obs.index].to_numpy()
        group_coordinates = coordinates.loc[group_adata.obs.index].to_numpy()
        group_reference_val = reference_val.loc[group_adata.obs.index].to_numpy()
        no_nan = ~np.isnan(group_annotation).any(axis=1)
	do_norm = norm
        anno = _anno_hist(group_coordinates[no_nan], intervals, group_reference_val[no_nan], group_annotation[no_nan], norm=do_norm)
        anno = pd.DataFrame(anno,index=midpoints,columns=annotation.columns)
        if stacked:
            for i,c in enumerate(colors.index):
                axs[i,gi].plot(anno[c],c=colors[c],label=c, linewidth=lw)
		for ld in layer_dists:
		    axs[i,gi].axvline(ld, color='lightgray', linestyle='--')
        else:
            for i,c in enumerate(colors.index):
                axs[0,gi].plot(anno[c],c=colors[c],label=c, linewidth=lw)
		for ld in layer_dists:
                    axs[0,gi].axvline(ld, color='lightgray', linestyle='--')
        group_string = f' in {group}' if group != '' else ''
        axs[0,gi].set_title(f'{annotation_key} VS distance from {coordinate_key}{group_string}')
	xpos = np.array([0])
	for i in range(len(layer_dists)):
            if i+1 == len(layer_dists):
	        break
	    xpos=np.append(xpos, layer_dists[[i, i+1]].sum()/2)
	xpos=np.append(xpos, np.array([layer_dists[-1],max_dist]).sum()/2)
	if custom_stick:
            axs[0,gi].set_xticks(xpos, layer_labels)
	axs[0,gi].set_xlabel('Distance')
    axs[0,-1].legend(handles=[mpatches.Patch(color=color, label=ind) for (ind, color) in colors.items() ],
        bbox_to_anchor=(1, 1), loc='upper left', ncol=1)
    return fig


def plot_ci_bootstrap(xs, ys, resid, nboot=500, ax=None):
    if ax is None:
        ax = plt.gca()
    bootindex = sp.random.randint
    for _ in range(nboot):
        resamp_resid = resid[bootindex(0, len(resid) - 1, len(resid))]
        # Make coeffs of for polys
        pc = sp.polyfit(xs, ys + resamp_resid, 1)                   
        # Plot bootstrap cluster
        ax.plot(xs, sp.polyval(pc, xs), "b-", linewidth=2, alpha=3.0 / float(nboot))
    return ax


def plot_ci_manual(t, s_err, n, x, x2, y2, ax=None):
    if ax is None:
        ax = plt.gca()
    ci = t * s_err * np.sqrt(1/n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))
    ax.fill_between(x2, y2 + ci, y2 - ci, color="gray", alpha=0.5)
    return ax


def equation(a, b):
    """Return a 1D polynomial."""
    return np.polyval(a, b) 


def plot_line(x, y, deg=1, nam='test'):
    p, cov = np.polyfit(x, y, deg, cov=True)
    y_model = equation(p, x)
    # Statistics
    n = val.size                                           # number of observations
    m = p.size                                                 # number of parameters
    dof = n - m                                                # degrees of freedom
    t = stats.t.ppf(0.975, n - m)                              # t-statistic; used for CI and PI bands
    # Estimates of Error in Data/Model
    resid = y - y_model                                        # residuals; diff. actual data from predicted values
    chi2 = np.sum((resid / y_model)**2)                        # chi-squared; estimates error in data
    chi2_red = chi2 / dof                                      # reduced chi-squared; measures goodness of fit
    s_err = np.sqrt(np.sum(resid**2) / dof)                    # standard deviation of the error
    # Plotting --------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))
    # Data
    ax.scatter(x, y, s=5, alpha=0.9, color='black')
    # Fit
    ax.plot(x, y_model, "-", color="0.1", linewidth=3, alpha=0.6, label="Fit")  
    x2 = np.linspace(np.min(x), np.max(x), 100)
    y2 = equation(p, x2)
    # Confidence Interval (select one)
    plot_ci_manual(t, s_err, n, x, x2, y2, ax=ax)
    #plot_ci_bootstrap(x, y, resid, ax=ax)
    # Prediction Interval
    #pi = t * s_err * np.sqrt(1 + 1/n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))   
    #ax.fill_between(x2, y2 + pi, y2 - pi, color="None", linestyle="--")
    #ax.plot(x2, y2 - pi, "--", color="0.5", label="95% Prediction Limits")
    #ax.plot(x2, y2 + pi, "--", color="0.5")
    ax.set_title(nam)
    #ax.set_xlabel(r'$Layer \hspace{1} 2 \hspace{5} \longrightarrow \hspace{5} WM$', color='black')
    ax.set_xlabel('Layers')
    ax.set_xticks([])
    ax.set_ylabel("Normalized gene expression")
    ax.set_xlim((0,2300))
    ax.set_ylim(bottom=0)
    ax.axhline(300, color='lightgray', linestyle='--') # 300 L1
    ax.axhline(850, color='lightgray', linestyle='--') # 550 L2-3
    ax.axhline(1100, color='lightgray', linestyle='--') # 250 L4
    ax.axhline(1700, color='lightgray', linestyle='--') # 600 L5
    ax.axhline(2000, color='lightgray', linestyle='--') # 300 L6
    ax.axhline(2100, color='lightgray', linestyle='--') # 100 L6b
    fig.savefig(f"supp_fig_7c_{nam}_line_fit.tiff", bbox_inches='tight', dpi=300)
    plt.close()
    plt.clf()



# plotting options
highres = True
default_dpi = 100.0
if highres:
    matplotlib.rcParams['figure.dpi'] = 648.0
    hr_ext = '_hd'
else:
    matplotlib.rcParams['figure.dpi'] = default_dpi
    hr_ext = ''

axsize = np.array([6,6])

ref_file = "/stanley/levin_dr/kwanho/projects/vahbiz/analysis_20211210/new_data_analysis/seurat/OL_ref_for_slideseq/ref_OL_types.h5ad"  # use our scRNA of oligos both reps : OPC COP NFOL MFOL MOL (use MOL subtypes)
puck_file = "../puck_final_regions.h5ad"  # oligo defined as oligo

ref = ad.read(ref_file)
puck = ad.read(puck_file)

# sanity check ref
#tc.pl.scatter(reference, keys='subclass_label', position_key='X_umap', joint=True, point_size=5, axsize=axsize, noticks=True, axes_labels=['UMAP 0','UMAP 1']);
#plt.savefig("ref_umap.png", bbox_inches='tight')

# subset reference class represented by fewer than 10 cells
#reference.var.index = ref.var.features
#keep_subclass = reference.obs.subclass_label.value_counts().index[reference.obs.subclass_label.value_counts()>10]
#ref = reference[reference.obs['subclass_label'].isin(keep_subclass)]

#puck.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'feature_names'}, inplace=True)
#puck.var.drop(columns='_index', inplace=True)
ref.var.set_index('features', inplace=True)
ref.var.drop(columns='_index', inplace=True)
ref.obs['final_name'] = ref.obs.final_name.astype('category')
#in order: OPC   COP  NFOL  MFOL MOL_A MOL_B MOL_C
ref.obs['final_name'] = ref.obs.final_name.cat.rename_categories({0: 'OPC', 1: 'COP', 2: 'NFOL', 3: 'MFOL', 4: 'MOL_A', 5: 'MOL_B', 6: 'MOL_C'})

# subset puck to those spots with >0 oligo proportions

# annotate puck
anno_nam = 'anno_OL'
tc.tl.annotate(puck,ref,annotation_key='final_name',result_key=anno_nam, max_annotation=5, min_log2foldchange=1)
puck.write_h5ad(f"puck_{anno_nam}.h5ad")

#genes=[]  # OL genes
#sc.pl.spatial(puck, color=genes, use_raw=False, spot_size=10, size=5, color_map="viridis_r")
#plt.savefig("oligo_gene_exp.png", dpi=300, bbox_inches='tight')

#puck.write_h5ad("puck_oligo_cortex.h5ad")

# distribute Oligo annotation into OL subtype annotation
#puck.obs['Oligo_frac'] = puck.obsm['anno_OT_max5']['Oligo']
puck.obsm['anno_OL'][puck.obsm['anno_OL'] < 0.1] = 0
puck.obsm['anno_OL'] = puck.obsm['anno_OL'].div(puck.obsm['anno_OL'].sum(axis=1), axis=0)
puck.obsm['anno_OL'] = puck.obsm['anno_OL'].fillna(0)
puck.obsm['OL_subtypes'] = puck.obsm['anno_OL'].mul(puck.obsm['clean_anno']['Oligo'], axis=0)
puck.obsm['final_name'] = puck.obsm['clean_anno'].drop(columns=['Oligo'])
puck.obsm['final_name'] = pd.concat([puck.obsm['final_name'], puck.obsm['OL_subtypes']], axis=1)
puck.obsm['final_name']['MOL'] = puck.obsm['final_name']['MOL_A'] + puck.obsm['final_name']['MOL_B'] + puck.obsm['final_name']['MOL_C']

# visualize cell types in slide-seq
tc.pl.scatter(puck, keys='OL_subtypes', position_key=['x','y'], joint=True, point_size=8, axsize=axsize, noticks=False, axes_labels=['X','Y']);
plt.savefig(f"slideseq_OL_subtypes.png", bbox_inches='tight')


puck.obsm['cell_type'] = puck.obsm['final_name'][['L2/3 IT CTX-1', 'L4/5 IT CTX', 'L5 IT CTX', 'L5 NP CTX','L5 PT CTX','L6 IT CTX', 'L6b CTX','OPC', 'COP', 'NFOL', 'MFOL','MOL','MOL_A', 'MOL_B', 'MOL_C']]


neu_type_simple = ['L2/3 IT CTX-1', 'L4/5 IT CTX', 'L5 IT CTX','L6 IT CTX', 'L6b CTX']
neu_type = ['L2/3 IT CTX-1', 'L4/5 IT CTX', 'L5 IT CTX', 'L5 NP CTX','L5 PT CTX','L6 CT CTX','L6b CTX']
ol_type = ['OPC','COP','NFOL','MFOL','MOL']
mols = ['MOL_A','MOL_B','MOL_C']
my_cols = {'L2/3 IT CTX-1':"#88B2D5", 'L4/5 IT CTX':"#A6A2E0", 'L5 IT CTX':"#ffc200", 'L5 NP CTX':"#DD3497",  'L5 PT CTX':"#FA9FB5", 'L6 CT CTX':"#D9F0A3", 'L6b CTX':"#8DD3C7", 'OPC':"#A6CEE3", 'COP':"#b452cd", 'NFOL':"#028F1E", 'MFOL':"#FF474C",'MOL':"#FFB16D", 'MOL_A': "#FFA500", 'MOL_B': "#FF5800", 'MOL_C': "#883000"}

num_plots = len(mols)
#nc = int(np.ceil(np.sqrt(num_plots)))
#nr = int(np.ceil(num_plots / nc))
nc=3
nr=1
fig,axs = tc.pl.subplots(n_x=nc, n_y=nr, x_padding=1.5)
tc.pl.scatter(puck, keys='cell_type', show_only=mols, position_key=['x','y'], joint=False, point_size=30, axsize=axsize, noticks=True, axes_labels=['X','Y'], ax=axs.reshape(nc*nr,1)[:num_plots,:], cmap='magma_r', legend=True)
plt.savefig(f"supp_fig_3i_slideseq_MOL_subtype_individual.tiff", bbox_inches='tight', dpi=300)

num_plots = len(ol_type)
nc=5
nr=1
fig,axs = tc.pl.subplots(n_x=nc, n_y=nr, x_padding=1.5)
tc.pl.scatter(puck, keys='cell_type', show_only=ol_type, position_key=['x','y'], joint=False, point_size=30, axsize=axsize, noticks=True, axes_labels=['x','y'], ax=axs.reshape(nc*nr,1)[:num_plots,:], cmap='magma_r', legend=True)
plt.savefig(f"supp_fig_2f_slideseq_oligo_type_individual.tiff", bbox_inches='tight', dpi=300)

num_plots = len(neu_type)
nc=7
nr=1
fig,axs = tc.pl.subplots(n_x=nc, n_y=nr, x_padding=1.5)
tc.pl.scatter(puck, keys='cell_type', show_only=neu_type, position_key=['x','y'], joint=False, point_size=30, axsize=axsize, noticks=True, axes_labels=['x','y'], ax=axs.reshape(nc*nr,1)[:num_plots,:], cmap='magma_r', legend=True)
plt.savefig(f"supp_fig_2d_slideseq_neuron_type_individual_with_L6CT.tiff", bbox_inches='tight', dpi=300)



puck.obsm['neu_type'] = puck.obsm['final_name'][neu_type]
puck.obsm['neu_type_simple'] = puck.obsm['final_name'][neu_type_simple]
puck.obsm['ol_type'] = puck.obsm['final_name'][ol_type]
puck.obsm['mol_subtype'] = puck.obsm['final_name'][mols]


# compute distance
#puck.obsm['anno_MOL'] = puck.obsm['OL_subtypes'][mols]
region_key = 'final_regions'
md = 2000
delta = 100
delta2 = 200
tc.tl.annotation_coordinate(puck,annotation_key=region_key,result_key='OL_dist', max_distance=md, delta_distance=delta);

region_colors = tc.pl.get_default_colors(puck.obs[region_key].cat.categories.values, offset=17)

# Plot of cell densities at various distances away from Pia
# all in one
#fig = plot_annotation_coordinate(puck,annotation_key='ol_type',coordinate_key=('OL_dist','Pia'),axsize=(4,3), max_coordinate=md, delta_coordinate=delta2, stacked=False, colors=dict((k, my_cols[k]) for k in ol_type), lw=3)
#fig.savefig(f"main_fig_1e_line_graph_OL_densities.tiff", bbox_inches='tight', dpi=300)
#
#fig = plot_annotation_coordinate(puck,annotation_key='mol_subtype',coordinate_key=('OL_dist','Pia'),axsize=(4,3), max_coordinate=md, delta_coordinate=delta2, stacked=False, colors=dict((k, my_cols[k]) for k in mols), lw=3)
#fig.savefig(f"main_fig_1g_line_graph_MOL_densities.tiff", bbox_inches='tight', dpi=300)
#
#fig = plot_annotation_coordinate(puck,annotation_key='neu_type_simple',coordinate_key=('OL_dist','Pia'),axsize=(4,3), max_coordinate=md, delta_coordinate=delta2, stacked=False, colors=dict((k, my_cols[k]) for k in neu_type_simple), lw=3)
#fig.savefig(f"main_fig_1d_line_graph_Neu_densities.tiff", bbox_inches='tight', dpi=300)
#
## Each cell type plotted separately
#fig = plot_annotation_coordinate(puck,annotation_key='ol_type',coordinate_key=('OL_dist','Pia'),axsize=(4,0.7), max_coordinate=md, delta_coordinate=delta, stacked=True, colors=dict((k, my_cols[k]) for k in ol_type), lw=3)
#fig.savefig(f"main_fig_1e_line_graph_stacked_OL_density_along_layers.tiff", bbox_inches='tight', dpi=300)
#
#fig = plot_annotation_coordinate(puck,annotation_key='mol_subtype',coordinate_key=('OL_dist','Pia'),axsize=(4,0.7), max_coordinate=md, delta_coordinate=delta, stacked=True, colors=dict((k, my_cols[k]) for k in mols), lw=3)
#fig.savefig(f"main_fig_1g_line_graph_stacked_MOL_density_along_layers.tiff", bbox_inches='tight', dpi=300)
#
#fig = plot_annotation_coordinate(puck,annotation_key='neu_type',coordinate_key=('OL_dist','Pia'),axsize=(4,0.7), max_coordinate=md, delta_coordinate=delta, stacked=True, colors=dict((k, my_cols[k]) for k in neu_type), lw=3)
#fig.savefig(f"main_fig_1d_line_graph_stacked_Neu_density_along_layers.tiff", bbox_inches='tight', dpi=300)


# distance without WM
wm_cells = puck.obs.final_regions.index.values[puck.obs.final_regions=='WM']
rest_cells = [x for x in puck.obs.index.values if x not in wm_cells]
spuck = puck[rest_cells]

md=2000
delta = 100
delta2 = 200
tc.tl.annotation_coordinate(spuck,annotation_key=region_key,result_key='OL_dist2', max_distance=md, delta_distance=delta);

max_dist = 2000
layer_dists = np.array([30, 550,800,1400,1700])+50
layer_labels=['Pia', 'L2-3','L4','L5','L6','L6b']

# all in one
fig = plot_annotation_coordinate(spuck,annotation_key='ol_type',coordinate_key=('OL_dist2','Pia'),axsize=(4,3), max_coordinate=md, delta_coordinate=delta2, stacked=False, colors=dict((k, my_cols[k]) for k in ol_type), lw=3, layer_dists=layer_dists, layer_labels=layer_labels,max_dist=max_dist)
fig.savefig(f"main_fig_1e_line_graph_OL_densities_no_WM.tiff", bbox_inches='tight', dpi=300)

#fig = plot_annotation_coordinate(spuck,annotation_key='mol_subtype',coordinate_key=('OL_dist2','Pia'),axsize=(4,3), max_coordinate=md, delta_coordinate=delta2, stacked=False, colors=dict((k, my_cols[k]) for k in mols), lw=3, layer_dists=layer_dists, layer_labels=layer_labels,max_dist=max_dist)
#fig.savefig(f"main_fig_1g_line_graph_MOL_densities_no_WM.tiff", bbox_inches='tight', dpi=300)
#
#fig = plot_annotation_coordinate(spuck,annotation_key='neu_type_simple',coordinate_key=('OL_dist2','Pia'),axsize=(4,3), max_coordinate=md, delta_coordinate=delta2, stacked=False, colors=dict((k, my_cols[k]) for k in neu_type_simple), lw=3, layer_dists=layer_dists, layer_labels=layer_labels,max_dist=max_dist)
#fig.savefig(f"main_fig_1d_line_graph_Neu_densities_no_WM.tiff", bbox_inches='tight', dpi=300)

# all in one - normalized cell densities by the maximum in each cell type
#fig = plot_annotation_coordinate(spuck,annotation_key='ol_type',coordinate_key=('OL_dist','Pia'),axsize=(4,3), max_coordinate=md, delta_coordinate=delta2, stacked=False, colors=dict((k, my_cols[k]) for k in ol_type), lw=3, norm=True, layer_dists=layer_dists, layer_labels=layer_labels,max_dist=max_dist)
#fig.savefig(f"main_fig_1e_line_graph_OL_normalized_density_no_WM.tiff", bbox_inches='tight', dpi=300)

#fig = plot_annotation_coordinate(spuck,annotation_key='mol_subtype',coordinate_key=('OL_dist','Pia'),axsize=(4,3), max_coordinate=md, delta_coordinate=delta2, stacked=False, colors=dict((k, my_cols[k]) for k in mols), lw=3, norm=True, layer_dists=layer_dists, layer_labels=layer_labels,max_dist=max_dist)
#fig.savefig(f"main_fig_1g_line_graph_MOL_normalized_density_no_WM.tiff", bbox_inches='tight', dpi=300)

#fig = plot_annotation_coordinate(spuck,annotation_key='neu_type_simple',coordinate_key=('OL_dist','Pia'),axsize=(4,3), max_coordinate=md, delta_coordinate=delta2, stacked=False, colors=dict((k, my_cols[k]) for k in neu_type_simple), lw=3, norm=True, layer_dists=layer_dists, layer_labels=layer_labels,max_dist=max_dist)
#fig.savefig(f"main_fig_1d_line_graph_Neu_normalized_density_no_WM.tiff", bbox_inches='tight', dpi=300)

# Each cell type plotted separately
#fig = plot_annotation_coordinate(spuck,annotation_key='ol_type',coordinate_key=('OL_dist2','Pia'),axsize=(4,0.7), max_coordinate=md, delta_coordinate=delta2, stacked=True, colors=dict((k, my_cols[k]) for k in ol_type), lw=3, layer_dists=layer_dists, layer_labels=layer_labels,max_dist=max_dist)
#fig.savefig(f"main_fig_1e_line_graph_stacked_OL_density_along_layers_no_WM.tiff", bbox_inches='tight', dpi=300)

fig = plot_annotation_coordinate(spuck,annotation_key='mol_subtype',coordinate_key=('OL_dist2','Pia'),axsize=(4,0.7), max_coordinate=md, delta_coordinate=delta2, stacked=True, colors=dict((k, my_cols[k]) for k in mols), lw=3, layer_dists=layer_dists, layer_labels=layer_labels,max_dist=max_dist)
fig.savefig(f"main_fig_1g_line_graph_stacked_MOL_density_along_layers_no_WM.tiff", bbox_inches='tight', dpi=300)

fig = plot_annotation_coordinate(spuck,annotation_key='neu_type',coordinate_key=('OL_dist2','Pia'),axsize=(4,0.7), max_coordinate=md, delta_coordinate=delta2, stacked=True, colors=dict((k, my_cols[k]) for k in neu_type), lw=3, layer_dists=layer_dists, layer_labels=layer_labels,max_dist=max_dist)
fig.savefig(f"main_fig_1d_line_graph_stacked_Neu_density_along_layers_no_WM_L6CT.tiff", bbox_inches='tight', dpi=300)




# _supp fig 7c
genes = ['Gas6','Fgf18','Ncam1','Sema5a','Rspo3']
axes = sc.pl.spatial(puck, color=genes, use_raw=False, spot_size=10, size=5, color_map="magma_r", ncols=5, return_fig=True)
for ax in axes:
    ax.axhline(300, color='lightgray', linestyle='--') # 300 L1
    ax.axhline(850, color='lightgray', linestyle='--') # 550 L2-3
    ax.axhline(1100, color='lightgray', linestyle='--') # 250 L4
    ax.axhline(1700, color='lightgray', linestyle='--') # 600 L5
    ax.axhline(2000, color='lightgray', linestyle='--') # 300 L6
    ax.axhline(2100, color='lightgray', linestyle='--') # 300 L6

plt.savefig("supp_fig_7c_gene_exp_in_slide_seq.tiff", dpi=300, bbox_inches='tight')

sc.pl.dotplot(puck, genes, groupby='final_regions')
plt.savefig("dotplot_genes.png", bbox_inches='tight', dpi=300)

sc.pl.stacked_violin(puck, genes, groupby='final_regions', row_palette=region_colors.values())
plt.savefig("violin_genes.png", bbox_inches='tight', dpi=300)


for g in genes:
    print(g)
    puck.obs['temp'] = puck[:, g].X.todense() > 0
    val = np.squeeze(np.asarray(puck[puck.obs['temp']==True, puck.var_names.str.match(g)].X.todense()))
    spatial_y = np.array(puck[puck.obs['temp']==True].obsm['spatial'][:,1])
    #spatial_y = np.array(puck[puck.obs['temp']==True].obs.y)
    df = pd.DataFrame({'x':spatial_y, 'y':val})
    plt.figure()
    plt.hist(df.x.tolist(), bins=[0,300,850,1100,1700,2000,2100,2300], alpha=0.5, color='steelblue', edgecolor='gray')
    plt.xticks([150,575,975,1400,1850,2050,2200], ['Pia','L2-3','L4','L5','L6','L6b','WM'])
    plt.title(g)
    plt.ylabel('Number of cells')
    plt.xlabel('Spatial location binned by layers')
    plt.savefig(f"hist_nCells_expressing_{g}.png", bbox_inches='tight', dpi=300)
    df = df.drop_duplicates(subset=['x'])
    df = df.sort_values(by=['x'])
    plot_line(df.x.tolist(), df.y.tolist(), deg=2, nam=g)



puck.write("puck_final.h5ad")



# Not used

cts = ['Astro', 'Endo', 'L2/3 IT CTX-1','L4/5 IT CTX', 'L5 IT CTX', 'L5 NP CTX', 'L5 PT CTX', 'L6 CT CTX','L6 IT CTX', 'L6b CTX', 'Lamp5', 'Micro-PVM', 'Pvalb', 'SMC-Peri','Sncg', 'Sst', 'Sst Chodl', 'Vip']
# Analyse co-occurrence and neighbourhips
puck.obsm['final_name'].drop(columns='MOL', inplace=True)
tc.tl.co_occurrence(puck, 'final_name', result_key='type-type',delta_distance=20,max_distance=1000,sparse=False,n_permutation=10)
#, ax=axs.reshape(nc*nr,1)[:num_plots,:])
tc.pl.co_occurrence(puck, 'type-type', show_only=cts, show_only_center=mols, log_base=2, wspace=0.25) 
plt.savefig("cell_type_cooccurrence1.png", dpi=300, bbox_inches='tight')
tc.pl.co_occurrence(puck, 'type-type', show_only=mols, show_only_center=cts, log_base=2, wspace=0.25, colors=dict((k, my_cols[k]) for k in mols))
plt.savefig("cell_type_cooccurrence2.png", dpi=300, bbox_inches='tight')







