import os
import sys
import matplotlib
from matplotlib import pyplot as plt

import pandas as pd
import numpy as np
import anndata as ad

import tacco as tc
import scanpy as sc


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

#labels_colors = pd.Series({'Epi': (0.00784313725490196, 0.24313725490196078, 1.0), 'B': (0.10196078431372549, 0.788235294117647, 0.2196078431372549), 'TNK': (1.0, 0.48627450980392156, 0.0), 'Mono': (0.5490196078431373, 0.03137254901960784, 0.0), 'Mac': (0.9098039215686274, 0.0, 0.043137254901960784), 'Gran': (0.34901960784313724, 0.11764705882352941, 0.44313725490196076), 'Mast': (0.23529411764705882, 0.23529411764705882, 0.23529411764705882), 'Endo': (0.8549019607843137, 0.5450980392156862, 0.7647058823529411), 'Fibro': (0.6235294117647059, 0.2823529411764706, 0.0)})
#region_colors = tc.pl.get_default_colors([f'region_{i}' for i in range(4)], offset=17)
#split_names = np.array([f'sub_{i}' for i in range(4)])
#split_colors = tc.pl.get_default_colors(split_names, offset=12)


ref_file = "../data/ref_counts_allen_p56.h5ad"
puck_file = "data/slideseq.h5ad"

reference = ad.read(ref_file)
puck = ad.read(puck_file)

# sanity check ref
#tc.pl.scatter(reference, keys='subclass_label', position_key='X_umap', joint=True, point_size=5, axsize=axsize, noticks=True, axes_labels=['UMAP 0','UMAP 1']);
#plt.savefig("ref_umap.png", bbox_inches='tight')

# subset reference class represented by fewer than 10 cells
reference.var.index = reference.var.features
keep_subclass = reference.obs.subclass_label.value_counts().index[reference.obs.subclass_label.value_counts()>10]
ref = reference[reference.obs['subclass_label'].isin(keep_subclass)]

# subset puck
puck = puck[tc.sum(puck.X,axis=1)>=50].copy() # restrict downstream analysis to "good" beads

# visualize puck
puck.obs['total_counts'] = tc.sum(puck.X,axis=1)
puck.obs['log10_counts'] = np.log10(1+puck.obs['total_counts'])
tc.pl.scatter(puck, 'log10_counts', cmap='viridis', cmap_vmin_vmax=[1,3]);
plt.savefig("slideseq_puck_filt.png", bbox_inches='tight')

#cluster2type = reference.obs[['class_label', 'subclass_label']].drop_duplicates().groupby('subclass_label')['class_label'].agg(lambda x: list(x.to_numpy()))


#puck.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'feature_names'}, inplace=True)
#puck.var.drop(columns='_index', inplace=True)
#ref.var.set_index('features', inplace=True)
#ref.var.drop(columns='_index', inplace=True)

#common_genes = [x for x in ref.var.index if x in puck.var.index]
#gene_idx = [puck.var_names.get_loc(gene) for gene in common_genes]
#puck[:,gene_idx].X.sum(axis=1)


# subset to cortex
x_min = 250
x_max = 1500
y_min = 4500
y_max = 6900

df = puck.obs
my_obs = df.index[(x_min < df['x']) & (df['x'] < x_max) & (y_min < df['y']) & (df['y'] < y_max)]
spuck = puck[my_obs]

del puck




# annotate puck
anno_nam = 'anno_OT_max5'
tc.tl.annotate(spuck,ref,annotation_key='subclass_label',result_key=anno_nam, max_annotation=5, min_log2foldchange=1)
spuck.write_h5ad(f"puck_{anno_nam}.h5ad")

#spuck = ad.read(f"puck_{anno_nam}.h5ad")

# acronyms
# PVM, perivascular macrophage; SMC, smooth muscle cell; VLMC, vascular lepotomeningeal cell; 
# IT, intratelencephalic; PT, pyramidal tract; NP, near-projecting; CT, corticothalamic

# visualize cell types in slide-seq
tc.pl.scatter(spuck, keys=anno_nam, position_key=['x','y'], joint=True, point_size=8, axsize=axsize, noticks=False, axes_labels=['X','Y']);
plt.savefig(f"slideseq_{anno_nam}.png", bbox_inches='tight')

cts = ['L2/3 IT CTX-1','L2/3 IT CTX-2','L4/5 IT CTX','L5 IT CTX','L5 NP CTX','L5 PT CTX','L6 IT CTX','L6 CT CTX','L6b CTX','Oligo']
num_plots = len(cts)
nc = int(np.ceil(np.sqrt(num_plots)))
nr = int(np.ceil(num_plots / nc))
fig,axs = tc.pl.subplots(n_x=nr, n_y=nc, x_padding=1.5)
tc.pl.scatter(spuck, keys=anno_nam, show_only=cts, position_key=['x','y'], joint=False, point_size=8, axsize=axsize, noticks=True, axes_labels=['X','Y'], ax=axs.reshape(nc*nr,1)[:num_plots,:])
plt.savefig(f"slideseq_{anno_nam}_select_CT.png", bbox_inches='tight')



# plot gene exp (y becomes flipped)
df = spuck.obs.loc[:,['x','y']]
df.rename(columns={'x':'x_coord','y':'y_coord'}, inplace=True)
spuck.obsm['spatial'] = df.values
spuck.obsm['spatial'][:,1] = max(spuck.obsm['spatial'][:,1]) - spuck.obsm['spatial'][:,1]

genes=['Neurod2','Slc17a7','Cux1','Cux2','Satb2','Rorb','Fezf2','Bcl11b','Tcerg1l','Tle4','Olig1','Olig2','Plp1','Pdgfra','Mog','Tmem119','Sall1','Aldh1l1','Flt1','Vtn']
sc.pl.spatial(spuck, color=genes, use_raw=False, spot_size=10, size=5, color_map="viridis_r")
plt.savefig("gene_exp.png", dpi=300, bbox_inches='tight')



# clean up annotation
spuck.obsm['clean_anno'] = spuck.obsm[anno_nam].copy()
spuck.obsm['clean_anno'][spuck.obsm['clean_anno'] < 0.1] = 0
spuck.obsm['clean_anno']= spuck.obsm['clean_anno'].div(spuck.obsm['clean_anno'].sum(axis=1), axis=0)
spuck.obsm['clean_anno']=spuck.obsm['clean_anno'].fillna(0)

anno_nam = 'clean_anno'
# visualize cell types in slide-seq
tc.pl.scatter(spuck, keys=anno_nam, position_key=['x','y'], joint=True, point_size=8, axsize=axsize, noticks=False, axes_labels=['X','Y']);
plt.savefig(f"slideseq_{anno_nam}.png", bbox_inches='tight')

#cts = ['L2/3 IT CTX-1','L2/3 IT CTX-2','L4/5 IT CTX','L5 IT CTX','L5 NP CTX','L5 PT CTX','L6 IT CTX','L6 CT CTX','L6b CTX','Oligo']
#num_plots = len(cts)
#nc = int(np.ceil(np.sqrt(num_plots)))
#nr = int(np.ceil(num_plots / nc))
#fig,axs = tc.pl.subplots(n_x=nr, n_y=nc, x_padding=1.5)
#tc.pl.scatter(spuck, keys=anno_nam, show_only=cts, position_key=['x','y'], joint=False, point_size=8, axsize=axsize, noticks=True, axes_labels=['X','Y'], ax=axs.reshape(nc*nr,1)[:num_plots,:])
#plt.savefig(f"slideseq_{anno_nam}_select_CT.png", bbox_inches='tight')




#spuck.obs.y = max(spuck.obs.y) - spuck.obs.y

spuck.write_h5ad("puck_cortex.h5ad")

#tc.pl.scatter(spuck, keys=anno_nam, position_key=['x','y'], joint=True, point_size=8, axsize=axsize, noticks=False, axes_labels=['X','Y']);
#plt.savefig(f"slideseq_merged_anno.png", bbox_inches='tight')

#cts = spuck.obsm['anno_OT_max5'].columns.values
#num_plots = len(cts)
#nc = int(np.ceil(np.sqrt(num_plots)))
#nr = int(np.ceil(num_plots / nc))
#fig,axs = tc.pl.subplots(n_x=nr, n_y=nc, x_padding=1.5)
#tc.pl.scatter(spuck, keys=anno_nam, show_only=cts, position_key=['x','y'], joint=False, point_size=10, axsize=axsize, noticks=True, axes_labels=['x','y'], ax=axs.reshape(nc*nr,1)[:num_plots,:])
#plt.savefig(f"slideseq_individual_anno.png", bbox_inches='tight')

# layer associated cts only
cts = ['L2/3 IT CTX-1','L4/5 IT CTX','L5 IT CTX','L6 IT CTX','L6b CTX','Oligo']
num_plots = len(cts)
nc = int(np.ceil(np.sqrt(num_plots)))
nr = int(np.ceil(num_plots / nc))
fig,axs = tc.pl.subplots(n_x=nr, n_y=nc, x_padding=1.5)
tc.pl.scatter(spuck, keys=anno_nam, show_only=cts, position_key=['x','y'], joint=False, point_size=15, axsize=axsize, noticks=False, axes_labels=['x','y'], ax=axs.reshape(nc*nr,1)[:num_plots,:])
plt.savefig(f"slideseq_anno_layer_associated_cell_types.png", bbox_inches='tight')

# other cell types
#cts = ['Oligo','Micro-PVM','Astro']
#num_plots = len(cts)
#nc = int(np.ceil(np.sqrt(num_plots)))
#nr = int(np.ceil(num_plots / nc))
#fig,axs = tc.pl.subplots(n_x=nr, n_y=nc, x_padding=1.5)
#tc.pl.scatter(spuck, keys=anno_nam, show_only=cts, position_key=['x','y'], joint=False, point_size=10, axsize=axsize, noticks=True, axes_labels=['x','y'], ax=axs.reshape(nc*nr,1)[:num_plots,:])
#plt.savefig(f"slideseq_anno_glial.png", bbox_inches='tight')


# simplify cell type deconvolution to identify layers
#cts = ['L2/3 IT CTX-1','L5 IT CTX','L6 IT CTX','L6b CTX','Oligo']

#spuck.obsm['anno_layers'] = spuck.obsm[anno_nam].copy()
#spuck.obsm['anno_layers'].drop(columns=['Astro', 'Car3', 'Endo', 'L2/3 IT CTX-2','Lamp5', 'Micro-PVM', 'Pvalb','SMC-Peri', 'Sncg', 'Sst', 'Sst Chodl', 'Vip'], inplace=True)
#spuck.obsm['anno_layers'] = spuck.obsm['anno_layers'][cts]
#spuck.obsm['anno_layers']= spuck.obsm['anno_layers'].div(spuck.obsm['anno_layers'].sum(axis=1), axis=0)
#spuck.obsm['anno_layers']=spuck.obsm['anno_layers'].fillna(0)

#anno_layer = 'anno_layers'
# Find Regions
pw=1.5
res=0.7
#pw=1.4
#res=1
#region_key = f'regions_pw{pw}_res{res}'
region_key = f'regions2'
#tc.tl.find_regions(spuck,key_added=region_key, position_weight=pw, resolution=res, annotation_key=anno_layer);
tc.tl.find_regions(spuck,key_added=region_key, position_weight=pw, resolution=res, annotation_key=None);
spuck.obs[region_key] = spuck.obs[region_key].map(lambda x: f'region_{x}').astype('category')
print(spuck.obs[region_key].value_counts())
region_colors = tc.pl.get_default_colors([f'region_{i}' for i in range(len(spuck.obs[region_key].value_counts().index))], offset=17)
tc.pl.scatter(spuck, region_key, colors=region_colors, axes_labels=['X','Y'], noticks=True)
plt.savefig(f"slideseq_{region_key}_.png", bbox_inches='tight')


# Using variable genes, tacco finds: pia, L2-4, L5, L6 and WM
# Using the deconvoluted cell types and kwown gene expression, manually identify L4 and L6b
spuck.obs[region_key] = spuck.obs[region_key].cat.add_categories(['region_5','region_6'])
dat = spuck.obs[region_key]
l6_cells = dat[dat=='region_2'].index.values
l6b_wm = spuck.obs.y[spuck.obs.y<4800].index.values
l6b_cells = [x for x in l6_cells if x in l6b_wm]
spuck.obs[region_key][l6b_cells] = 'region_5'

print(spuck.obs[region_key].value_counts())

#def rotate_matrix (x, y, angle, x_shift=0, y_shift=0, units="DEGREES"):
#    """
#    Rotates a point in the xy-plane counterclockwise through an angle about the origin
#    https://en.wikipedia.org/wiki/Rotation_matrix
#    :param x: x coordinate
#    :param y: y coordinate
#    :param x_shift: x-axis shift from origin (0, 0)
#    :param y_shift: y-axis shift from origin (0, 0)
#    :param angle: The rotation angle in degrees
#    :param units: DEGREES (default) or RADIANS
#    :return: Tuple of rotated x and y
#    """
#    # Shift to origin (0,0)
#    x = x - x_shift
#    y = y - y_shift
#    # Convert degrees to radians
#    if units == "DEGREES":
#        angle = math.radians(angle)
#    # Rotation matrix multiplication to get rotated x & y
#    xr = (x * math.cos(angle)) - (y * math.sin(angle)) + x_shift
#    yr = (x * math.sin(angle)) + (y * math.cos(angle)) + y_shift
#    return xr, yr
#
#
#l2_cells = dat[dat=='region_0'].index.values
#adata = spuck[l2_cells]
#x_rot, y_rot = rotate_matrix(adata.obs[['x']].values, adata.obs[['y']].values, angle=10)
#adata.obs['x_rot'] = x_rot
#adata.obs['y_rot'] = y_rot
#
#region_colors = tc.pl.get_default_colors([f'region_{i}' for i in range(len(adata.obs[region_key].value_counts().index))], offset=17)
#tc.pl.scatter(adata, region_key, position_key=('x_rot','y_rot'), colors=region_colors, axes_labels=['X','Y'], noticks=False)
#plt.savefig(f"slideseq_rotated_L2-4.png", bbox_inches='tight')
#
#l4_cells = adata.obs.y_rot[adata.obs.y_rot<6000].index.values
#spuck.obs[region_key][l4_cells] = 'region_6'

l2_cells = dat[dat=='region_0'].index.values
l4_and_below = spuck.obs.y[spuck.obs.y<5950].index.values
l4_cells = [x for x in l2_cells if x in l4_and_below]
spuck.obs[region_key][l4_cells] = 'region_6'

print(spuck.obs[region_key].value_counts())

region_colors = tc.pl.get_default_colors([f'region_{i}' for i in range(len(spuck.obs[region_key].value_counts().index))], offset=17)
tc.pl.scatter(spuck, region_key, position_key=('x','y'), colors=region_colors, axes_labels=['X','Y'], noticks=False)
plt.savefig(f"slideseq_refined_regions.png", bbox_inches='tight')

# rename categories for the region
spuck.obs[region_key] = spuck.obs[region_key].cat.rename_categories({'region_0': 'L2-3', 'region_1': 'L5', 'region_2': 'L6', 'region_3': 'Pia', 'region_4': 'WM', 'region_5': 'L6b', 'region_6': 'L4'})
spuck.obs[region_key] = spuck.obs[region_key].cat.reorder_categories(['Pia','L2-3','L4','L5','L6','L6b','WM'])

print(spuck.obs[region_key].value_counts())

# additional clean up
spuck.obs.loc[(spuck.obs[region_key]=='WM')&(spuck.obs.y>5500),region_key] = 'L4'  # 1 bead
spuck.obs.loc[(spuck.obs[region_key]=='WM')&(spuck.obs.y>4790),region_key] = 'L6b'  # 2 beads
spuck.obs.loc[(spuck.obs[region_key]=='L6')&(spuck.obs.y>5500),region_key] = 'L5'  # 1 bead
spuck.obs.loc[(spuck.obs[region_key]=='L4')&(spuck.obs.y<5500),region_key] = 'L6'  # 1 bead
spuck.obs.loc[(spuck.obs[region_key]=='L5')&(spuck.obs.y<5150)&(700<spuck.obs.x)&(spuck.obs.x<900),region_key] = "L6"  # 1 bead

spuck.obs['final_regions'] = spuck.obs[region_key].copy()
region_key='final_regions'
region_colors = tc.pl.get_default_colors(spuck.obs.final_regions.cat.categories.values, offset=17)
tc.pl.scatter(spuck, region_key, position_key=('x','y'), colors=region_colors, axes_labels=['X','Y'], noticks=False)
plt.savefig(f"slideseq_final_regions.png", bbox_inches='tight')


spuck.write_h5ad("puck_final_regions.h5ad")


# scanpy uses y in puck.obsm['spatial']
genes = ['Satb2','Cux1','Cux2','Rorb','Fezf2','Bcl11b','Tle4','Foxp2','Ctgf','Mbp']
axes = sc.pl.spatial(spuck, color=genes, use_raw=False, spot_size=10, size=3, color_map="magma_r", ncols=4, return_fig=True, wspace=0)
for ax in axes:
    ax.axhline(300, color='lightgray', linestyle='--') # 300 L1
    ax.axhline(850, color='lightgray', linestyle='--') # 550 L2-3
    ax.axhline(1100, color='lightgray', linestyle='--') # 250 L4
    ax.axhline(1700, color='lightgray', linestyle='--') # 600 L5
    ax.axhline(2000, color='lightgray', linestyle='--') # 300 L6
    ax.axhline(2100, color='lightgray', linestyle='--') # 300 L6

plt.savefig("layer_gene_exp.png", dpi=300, bbox_inches='tight')



fig = tc.pl.scatter(spuck, region_key, colors=region_colors, axes_labels=['X','Y'], noticks=False)
for ax in fig.axes:
    ax.axhline(6500, color='lightgray', linestyle='--') # 300 L1
    ax.axhline(5950, color='lightgray', linestyle='--') # 550 L2-3
    ax.axhline(5700, color='lightgray', linestyle='--') # 250 L4
    ax.axhline(5100, color='lightgray', linestyle='--') # 600 L5
    ax.axhline(4800, color='lightgray', linestyle='--') # 300 L6
    ax.axhline(4700, color='lightgray', linestyle='--') # 100 L6b
    ax.axhline(4500, color='lightgray', linestyle='--') # 200 WM


fig.savefig("regions_with_labels.png", bbox_inches='tight')
plt.clf()
plt.close()











#tc.tl.annotation_coordinate(spuck,annotation_key=region_key,result_key='region_dist',max_distance=1000,delta_distance=20,sparse=False);
#tc.pl.annotation_coordinate(spuck,annotation_key=region_key,coordinate_key=('region_dist','region_3'),max_coordinate=1000,delta_coordinate=20, axsize=(3,0.45), colors=region_colors);
md = 2800
tc.tl.annotation_coordinate(spuck,annotation_key=region_key,result_key='region_dist', max_distance=md, delta_distance=20);
cts = ['L2/3 IT CTX-1','L4/5 IT CTX','L5 IT CTX','L6 CT CTX','L6b CTX']
spuck.obsm['anno_neu'] = spuck.obsm[anno_nam].copy()
spuck.obsm['anno_neu'] = spuck.obsm['anno_neu'][cts]
spuck.obsm['anno_neu']= spuck.obsm['anno_neu'].div(spuck.obsm['anno_neu'].sum(axis=1), axis=0)
spuck.obsm['anno_neu']=spuck.obsm['anno_neu'].fillna(0)
tc.pl.annotation_coordinate(spuck,annotation_key='anno_neu',coordinate_key=('region_dist','L2-3'),axsize=(3,0.45), max_coordinate=md, delta_coordinate=20);
plt.savefig(f"slideseq_distances.png", bbox_inches='tight')


spuck.write_h5ad("puck_regions_dist.h5ad")














# not used



# Analyse co-occurrence and neighbourhips
tc.tl.co_occurrence(puck, anno_layer, result_key='type-type',delta_distance=20,max_distance=1000,sparse=False,n_permutation=10)

num_plots = len(puck.obsm[anno_layer].columns)
nc = int(np.ceil(np.sqrt(num_plots)))
nr = int(np.ceil(num_plots / nc))
fig,axs = tc.pl.subplots(n_x=nr, n_y=nc, x_padding=1.5, title="co-occurrence")
tc.pl.co_occurrence(puck, 'type-type', show_only=cts, show_only_center=['Oligo', 'L5 PT CTX', 'L6b CTX'], log_base=2, wspace=0.25) #, ax=axs.reshape(nc*nr,1)[:num_plots,:])
plt.savefig("cell_type_cooccurrence_.png", dpi=300, bbox_inches='tight')

#tc.pl.co_occurrence_matrix(puck, 'type-type',show_only=cts, show_only_center=['Oligo', 'L5 PT CTX', 'L6b CTX'], score_key='z',cmap_vmin_vmax=[-5,5], value_cluster=True, group_cluster=True);
#plt.savefig("cell_type_cooccurrence_matrix.png", dpi=300, bbox_inches='tight')

# Find spatially contiguous regions of comparable expression patterns
pws = [1,1.25,1.5,1.75,2]
ress = [0.75,1,1.25,1.5]
for pw in pws:
    for res in ress:
        region_key = f'regions_pw{pw}_res{res}'
        print(region_key)
        tc.tl.find_regions(puck,key_added=region_key, position_weight=pw, resolution=res);
        puck.obs[region_key] = puck.obs[region_key].map(lambda x: f'region_{x}').astype('category')
        print(puck.obs.regions.value_counts())
        region_colors = tc.pl.get_default_colors([f'region_{i}' for i in range(len(puck.obs[region_key].value_counts().index))], offset=17)
        tc.pl.scatter(puck, region_key, colors=region_colors, axes_labels=['X','Y'], noticks=True)
        plt.savefig(f"slideseq_{region_key}.png", bbox_inches='tight')

# Get regularized distances from these regions
tc.tl.annotation_coordinate(puck,annotation_key='regions',result_key='region_dist',max_distance=500,delta_distance=20,sparse=False);
num_plots = len(spuck.obs.regions.value_counts().index)
nc = int(np.ceil(np.sqrt(num_plots)))
nr = int(np.ceil(num_plots / nc))

fig,axs=tc.pl.subplots(nc,nr,axsize=axsize,x_padding=0.5,y_padding=0.5)
axs=axs.flatten()[:,None]
fig = tc.pl.scatter(puck,'region_dist',cmap='jet', joint=False,axsize=axsize, point_size=1, noticks=True, axes_labels=['X','Y'], ax=axs);
for i in [-4,-2,-1]:
    fig.axes[i].remove()
plt.savefig(f'spatial_distance_from_regions.png', dpi=300, bbox_inches='tight')

fig = tc.pl.annotation_coordinate(puck,annotation_key='labels',coordinate_key=('region_dist','region_2'),colors=labels_colors,max_coordinate=500,delta_coordinate=20, axsize=(3,0.45));
plt.savefig(f'linegraph_distance_from_region2.png', dpi=300, bbox_inches='tight')

fig = tc.pl.compositions(puck, 'labels', 'regions', colors=labels_colors, axsize=(2.4,2.5));
plavefig(f'barchart_cell_type_composition_in_the_regions.png', dpi=300, bbox_inches='tight')

