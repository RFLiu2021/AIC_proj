
import warnings
warnings.filterwarnings("ignore")
import pickle
from SCCAF import *
from NaiveDE import *
import scanpy as sc
import harmonypy as hm
from sklearn.preprocessing import scale
from numpy import unique
import os


import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D 

#os.chdir(r"C:/Users/admin/Documents/Kingston480/MonkeyV1/")
#os.chdir("/media/zhe/Kingston480/MonkeyV1/")
#os.chdir("/media/zhe/A87EE4367EE3FB48/MonkeyV1/")
#os.chdir("/media/zhe/zhe3/MonkeyV1/")

 
def test_code(aaa,ad):
    print('working')
    print(aaa)
    print(ad)
    
##########################
## cluster tools        ##
########################## 
def plot_rank_genes(ad, level_key = 'louvain',plot_group = 0,run_rank = False,  \
                    min_in_group_fraction=0.5, max_out_group_fraction=0.3, min_fold_change=1.5):
# plot all tentative genes to select good marker genes
# input: ad, scanpy dataset; run_rank, True or False, whether rerun select gene function
#input: level_key, the groups to run select genes
#input: plot_group, the cluster that pick genes from

    if run_rank:
        print('redo rank gene')
        sc.tl.rank_genes_groups(ad, groupby=level_key, method='wilcoxon',n_genes = 1000)
        sc.tl.filter_rank_genes_groups(ad,
                               min_in_group_fraction=min_in_group_fraction, #Hodge
                               max_out_group_fraction=max_out_group_fraction,#Hodge
                               min_fold_change=min_fold_change) #Hodge, missing p-value

    all_groups = list(unique(ad.obs[level_key]))
    
    marker_genes_list = ad.uns['rank_genes_groups']['names']
    gene_list = list()
    for gene in marker_genes_list:
        gene0 = gene[plot_group]
        #print (gene)  
        if (str(gene0)=='nan') ==  False:
            gene_list.append(gene0)
            


    plot_gene = gene_list[0:60]
    ad.layers['CPM'] = np.expm1(ad.layers['raw'])
    sc.pl.stacked_violin(ad, groupby=level_key,use_raw = False,\
                     var_names= plot_gene, jitter=False,layer='CPM',
                     swap_axes = True)
    
    return plot_gene

def show_one_cluster(ad,key_level, key):
    #high light selected cluster
    ad.obs['class_show'] = (ad.obs[key_level]== key)
    ad.obs['class_show'] = ad.obs['class_show'].astype(str)
    ad.uns['class_show_colors'] = ['#d3d3d3','#0000FF']
    sc.pl.umap(ad, color = [key_level,'class_show'], title = [key_level,key])

def check_data(ad,plot_type ='umap', plot_gene = ['SNAP25','GAD1','SCL17A7']):
    #show basic information of the dataset ad
    print(ad)
    if plot_type == 'umap':
        sc.pl.umap(ad, color = 'class')
        sc.pl.umap(ad, color = 'subclass')
        sc.pl.umap(ad, color = 'cell_label')
        sc.pl.umap(ad, color = plot_gene)
    
    elif plot_type =='tsne':
        sc.pl.tsne(ad, color = 'class')
        sc.pl.tsne(ad, color = 'subclass')
        sc.pl.tsne(ad, color = 'cell_label')
        sc.pl.tsne(ad, color = plot_gene)
   
###########################
## projection and SCCAF ##
##########################    
def match_gene(adTr, adLn,ref_csv_path = "data/python/human2monkey.csv",
                Ln_col_name = 'Crab-eating macaque gene name',
                Ln_col_ID = 'Crab-eating macaque gene stable ID',
                Tr_species = 'Human',Ln_species = 'Monkey'):
    df = pd.read_csv(ref_csv_path)
    sele = df[Ln_col_name].isna()
    df[Ln_col_name][sele] = df[Ln_col_ID][sele]
    genes = df[['Gene name',Ln_col_name]]
    genes = genes.drop_duplicates('Gene name')
    genes = genes.drop_duplicates(Ln_col_name)
    genes.index = genes[Ln_col_name]
    print(genes.shape)

    genes.loc[genes.index.intersection(adLn.var_names)]['Gene name']

    adLn.var[Tr_species] = genes.loc[genes.index.intersection(adLn.var_names)]['Gene name']

    adLn = adLn[:,adLn.var_names.isin(genes[Ln_col_name])]
    adTr = adTr[:,adTr.var_names.isin(adLn.var[Tr_species])]
    adLn = adLn[:,adLn.var[Tr_species].isin(adTr.var_names)]

    import copy
    df = copy.copy(adLn.var)
    df['MonkeyGene'] = df.index
    df.index = df[Tr_species]

    adLn = adLn[:,df.loc[adTr.var_names]['MonkeyGene'].tolist()]

    
    return [adTr, adLn]


def sele_dataset(adTr, adLn, n_HVG = 2000, Ln_species = 'Monkey', Tr_species = 'Human'):
#generate new dataset with HVG only
#input adTr: trainging set, full dataset; adLn:testing set, full dataset; n_HVG, number of HVG
#output Tr_sub, Ln_sub: splited Training and testing dataset, with combine calculated PCA
    adLn.var[Ln_species] = adLn.var_names
    adLn.var_names = adLn.var[Tr_species]
    
    if scipy.sparse.issparse(adTr.layers['raw']):
        adTr.X = adTr.layers['raw'].todense()
    else: adTr.X = adTr.layers['raw']     
        
    if scipy.sparse.issparse(adLn.layers['raw']):
        adLn.X = adLn.layers['raw'].todense()
    else: adLn.X = adLn.layers['raw']
    

    adTr_selected_gene = adTr.var_names[np.array(np.std(adTr.X, axis=0).argsort())[-n_HVG:][::-1]]
        

    adLn_selected_gene = adLn.var_names[np.array(np.std(adLn.X, axis=0).argsort())[-n_HVG:][::-1]]
    
    select_gene  = unique(list(adTr_selected_gene)+list(adLn_selected_gene))
    
    adTr_sele = adTr[:,select_gene]
    adTr_sele.obs['source']= Tr_species

    adLn_sele = adLn[:,select_gene]
    adLn_sele.obs['source']= Ln_species

    TrLn = adTr_sele.concatenate(adLn_sele)
    

    #sc.pp.scale(TrLn)
    #sc.tl.pca(TrLn)
    #sc.pl.pca(TrLn, color = ['marker1','source'])
    
    #Tr_sub = TrLn[TrLn.obs['source']==Tr_source,:]
    #Ln_sub = TrLn[TrLn.obs['source']==Ln_source,:]
     
    return TrLn

def run_PCA_Harmony(TrLn, run_Harmony = True, batch_key = 'source',theta = 5,rep = 3):
#run PCA and Harmony (optional) 
#input TrLn, combined dataset from Training set and testing set, generated by function sele_dataset
#output Tr_sub, Ln_sub: splited Training and testing dataset, with combine calculated PCA and Harmony, in layer obsm
    
    TrLn.X = TrLn.layers['raw']
    if scipy.sparse.issparse(TrLn.X): 
        X = TrLn.X.todense().T
    else:
        X = TrLn.T
    #TrLn.obsm['X_pca'] = sc.tl.pca(np.array(scale(X.T)), n_comps=100)
    sc.pp.scale(TrLn)
    sc.tl.pca(TrLn)
    #sc.pl.pca(TrLn, color = ['marker1','source'])
    
    #Tr_sub = TrLn[TrLn.obs['source']==Tr_source,:]
    #Ln_sub = TrLn[TrLn.obs['source']==Ln_source,:]
    
    if run_Harmony:
        
        data_mat = TrLn.obsm['X_pca']
        meta_data = TrLn.obs

        TrLn.obsm['X_pca'].shape
        #ho = hm.run_harmony(data_mat, meta_data, 'species')
        ho = hm.run_harmony(data_mat, meta_data, batch_key,theta = theta)
        res = ho.Z_corr
        TrLn.obsm['X_harmony'] = res.T
        sc.pp.neighbors(TrLn)
        sc.tl.umap(TrLn)
        #sc.tl.tsne(humk, n_jobs = 20)
        sc.pl.umap(TrLn, color = [batch_key,'subclass'])
        #sc.pl.tsne(humk,color = ['species','marker1','marker2'])

        for i in range(1,rep):
        
            print("this is the " + str(i) + " time of harmony" )
    
            data_mat = TrLn.obsm['X_harmony']
            meta_data = TrLn.obs

            #ho = hm.run_harmony(data_mat, meta_data, 'species')
            ho = hm.run_harmony(data_mat, meta_data, batch_key, theta = 5)
            res = ho.Z_corr
            TrLn.obsm['X_harmony'] = res.T
            sc.pp.neighbors(TrLn,use_rep = 'X_harmony')
            sc.tl.umap(TrLn)
            #sc.tl.tsne(humk, n_jobs = 20)
            sc.pl.umap(TrLn, color = ['source','subclass'])
            #sc.pl.tsne(humk,color = ['species','marker1','marker2'])
 
    return TrLn

def data_subset(adTr,key = 'marker1', n=500, frac=1):
    X_train, X_test, y_train, y_test = train_test_split_per_type(adTr.X, adTr.obs[key], n=n, frac=frac)
    adTr_sub = adTr[y_train.index,:]
    return adTr_sub

def label_projection (adTr, adLn, key = 'marker1', mode = 'X', n=1000, frac=0.8, plot_tsne = False):
# SCCAF projectin
#input: adTr: training set; adLn, testing set;key, label for learning; mode, X or PCA or harmony PCA
#output: dataframe of probablity of each cell, each label
    
    if mode == 'X':
        Tr_matrix = adTr.X
        Ln_matrix = adLn.X
        
    elif mode == 'X_pca':
            
        if scipy.sparse.issparse(adTr.obsm['X_pca']): 
            Tr_matrix  = adTr.obsm['X_pca'].todense()
        else:
            Tr_matrix  = adTr.obsm['X_pca']
            
        if scipy.sparse.issparse(adLn.obsm['X_pca']): 
            Ln_matrix  = adLn.obsm['X_pca'].todense()
        else:
            Ln_matrix  = adLn.obsm['X_pca']
            
        
    elif mode == 'X_harmony':    
        if scipy.sparse.issparse(adTr.obsm['X_harmony']): 
            Tr_matrix  = adTr.obsm['X_harmony'].todense()
        else:
            Tr_matrix  = adTr.obsm['X_harmony']
            
        if scipy.sparse.issparse(adLn.obsm['X_harmony']): 
            Ln_matrix  = adLn.obsm['X_harmony'].todense()
        else:
            Ln_matrix  = adLn.obsm['X_harmony']
        
        
    X_train, X_test, y_train, y_test = train_test_split_per_type(Tr_matrix, adTr.obs[key], n=n, frac=frac)
    clf = LogisticRegression(random_state=1, penalty='l1', C=0.5, solver = 'saga')
    clf.fit(X_train, y_train)

    if plot_tsne:
        sc.pl.umap(adTr, color = ['marker1','marker2'])

    adTr.obs[key+'_predict'] = clf.predict(Tr_matrix)
    if plot_tsne:
        sc.pl.tsne(adTr, color = [key,key+'_predict'])

    adLn.obs[key + '_predict'] = clf.predict(Ln_matrix)
    #sc.pl.tsne(adLn, color = [key,key+'_predict'])
    df_prob = pd.DataFrame(clf.predict_proba(Ln_matrix),index = adLn.obs_names, columns = clf.classes_)
    
    return df_prob, adLn


def get_prob_median(df_prob_rep):
#caluculating probability median from the big repeat df
#input: pd.Series each obj is a prob dataframe of each repeat
#output: df with median of each cell, each prob 
    df_prob_median = pd.DataFrame()
    cell_type_list = list(df_prob_rep[0].columns)

    for cell_type in cell_type_list:
        print('calculating ' + cell_type + '...')
        for rep_ind in list(df_prob_rep.index):
            #print(rep_ind)
            df_prob_type = pd.DataFrame()
            df_prob_type_median = pd.DataFrame()
            try: df_prob_type[rep_ind]=df_prob_rep[rep_ind][cell_type]
            except: 
                rep_ind_sub = list(df_prob_rep.index)[0]
                df_prob_type[rep_ind]=df_prob_rep[rep_ind_sub][cell_type]
            df_prob_type_median = np.median(df_prob_type, axis = 1)
        df_prob_median[cell_type]=df_prob_type_median
    return df_prob_median

########################
## plot hodge heatmap##
#######################
def two_species_heatmap(ad, species_1 = 'FM', species_2 = 'HS',species_1_key = 'subclass', species_2_key = 'marker2',\
                        louvain = 0,figure_path = 'test_heatmap.png'):
#generate hodge fig 5d heatmap
#input: ad: pre_merged, harmony corrected two species data, with .obs['species'] mark the species
#input : species_key, which .observe to generate the heatmap
#input: louvain, the louvain resolution for cluster the merged data
#input: figure_path, the path to save the heatplot, if figure_path = 0, do not save figure

    #prepare data
    if not louvain ==0:
        sc.pp.neighbors(ad,  metric='euclidean',use_rep = 'X_harmony' )
        sc.tl.louvain(ad, resolution = louvain,key_added = 'louvain')

    sc.pl.umap(ad, color = ['species','louvain'])
    
    ad_1 = ad[ad.obs['species'].isin([species_1]),:] 
    sc.pl.umap(ad_1, color = [species_1_key,'louvain'], legend_loc = 'on data')
    
    ad_2 = ad[ad.obs['species'].isin([species_2]),:]
    sc.pl.umap(ad_2, color = [species_2_key,'louvain'], legend_loc = 'on data')
    
    df_1 = pd.crosstab(ad_1.obs[species_1_key], ad_1.obs['louvain'], normalize ='index')
    df_2 = pd.crosstab(ad_2.obs[species_2_key], ad_2.obs['louvain'], normalize ='index')


    #generate heatmap matrix
    mk_cluster_all = df_1.index
    hu_cluster_all = df_2.index 
    low_sum_matrix = pd.DataFrame(index=mk_cluster_all, columns = hu_cluster_all)

    low_sum_matrix.index.name = species_1 +'_cluster'
    low_sum_matrix.columns.name = species_2 +'_cluster'
    
    
    low_sum_matrix = pd.DataFrame(index=mk_cluster_all, columns = hu_cluster_all)
    for mk_cluster in mk_cluster_all:
        for hu_cluster in hu_cluster_all:
        
            two_row = np.column_stack([df_1.loc[mk_cluster],df_2.loc[hu_cluster]])
            low_sum = sum(two_row.min(axis=1))
            low_sum_matrix.loc[mk_cluster, hu_cluster] = low_sum
            
    low_sum_matrix.to_csv('data/python/Homology_map/temp.csv')
    del(low_sum_matrix)
    low_sum_matrix =pd.read_csv('data/python/Homology_map/temp.csv')
    low_sum_matrix.set_index((species_1+'_cluster'), inplace = True)
    # for some reason, direct use matrix does not work.... save and read did the trick...
    
    low_sum_matrix_sort = low_sum_matrix.copy()

    ss = sns.clustermap(low_sum_matrix_sort, cmap = 'Greys',  cbar = True, col_cluster = False,  xticklabels=1, yticklabels=1,  method =  'single')
    
    hu_cluster_all = low_sum_matrix_sort.columns
    low_sum_matrix_sort1 = low_sum_matrix_sort
    low_sum_matrix_sort2 = low_sum_matrix_sort1.iloc[ss.dendrogram_row.reordered_ind]  
    
    low_sum_matrix1 =  low_sum_matrix_sort2
    for hu_cluster in hu_cluster_all:
        line = list(low_sum_matrix1.loc[:,hu_cluster])
        ss = np.array(line)
        tmp = []
        for kk in  range(len(ss)):
            if kk < len(ss) - 5:
                tmp.append( ss[kk]+ss[kk+1]/2+ss[kk+2]/3+ss[kk+3]/4+ss[kk+4]/5)
            elif kk < len(ss) - 4: 
                tmp.append( ss[kk]+ss[kk+1]/2+ss[kk+2]/3+ss[kk+3]/4)
            elif kk < len(ss) - 3: 
                tmp.append( ss[kk]+ss[kk+1]/2+ss[kk+2]/3 )
            elif kk < len(ss) - 2: 
                tmp.append( ss[kk]+ss[kk+1]/2 )
            elif kk < len(ss) - 1: 
                tmp.append( ss[kk] )
            elif kk < len(ss): 
                tmp.append(ss[kk]+ss[kk]/2+ss[kk]/3+ss[kk]/4+ss[kk]/5+ss[kk]/6)
        low_sum_matrix_sort2.loc['max_order',hu_cluster]= float(tmp.index(max(tmp)))
        del tmp
    low_sum_matrix_sort2 = low_sum_matrix_sort2.sort_values(by = ['max_order'],axis=1)
    plt.figure(figsize = (10,10))
    graph = sns.heatmap(low_sum_matrix_sort2[0:-1], cmap="Greys", cbar=True, xticklabels=1,yticklabels=1)
    
    
    if not figure_path == 0:
        fig = graph.get_figure()
        fig.savefig(figure_path)
  
      
########################
###  prepare figure ####
#######################

def plot_color(color_list):
    """
    Helper function to plot data with associated colormap.
    """
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    colormaps = ListedColormap(color_list)
    np.random.seed(19680801)
    data = np.random.randn(30, 30)
    n = len(colormaps)
    fig, axs = plt.subplots(1, n, figsize=(n * 2 + 2, 3),
                            constrained_layout=True, squeeze=False)
    for [ax, cmap] in zip(axs.flat, colormaps):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=-4, vmax=4)
        fig.colorbar(psm, ax=ax)
    plt.show()  
    


def ad_3D_plot (ad, basis = 'diffmap',groupby = 'subclass',\
                components = '1,2,3', dotsize = 10,title = 'title'):
# scanpy 3D scatter seems broken, write one for my own
# using matpltlib scatter3D
# input: ad: annadata, basis: type of plot, can be umap or tsne
# input: dotsize: dotsize for the scatter
# return plot ax
    plot_df = pd.DataFrame(ad.obsm[('X_'+basis)], index = ad.obs_names)
    plot_groups = unique(ad.obs[groupby])
    
    comp_list = components.split(',')
    comp_list = [int(i) for i in comp_list]
    color_list = ad.uns[(groupby+'_colors')]
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    for i in range(0,len(plot_groups)):
        plot_group = plot_groups[i]
        cell_list = list(ad[ad.obs[groupby]==plot_group].obs_names)
        zdata = plot_df.loc[cell_list,comp_list[0]]
        ydata = plot_df.loc[cell_list,comp_list[1]]
        xdata = plot_df.loc[cell_list,comp_list[2]]
        ax.scatter3D(xdata,ydata,zdata,c = color_list[i], s = dotsize,label = plot_group, alpha = 0.5)
    plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8, bbox_to_anchor=(0, 0))
    plt.title = title
    ax.set_xlabel((basis +' ' + str(comp_list[2])))
    ax.set_ylabel((basis +' ' + str(comp_list[1])))
    ax.set_zlabel((basis +' ' + str(comp_list[0])))

    ax.set_xlim(min(plot_df.loc[:,comp_list[2]]), max(plot_df.loc[:,comp_list[2]]))
    ax.set_ylim(min(plot_df.loc[:,comp_list[1]]), max(plot_df.loc[:,comp_list[1]]))
    ax.set_zlim(min(plot_df.loc[:,comp_list[0]]), max(plot_df.loc[:,comp_list[0]]))

    ax.set_title(title)
    return ax






def ad_3D_plot_intensity (ad, basis = 'diffmap',color = 'latent_time',\
                components = '1,2,3', dotsize = 1,title = 'title'):
# scanpy 3D scatter seems broken, write one for my own
# using matpltlib scatter3D
# input: ad: annadata, basis: type of plot, can be umap or tsne
# input: dotsize: dotsize for the scatter
# return plot ax
    plot_df = pd.DataFrame(ad.obsm[('X_'+basis)], index = ad.obs_names)
    
    comp_list = components.split(',')
    comp_list = [int(i) for i in comp_list]
    
    fig_brightness = ad.obs[color]
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    zdata = plot_df.loc[:,comp_list[2]]
    ydata = plot_df.loc[:,comp_list[1]]
    xdata = plot_df.loc[:,comp_list[0]]
    
    p = ax.scatter3D(xdata,ydata,zdata,s = 1,c = fig_brightness,cmap='gnuplot')

    plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8, bbox_to_anchor=(0, 0))
    plt.title = title
    ax.set_xlabel((basis +' ' + str(comp_list[0])))
    ax.set_ylabel((basis +' ' + str(comp_list[1])))
    ax.set_zlabel((basis +' ' + str(comp_list[2])))

    ax.set_xlim(min(plot_df.loc[:,comp_list[0]]), max(plot_df.loc[:,comp_list[0]]))
    ax.set_ylim(min(plot_df.loc[:,comp_list[1]]), max(plot_df.loc[:,comp_list[1]]))
    ax.set_zlim(min(plot_df.loc[:,comp_list[2]]), max(plot_df.loc[:,comp_list[2]]))

    ax.set_title(title)
    
    #plt.colorbar(p)
    return ax


def plot_one_GO_Corr (res_mk, res_hu, res_mm, key1, key2, gene_list, saveFig = False):
    #key1 = 'PVALB'
    #key2 = 'PVALB'

    gene_list_overlap = (list(set(gene_list) & set(adhu.var_names)))
    #len(gene_list_overlap) 

    res_mk_sele = res_mk.loc[key1,gene_list_overlap]
    res_hu_sele = res_hu.loc[key2,gene_list_overlap]
    res_mm_sele = res_mm.loc[key2,gene_list_overlap]

    rho_humk, pval_humk = stats.spearmanr(res_mk_sele, res_hu_sele)
    rho_mmmk, pval_mmmk = stats.spearmanr(res_mk_sele, res_mm_sele)
    rho_mmhu, pval_mmhu = stats.spearmanr(res_hu_sele, res_mm_sele)

    figsize(9,3)
    matplotlib.rcParams.update({'font.size': 7})
    fig = plt.figure(constrained_layout=True)
    gs = GridSpec(1, 3, figure=fig)
    #sns.despine()

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    ax1.scatter(res_mk_sele, res_hu_sele, color = 'blue', s = 5, alpha = 0.5)
    ax2.scatter(res_mk_sele, res_mm_sele, color = 'black', s = 5, alpha = 0.5)
    ax3.scatter(res_hu_sele, res_mm_sele, color = 'green', s = 5, alpha = 0.5)

    
    ax1.text (0.5*max(res_mk_sele),0.8*max(res_hu_sele),'rho = '+str(format(rho_humk,'.3f')), color = 'black')
    ax2.text (0.5*max(res_mk_sele),0.8*max(res_hu_sele),'rho = '+str(format(rho_mmmk,'.3f')), color = 'black')
    ax3.text (0.5*max(res_mk_sele),0.8*max(res_hu_sele),'rho = '+str(format(rho_mmhu,'.3f')), color = 'black')

    ax1.title.set_text(f'axon geneset of {key1} subclass')
    ax2.title.set_text(f'axon geneset of {key1} subclass')
    ax3.title.set_text(f'axon geneset of {key1} subclass')

    x1 = np.array(res_mk_sele,dtype=float)
    y1 = np.array(res_hu_sele,dtype=float)
    z1 = numpy.polyfit(x1, y1, 1,)
    p1 = numpy.poly1d(z1)
    ax1.plot(x1,p1(x1),"b----")

    x2 = np.array(res_mk_sele,dtype=float)
    y2 = np.array(res_mm_sele,dtype=float)
    z2 = numpy.polyfit(x2, y2, 1,)
    p2 = numpy.poly1d(z2)
    ax2.plot(x2,p2(x2),"k----")

    x3 = np.array(res_hu_sele,dtype=float)
    y3 = np.array(res_mm_sele,dtype=float)
    z3 = numpy.polyfit(x3, y3, 1,)
    p3 = numpy.poly1d(z3)
    ax3.plot(x3,p3(x3),"g----")

    ax1.spines['bottom'].set_color('#000000')
    ax1.spines['left'].set_color('#000000') 
    ax1.spines['top'].set_color('none')
    ax1.spines['right'].set_color('none')

    ax2.spines['bottom'].set_color('#000000')
    ax2.spines['left'].set_color('#000000') 
    ax2.spines['top'].set_color('none')
    ax2.spines['right'].set_color('none')

    ax3.spines['bottom'].set_color('#000000')
    ax3.spines['left'].set_color('#000000') 
    ax3.spines['top'].set_color('none')
    ax3.spines['right'].set_color('none')

    ax1.set_xlabel('Monkey PVALB')
    ax1.set_ylabel('Human PVALB')
    ax2.set_xlabel('Monkey PVALB')
    ax2.set_ylabel('Mouse PVALB')
    ax3.set_xlabel('Houman PVALB')
    ax3.set_ylabel('Mouse PVALB')

    #ax1.set_xlim([-1,1])
    #ax1.set_ylim([-1,1])

    #ax2.set_xlim([-1,1])
    #ax2.set_ylim([-1,1])

    #ax3.set_xlim([-1,1])
    #ax3.set_ylim([-1,1])

    if saveFig:
        plt.savefig(f'Figure/Fig6/glut_{key1}_{key2}_axon_scatter.png', dpi = 600)
    return fig


######################
## plot stacked bar###
######################

def sc_plot_stacke_bar(ad, \
                       y_key="predict_plot_direct",\
                       x_key="cell_type",\
                       cluster_palette=None,xlabel_rotation=0,\
                       savePath = 'Figure/Fig_NG_SCCAF/stackedbar_',
                       rotation = 0, fig_size = (8,5)
                      ):
    
    #plot stacked bar plot for SCCAF projection results
    #implement from  gist.github.com/wflynny 
    #adata : AnnData object
    #cluster_key : original cell type key
    #predict_key : SCCAF projected key
    # cluster_palette: list of color, same sequence as category
    
    sizes = ad.obs.groupby([y_key, x_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index() 
    props = props.pivot(columns=x_key, index=y_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)
    
    #plot figure
    
    fig, ax = plt.subplots(dpi=600, figsize = fig_size)
    fig.patch.set_facecolor("white")
    
    cmap = None
    
    if cluster_palette is not None:
        cmap = sns.palettes.blend_palette(
            cluster_palette, 
            n_colors=len(cluster_palette), 
            as_cmap=True)
   
    cluster_props = props.copy()
    cluster_props.plot(
        kind="bar", 
        stacked=True, 
        ax=ax, 
        legend=None, 
        colormap=cmap
    )
    
    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False, title=" ")
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=rotation)
    ax.set_xlabel(cluster_props.index.name.capitalize())
    ax.set_ylabel("Proportion")
    fig.tight_layout() 
    
    plt.savefig(savePath+y_key+'.pdf', dpi = 600)
    plt.savefig(savePath+y_key+'.png', dpi = 600)
   
    
    return fig
