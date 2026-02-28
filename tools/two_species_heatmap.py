
########################
## plot hodge heatmap##
#######################

import numpy as np
import pandas as pd
import scanpy as sc
#import scvi
import seaborn as sns
import matplotlib.pyplot as plt



def two_species_heatmap(ad, species_1 = 'FM', species_2 = 'HS',species_1_key = 'subclass', species_2_key = 'marker2',\
                        louvain = 0,figure_path = 'test_heatmap.png'):
#generate hodge fig 5d heatmap
#input: ad: pre_merged, harmony corrected two species data, with .obs['species'] mark the species
#input : species_key, which .observe to generate the heatmap
#input: louvain, the louvain resolution for cluster the merged data
#input: figure_path, the path to save the heatplot, if figure_path = 0, do not save figure

    #prepare data
    if not louvain ==0:
        sc.pp.neighbors(ad,  metric='euclidean',use_rep = 'X_harmonypca' )
        sc.tl.louvain(ad, resolution = louvain,key_added = 'louvain')

    sc.pl.umap(ad, color = ['species','louvain'])
    
    ad_1 = ad[ad.obs['species'].isin([species_1]),:] 
    sc.pl.umap(ad_1, color = [species_1_key,'louvain'], legend_loc = 'on data')
    
    ad_2 = ad[ad.obs['species'].isin([species_2]),:]
    sc.pl.umap(ad_2, color = [species_2_key,'louvain'], legend_loc = 'on data')
    
    df_1 = pd.crosstab(ad_1.obs[species_1_key], ad_1.obs['louvain'], normalize ='index')
    df_2 = pd.crosstab(ad_2.obs[species_2_key], ad_2.obs['louvain'], normalize ='index')

    df_1.columns = df_1.columns.tolist()
    df_2.columns = df_2.columns.tolist()

    #add missing louvain
    dif_list = list(set(list(df_2.columns)) ^set(list(df_1.columns)))
    for dif in dif_list:
        if dif not in list(df_1.columns):
            df_1[dif] = 0
        if dif not in list(df_2.columns):
            df_2[dif] = 0  
    df_2 = df_2[df_1.columns] #reorder column sequence, important!

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
            
    low_sum_matrix.to_csv('data_tem/temp.csv')
    del(low_sum_matrix)
    low_sum_matrix =pd.read_csv('data_tem/temp.csv')
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
    plt.figure(figsize = (15,15))
    graph = sns.heatmap(low_sum_matrix_sort2[0:-1], cmap="Greys", cbar=True, xticklabels=1,yticklabels=1, linewidth = 0.01, linecolor = 'gray')
    
    plot_matrix = low_sum_matrix_sort2[0:-1]
    
    if not figure_path == 0:
        fig = graph.get_figure()
        fig.savefig(figure_path)
        
    return(plot_matrix)    
      
