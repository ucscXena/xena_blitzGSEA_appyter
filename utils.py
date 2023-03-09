from rpy2 import robjects
from rpy2.robjects import r, pandas2ri

# Basic libraries
import pandas as pd
import math
import os
import urllib3, certifi
import urllib.parse
import requests, json
import sys
import geode
import random
from time import sleep
import time
import numpy as np
import warnings
import base64  

# Visualization
from dash import Dash, dash_table, dcc, html
from dash.dependencies import Input, Output
import plotly
from plotly import tools
import plotly.express as px
import plotly.graph_objs as go
from jupyter_dash import JupyterDash
plotly.offline.init_notebook_mode() # To embed plots in the output cell of the notebook

import IPython
from IPython.display import HTML, display, Markdown, IFrame

import chart_studio
import chart_studio.plotly as py

# Data analysis
from itertools import combinations
import scipy.spatial.distance as dist
import scipy.stats as ss
from sklearn.decomposition import PCA
from maayanlab_bioinformatics.normalization.quantile import quantile_normalize
from maayanlab_bioinformatics.dge.characteristic_direction import characteristic_direction

# Import umap
from sklearn.manifold import TSNE

# blitzGSEA
import blitzgsea as blitz

http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())

def check_files(fname):
    if fname == "":
        raise IOError
    if fname.endswith(".txt") == False and fname.endswith(".csv") ==False and fname.endswith(".tsv")==False:
        raise IOError

def find_remote_file(fname):
    r = http.request('HEAD', fname)
    http.clear()
    if r.status != 200:
        fname = fname +".gz"
        r = http.request('HEAD', fname)
        http.clear()
        if r.status != 200:
            raise IOError
        else:
            return fname
    else:
        return fname

def to_float(x):
    try:
        return float(x)
    except:
        return float('nan')

def parse_xena_expr(fin, samples):
    columns = fin.readline()[:-1].split('\t')
    # find the columns to keep
    keep_pos =[]
    keep_samples =[]
    ids =[]
    for i in range(1, len(columns)):
        if columns[i] in samples:
            keep_pos.append(i)
            keep_samples.append(columns[i])
    big_list = []
    while 1:
        line = fin.readline()
        if line == '':
            break
        data = line[:-1].split('\t')
        id = data[0]
        ids.append(id)
        big_list.append(list(map(to_float, map(data.__getitem__, keep_pos))))
    fin.close()
    df = pd.DataFrame(big_list, index = ids, columns=keep_samples)
    return df

def check_subgroups (control_group, case_group):
    if len(case_group) == 0:
        return 'Subgroup 1 has no samples.'
    if len(control_group) == 0:
        return 'Subgroup 2 has no samples.'
    if control_group == case_group:
        return 'You selected the same categories \"' + '_'.join(control_group)+ '\" for both subgroups.'
    if len(list(set(control_group) & set(case_group))) != 0:
        return 'There can not be shared samples between the two subgroups.'
    return 0

def check_probe_level(probemap_df):
    counter = 0
    counter_gene_level = 0
    gene_pos = 1
    for probe in probemap_df.index:
        if type(probemap_df.loc[probe][gene_pos]) == str and probe == probemap_df.loc[probe][gene_pos]:
            counter_gene_level += 1
        counter += 1
        if counter == 1000:
            break
    if counter_gene_level/counter > 0.5:
        probe_level = "HUGO"
    else:
        probe_level = "probe"
    return probe_level

def convert_to_hugo(expr_df, probemap_df):
    '''
    convert expr_df matrix to HUGO gene_level expr_df using probemap_df
    '''
    big_list= []
    gene_list = []
    gene_pos = 1
    for probe in expr_df.index:
        try:
            if type(probemap_df.loc[probe][gene_pos]) == str:
                genes = probemap_df.loc[probe][gene_pos].split(",")
                gene_list += genes
                for gene in genes:
                    big_list.append(list(expr_df.loc[probe]))
        except:
            pass
    df_copy = pd.DataFrame(big_list, index = gene_list, columns=expr_df.columns)
    df_copy = df_copy.groupby(df_copy.index).mean()
    return df_copy

def build_meta_df(samples, values, codes):
    big_list= []
    for i in range (0, len(samples)):
        sample = samples[i]
        data = values[i]
        if data is None:
            value = ''
        else:
            value = codes[values[i]]
        big_list.append([sample, value])
    return pd.DataFrame(big_list, columns=['sample', 'category'])

def check_df(df, col):
    if col not in df.columns:
        raise IOError
        
def CPM(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = (data/data.sum())*10**6
        data = data.fillna(0)
        
    return data
def logCPM(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = (data/data.sum())*10**6
        data = data.fillna(0)
        data = np.log10(data+1)

    # Return
    return data
def log(data):

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data = data.fillna(0)
        data = np.log10(data+1)

    return data
def qnormalization(data):
  
    X_quantile_norm = quantile_normalize(data)
    return X_quantile_norm  

def normalize(dataset, current_dataset, logCPM_normalization, log_normalization, z_normalization, q_normalization):
    normalization = current_dataset
    if logCPM_normalization == True:  
        data = dataset[normalization]
        normalization += '+logCPM'
        dataset[normalization] = logCPM(data)
        
    if log_normalization == True:    
        data = dataset[normalization]
        normalization += '+log'
        dataset[normalization] = log(data)
        
    if z_normalization == True:
        data = dataset[normalization]
        normalization += '+z_norm'    
        dataset[normalization] = data.T.apply(ss.zscore, axis=0).T.dropna()

    if q_normalization == True:
        data = dataset[normalization]
        normalization += '+q_norm'
        dataset[normalization] = qnormalization(data)
    return dataset, normalization

def create_download_link(df, title = "Download CSV file: {}", filename = "data.csv"):  
    df.to_csv(filename)
    html = "<a href=\"./{}\" target='_blank'>{}</a>".format(filename, title.format(filename))
    
    return HTML(html)

def display_link(url):
    raw_html = '<a href="%s" target="_blank">%s</a>' % (url, url)
    return display(HTML(raw_html))

def display_object(counter, caption, df=None, istable=True):
    if df is not None:
        display(df)
    if istable == True:
        display(Markdown("*Table {}. {}*".format(counter, caption)))
    else:
        display(Markdown("*Figure {}. {}*".format(counter, caption)))
    counter += 1
    return counter

def run_dimension_reduction(dataset, meta_id_column_name, method='PCA', normalization='logCPM', nr_genes=2500, filter_samples=True, plot_type='interactive'):
    # Get data
    before_norm = normalization.replace("+z_norm", "").replace("+q_norm", "")
    top_genes = dataset[before_norm].var(axis=1).sort_values(ascending=False)
    
    expression_dataframe = dataset[normalization]
    
    # Filter columns
    if filter_samples and dataset.get('signature_metadata'):
        selected_samples = [sample for samples in list(dataset['signature_metadata'].values())[0].values() for sample in samples]
        expression_dataframe = expression_dataframe[selected_samples]

    # Filter rows
    expression_dataframe = expression_dataframe.loc[top_genes.index[:nr_genes]]
    result=None
    axis = None
    if method == 'PCA':
        # Run PCA
        pca=PCA(n_components=3)
        result = pca.fit(expression_dataframe)

        # Get Variance
        axis = ['PC'+str((i+1))+'('+str(round(e*100, 1))+'% var. explained)' for i, e in enumerate(pca.explained_variance_ratio_)]
    #elif method == "UMAP":
    #    u = umap.UMAP(n_components=3)
    #    result = u.fit_transform(expression_dataframe.T)
    #    axis = ['UMAP 1', 'UMAP 2', 'UMAP 3']
    elif method == "t-SNE":
        u = TSNE(n_components=3)
        result = u.fit_transform(expression_dataframe.T)
        axis = ['t-SNE 1', 't-SNE 2', 't-SNE 3']

    # Add signature groups
    if dataset.get('signature_metadata'):
        A_label, B_label = list(dataset.get('signature_metadata').keys())[0].split(' vs ')
        col = []
        group_dict = list(dataset.get('signature_metadata').values())[0]
        for gsm in dataset['sample_metadata'].index:
            if gsm in group_dict['A']:
                col.append(A_label)
            elif gsm in group_dict['B']:
                col.append(B_label)
            else:
                col.append('Other')
        dataset['dataset_metadata']['Sample Group'] = col       

    # Return
    results = {'method': method, 'result': result, 'axis': axis, 
                   'dataset_metadata': dataset['dataset_metadata'][dataset['dataset_metadata'][meta_id_column_name] == expression_dataframe.columns], 
                   'nr_genes': nr_genes, 
                   'normalization': normalization, 'signature_metadata': dataset.get('signature_metadata'), 
                   'plot_type': plot_type}
    return results

def plot_samples(pca_results, meta_id_column_name, meta_class_column_name, counter, plot_type='interactive',):
    pca_transformed = pca_results['result']
    axis = pca_results['axis']
    meta_df = pca_results['dataset_metadata']
    
    if pca_results["method"] == "PCA":
        meta_df['x'] = pca_transformed.components_[0]
        meta_df['y'] = pca_transformed.components_[1]
        meta_df['z'] = pca_transformed.components_[2]
    else:
        meta_df['x'] = [x[0] for x in pca_transformed]
        meta_df['y'] = [x[1] for x in pca_transformed]
        meta_df['z'] = [x[2] for x in pca_transformed]
        

    caption = '3D {} plot for samples using {} genes having largest variance. \
    The figure displays an interactive, three-dimensional scatter plot of the data. \
    Each point represents an gene expression sample. \
    Samples with similar gene expression profiles are closer in the three-dimensional space. \
    If provided, sample groups are indicated using different colors, allowing for easier interpretation of the results.'.format(pca_results["method"], pca_results['nr_genes'])
    display(IPython.core.display.HTML('''
            <script src="/static/components/requirejs/require.js"></script>
            <script>
              requirejs.config({
                paths: {
                  base: '/static/base',
                  plotly: 'https://cdn.plot.ly/plotly-latest.min.js?noext',

                },
              });
            </script>
            '''))

    classes = meta_df[meta_class_column_name].unique().tolist()
    SYMBOLS = ['circle', 'square']
    
    if len(classes) > 10:
        def r(): return random.randint(0, 255)
        COLORS = ['#%02X%02X%02X' % (r(), r(), r())
                          for i in range(len(classes))]
    else:
        COLORS = [
            '#1f77b4',
            '#ff7f0e',
            '#2ca02c',
            '#d62728',
            '#9467bd',
            '#8c564b',
            '#e377c2',
            '#7f7f7f',
            '#bcbd22',
            '#17becf',
            ]

    data = [] # To collect all Scatter3d instances
    for (cls), meta_df_sub in meta_df.groupby([meta_class_column_name]):
        # Iteratate through samples grouped by class
        display_name = '%s' % (cls)
        # Initiate a Scatter3d instance for each group of samples specifying their coordinates
        # and displaying attributes including color, shape, size and etc.
        trace = go.Scatter3d(
            x = meta_df_sub['x'],
            y = meta_df_sub['y'],
            z = meta_df_sub['z'],

            text=meta_df_sub[meta_id_column_name],
            mode='markers',
            marker=dict(
                size=10,
                color=COLORS[classes.index(cls)], # Color by infection status
                opacity=.8,
            ),
#             name=meta_df_sub[meta_id_column_name]
            name=display_name,
        )

        data.append(trace)

    # Configs for layout and axes
    
    layout=dict(height=1000, width=1000, 
                title='3D {} plot for samples'.format(pca_results["method"]),
                scene=dict(
                    xaxis=dict(title=axis[0]),
                    yaxis=dict(title=axis[1]),
                    zaxis=dict(title=axis[2])
                    )
    )
    

    fig=dict(data=data, layout=layout)
    if plot_type == "interactive":
        plotly.offline.iplot(fig)
    else:
        py.image.ishow(fig)
    display(Markdown("*Figure {}. {}*".format(counter, caption)))
    counter += 1
    return counter
    
#############################################
########## 2. Plot
#############################################

def plot_2D_scatter(x, y, text='', title='', xlab='', ylab='', hoverinfo='text', color='black', colorscale='Blues', size=8, showscale=False, symmetric_x=False, symmetric_y=False, pad=0.5, hline=False, vline=False, return_trace=False, labels=False, plot_type='interactive', de_type='ma'):
    range_x = [-max(abs(x))-pad, max(abs(x))+pad]if symmetric_x else None
    range_y = [-max(abs(y))-pad, max(abs(y))+pad]if symmetric_y else None
    trace = go.Scattergl(x=x, y=y, mode='markers', text=text, hoverinfo=hoverinfo, marker={'color': color, 'colorscale': colorscale, 'showscale': showscale, 'size': size})
    if return_trace:
        return trace
    else:
        if de_type == 'ma':
            annotations = [
                {'x': 1, 'y': 0.1, 'text':'<span style="color: blue; font-size: 10pt; font-weight: 600;">Down-regulated in '+labels[0]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'right', 'yanchor': 'top'},
                {'x': 1, 'y': 0.9, 'text':'<span style="color: red; font-size: 10pt; font-weight: 600;">Up-regulated in '+labels[0]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'right', 'yanchor': 'bottom'}
            ] if labels else []
        elif de_type == 'volcano':
            annotations = [
                {'x': 0.25, 'y': 1.07, 'text':'<span style="color: blue; font-size: 10pt; font-weight: 600;">Down-regulated in '+labels[0]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'},
                {'x': 0.75, 'y': 1.07, 'text':'<span style="color: red; font-size: 10pt; font-weight: 600;">Up-regulated in '+labels[0]+'</span>', 'showarrow': False, 'xref': 'paper', 'yref': 'paper', 'xanchor': 'center'}
            ] if labels else []
        layout = go.Layout(title=title, xaxis={'title': xlab, 'range': range_x}, yaxis={'title': ylab, 'range': range_y}, hovermode='closest', annotations=annotations)
        fig = go.Figure(data=[trace], layout=layout)
    
    if plot_type=='interactive':
        plotly.offline.iplot(fig)
    else:
        py.image.ishow(fig)
        
robjects.r('''limma_voom <- function(rawcount_dataframe, design_dataframe, adjust="BH") {
    # Load packages
    suppressMessages(require(limma))
    suppressMessages(require(edgeR))
    # Convert design matrix
    design <- as.matrix(design_dataframe)

    # Create DGEList object
    dge <- DGEList(counts=rawcount_dataframe)

    # Filter genes
    keep <- filterByExpr(dge, design)
    dge <- dge[keep,]

    # Calculate normalization factors
    dge <- calcNormFactors(dge)
    # Run VOOM
    v <- voom(dge, plot=FALSE)
    # Fit linear model
    fit <- lmFit(v, design)
    # Make contrast matrix
    cont.matrix <- makeContrasts(de=B-A, levels=design)
    # Fit
    fit2 <- contrasts.fit(fit, cont.matrix)
    # Run DE
    fit2 <- eBayes(fit2)
    # Get results
    limma_dataframe <- topTable(fit2, adjust=adjust, number=nrow(rawcount_dataframe))
    
    # Return
    results <- list("limma_dataframe"= limma_dataframe, "rownames"=rownames(limma_dataframe))
    return (results)
}
''')

robjects.r('''limma <- function(rawcount_dataframe, design_dataframe, adjust="BH") {
    # Load packages
    suppressMessages(require(limma))
    suppressMessages(require(edgeR))
    # Convert design matrix
    design <- as.matrix(design_dataframe)

    # Create DGEList object
    dge <- DGEList(counts=rawcount_dataframe)

    # Filter genes
    keep <- filterByExpr(dge, design)
    dge <- dge[keep,]

    # Calculate normalization factors
    dge <- calcNormFactors(dge)

    # Fit linear model
    fit <- lmFit(normalizeMedianValues(dge$counts), design)  # Normalization so that each column has the same median value
    # Make contrast matrix
    cont.matrix <- makeContrasts(de=B-A, levels=design)
    # Fit
    fit2 <- contrasts.fit(fit, cont.matrix)
    # Run DE
    fit2 <- eBayes(fit2)
    # Get results
    limma_dataframe <- topTable(fit2, adjust=adjust, number=nrow(rawcount_dataframe))

    # Return
    results <- list("limma_dataframe"= limma_dataframe, "rownames"=rownames(limma_dataframe))
    return (results)
}
''')

robjects.r('''edgeR <- function(rawcount_dataframe, g1, g2) {
    # Load packages
    suppressMessages(require(limma))
    suppressMessages(require(edgeR))
    suppressMessages(require(dbplyr))

    colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
    rownames(colData) <- c(g1,g2)
    colnames(colData) <- c("group")
    colData$group = relevel(as.factor(colData$group), "Control")
    
    target <- colnames(rawcount_dataframe)
    DT <- colData %>% dplyr::slice(match(target, rownames(colData)))
    rownames(DT) <- target
    colData <- DT

    y <- DGEList(counts=rawcount_dataframe, group=colData$group)
    # Filter genes
    keep <- filterByExpr(y)
    y <- y[keep,]

    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y)
    res <- topTags(et, n=Inf)
    # Return
    res <- as.data.frame(res)
    results <- list("edgeR_dataframe"= res, "rownames"=rownames(res))
    return (results)
}


''')

robjects.r('''deseq2 <- function(rawcount_dataframe, g1, g2) {
    # Load packages
    suppressMessages(require(DESeq2))
    suppressMessages(require(edgeR))
    suppressMessages(require(dbplyr))
    colData <- as.data.frame(c(rep(c("Control"),length(g1)),rep(c("Condition"),length(g2))))
    rownames(colData) <- c(g1,g2)
    colnames(colData) <- c("group")
    colData$group = relevel(as.factor(colData$group), "Control")

    dge <- DGEList(counts=rawcount_dataframe, group=colData$group)
    # Filter genes
    keep <- filterByExpr(dge)
    dge <- dge[keep,]

    target <- colnames(rawcount_dataframe)
    DT <- colData %>% dplyr::slice(match(target, rownames(colData)))
    rownames(DT) <- target
    colData <- DT
    dds <- DESeqDataSetFromMatrix(countData = dge$counts +1 , colData = colData, design=~(group)) # add pseudocount =1

    dds <- DESeq(dds)
    res <- results(dds)
    
    res[which(is.na(res$padj)),] <- 1
    res <- as.data.frame(res)
    
    results <- list("DESeq_dataframe"= res, "rownames"=rownames(res))
    return(results)
}
''')

robjects.r('''preprocessing <- function(expr_matrix,  max_gene) {
    # normalize, most variable genes
    # Load packages
    suppressMessages(require(statmod))
    suppressMessages(require(limma))
    
    # http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html

    norm_matrix <- normalizeQuantiles(expr_matrix)
    
    if (nrow(expr_matrix) > max_gene) {
        winsorize <- function (x, fraction=0.05) {
            if(length(fraction) != 1 || fraction < 0 || fraction > 0.5) {
                stop("bad value for 'fraction'")
            }
            lim <- quantile(x, probs=c(fraction, 1-fraction), na.rm = TRUE,)
            x[ x < lim[1] ] <- lim[1]
            x[ x > lim[2] ] <- lim[2]
            x
        }
        # winsorize to remove 2 most extreme cells (from each side)
        wed <- t(apply(norm_matrix, 1, winsorize, fraction=2/ncol(norm_matrix)))

        means <- rowMeans(wed, na.rm = TRUE)
        vars <- apply(wed, 1, var, na.rm = TRUE)
        cv2 <- vars/means^2
        minMeanForFit <- unname( quantile( means[ which( cv2 > .3 ) ], .95 ) )
        useForFit <- means >= minMeanForFit # & spikeins
        fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] )
        a0 <- unname( fit$coefficients["a0"] )
        a1 <- unname( fit$coefficients["a1tilde"])
        afit <- a1/means+a0
        varFitRatio <- vars/(afit*means^2)
        varorder <- order(varFitRatio,decreasing=T)
        return(as.data.frame(norm_matrix[varorder[0:max_gene],]))
    } else {
        return(as.data.frame(norm_matrix))
    }
}
''')

def get_signatures(classes, dataset, normalization, method, meta_class_column_name, meta_id_column_name, filter_genes, logData, pseudocount):
    expr_df = dataset['rawdata']
    if filter_genes == True:
        expr_df = dataset['rawdata+filter_genes']
        
    signatures = dict()

    for cls1, cls2 in combinations(classes, 2):
        cls1_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"][meta_class_column_name]==cls1, meta_id_column_name].tolist() #control
        cls2_sample_ids = dataset["dataset_metadata"].loc[dataset["dataset_metadata"][meta_class_column_name]==cls2, meta_id_column_name].tolist() #case
        
        signature_label = " vs. ".join([cls2, cls1])
        print (signature_label)

        # check values all non-negative for limma_voom, edgeR and DESeq2
        if method == "limma_voom" or method == "edgeR" or method == "DESeq2":
            if logData: # transform back to non-log data
                unlog_expr_df = np.exp2(expr_df) - pseudocount
                print(f"Info. Log transformed data. Base 2 exponentiation is applied")
                if (unlog_expr_df < 0).any().any(): # requires non-negative values
                    unlog_expr_df = np.exp2(expr_df) # the entire matrix is shifted by pseudocount
                    print(f"Warning! {method} requires all non-negative values. Negative values detected, only apply base2 exponentiation to log data, pseudocount ignored.")
                expr_df = unlog_expr_df

            if (expr_df < 0).any().any():
                print(f"Error! {method} requires non-negative gene expression values. Negative values detected")
                sleep(10)
                raise SystemExit(1)

        # check values all integer for DESeq2
        if method == "DESeq2":
            if (expr_df.fillna(-999) % 1  == 0).all().all():
                pass
            else:
                expr_df = expr_df.round() # requires all non-negative integer counts
                print(f"Warning! {method} requires all non-negative integer count values. Non integer values detected, round to integer values.")

        if method == "limma_voom":
            limma_voom = robjects.r['limma_voom']
            design_dataframe = pd.DataFrame([{'index': x, 'A': int(x in cls1_sample_ids), 'B': int(x in cls2_sample_ids)} for x in expr_df.columns]).set_index('index')
            processed_data = {"expression": expr_df, 'design': design_dataframe}
            limma_results = pandas2ri.conversion.rpy2py(limma_voom(pandas2ri.conversion.py2rpy(processed_data['expression']), pandas2ri.conversion.py2rpy(processed_data['design'])))
            
            # signature = pd.DataFrame(limma_results[0])
            # signature.index = limma_results[1]
            signature = pd.DataFrame(limma_results[0], index = limma_results[0].colnames, columns = limma_results[1]).T
            signature = signature.sort_values("t", ascending=False)

        elif method == "characteristic_direction":
            max_gene = 25 * 1000
            if len(expr_df) > max_gene:
                max_gene = int(len(expr_df) * 0.50) # top 50% most variable genes
                print (f'Using the most variable {max_gene} genes')
            preprocessing = robjects.r['preprocessing']
            processed_data = pandas2ri.conversion.rpy2py(preprocessing(pandas2ri.conversion.py2rpy(expr_df), pandas2ri.conversion.py2rpy(max_gene)))
            signature = characteristic_direction(processed_data.loc[:, cls1_sample_ids], processed_data.loc[:, cls2_sample_ids], calculate_sig=True)
            signature = signature.sort_values("CD-coefficient", ascending=False)

        elif method == "edgeR":
            edgeR = robjects.r['edgeR']
            edgeR_results = pandas2ri.conversion.rpy2py(edgeR(pandas2ri.conversion.py2rpy(expr_df), pandas2ri.conversion.py2rpy(cls1_sample_ids), pandas2ri.conversion.py2rpy(cls2_sample_ids)))
            
            signature = pd.DataFrame(edgeR_results[0])
            signature.index = edgeR_results[1]
            signature = signature.sort_values("logFC", ascending=False)

        elif method == "DESeq2":
            DESeq2 = robjects.r['deseq2']
            DESeq2_results = pandas2ri.conversion.rpy2py(DESeq2(pandas2ri.conversion.py2rpy(expr_df), pandas2ri.conversion.py2rpy(cls1_sample_ids), pandas2ri.conversion.py2rpy(cls2_sample_ids)))
            
            signature = pd.DataFrame(DESeq2_results[0])
            signature.index = DESeq2_results[1]
            signature = signature.sort_values("log2FoldChange", ascending=False)
                        
        signatures[signature_label] = signature

    return signatures

def run_volcano(signature, signature_label, dataset, pvalue_threshold, logfc_threshold, plot_type):
    color = []
    text = []
    for index, rowData in signature.iterrows():
        if "AveExpr" in rowData.index: # limma limma_voom
            expr_colname = "AveExpr"
            pval_colname = "P.Value"
            logfc_colname = "logFC"
        elif "logCPM" in rowData.index: #edgeR
            expr_colname = "logCPM"
            pval_colname = "PValue"
            logfc_colname = "logFC"
        elif "baseMean" in rowData.index: #DESeq2
            expr_colname = "baseMean"
            pval_colname = "pvalue"
            logfc_colname = "log2FoldChange"
        # Text
        text.append('<b>'+index+'</b><br>Avg Expression = '+str(round(rowData[expr_colname], ndigits=2))+'<br>logFC = '+str(round(rowData[logfc_colname], ndigits=2))+'<br>p = '+'{:.2e}'.format(rowData[pval_colname])+'<br>FDR = '+'{:.2e}'.format(rowData[pval_colname]))

        # Color
        if rowData[pval_colname] < pvalue_threshold:
            if rowData[logfc_colname] < -logfc_threshold:
                color.append('blue')
            elif rowData[logfc_colname] > logfc_threshold:
                color.append('red')
            else:
                color.append('grey')
        else:
            color.append('grey')

    volcano_plot_results = {'x': signature[logfc_colname], 'y': -np.log10(signature[pval_colname]), 'text':text, 'color': color, 'signature_label': signature_label, 'plot_type': plot_type}
    return volcano_plot_results
        

def plot_volcano(volcano_plot_results):
    spacer = ' '*50
    plot_2D_scatter(
        x=volcano_plot_results['x'],
        y=volcano_plot_results['y'],
        text=volcano_plot_results['text'],
        color=volcano_plot_results['color'],
        symmetric_x=True,
        xlab='log2FC',
        ylab='-log10P',
        title='<b>{volcano_plot_results[signature_label]} Signature | Volcano Plot</b>'.format(**locals()),
        labels=volcano_plot_results['signature_label'].split(' vs. '),
        plot_type=volcano_plot_results['plot_type'],
        de_type='volcano'
    )        

def run_maplot(signature, signature_label='', pvalue_threshold=0.05, logfc_threshold=1.5, plot_type='interactive'):

    # Loop through signature
    color = []
    text = []
    for index, rowData in signature.iterrows():
        if "AveExpr" in rowData.index: # limma
            expr_colname = "AveExpr"
            pval_colname = "adj.P.Val"
            logfc_colname = "logFC"
        elif "logCPM" in rowData.index: #edgeR
            expr_colname = "logCPM"
            pval_colname = "PValue"
            logfc_colname = "logFC"
        elif "baseMean" in rowData.index: #DESeq2
            expr_colname = "baseMean"
            pval_colname = "padj"
            logfc_colname = "log2FoldChange"
        # Text
        text.append('<b>'+index+'</b><br>Avg Expression = '+str(round(rowData[expr_colname], ndigits=2))+'<br>logFC = '+str(round(rowData[logfc_colname], ndigits=2))+'<br>p = '+'{:.2e}'.format(rowData[pval_colname])+'<br>FDR = '+'{:.2e}'.format(rowData[pval_colname]))

        # Color
        if rowData[pval_colname] < pvalue_threshold:
            if rowData[logfc_colname] < -logfc_threshold:
                color.append('blue')
            elif rowData[logfc_colname] > logfc_threshold:
                color.append('red')
            else:
                color.append('black')
        else:
            color.append('black')
    
    # Return 
    volcano_plot_results = {'x': signature[expr_colname], 'y': signature[logfc_colname], 'text':text, 'color': color, 'signature_label': signature_label, 'plot_type': plot_type}
    return volcano_plot_results

def plot_maplot(volcano_plot_results):    
    plot_2D_scatter(
        x=volcano_plot_results['x'],
        y=volcano_plot_results['y'],
        text=volcano_plot_results['text'],
        color=volcano_plot_results['color'],
        symmetric_y=True,
        xlab='Average Expression',
        ylab='logFC',
        title='<b>{volcano_plot_results[signature_label]} Signature | MA Plot</b>'.format(**locals()),
        labels=volcano_plot_results['signature_label'].split(' vs. '),
        plot_type=volcano_plot_results['plot_type']
    )
        
def xenaFileDownloadLink(host, dataset_name):
    if dataset_name == None:
        return None
    if host.endswith('/'):
        host = host[:-1]
    head, tail = os.path.split(dataset_name)
    tail = urllib.parse.quote_plus(tail)
    return host + '/download/' + os.path.join(head, tail)

def reformat_gmt(file):
    geneset_dict = {} 
    with open(file) as gmt:
        for line in gmt: 
            gene = line.split()
            del gene[1] 
            geneset_dict[gene[0]] = [gene[1]] 
            for item in range(2, len(gene)):
                geneset_dict[gene[0]].append(gene[item])
    return geneset_dict

def reformat_signature(signatures, diff_gex_method):
    for label, signature in signatures.items():
        if diff_gex_method == "characteristic_direction":
            signature_df = signature[['CD-coefficient']].copy()
            signature_df.reset_index(inplace=True)
            signature_df.rename(columns={'index': '0', 'CD-coefficient': '1'}, inplace=True)
        elif diff_gex_method == "limma" or diff_gex_method == "limma_voom" or diff_gex_method == "edgeR":
            signature_df = signature[['logFC']].copy()
            signature_df.reset_index(inplace=True)
            signature_df.rename(columns={'index': '0', 'logFC': '1'}, inplace=True)
        elif diff_gex_method == "DESeq2":
            signature_df = signature[['log2FoldChange']].copy()
            signature_df.reset_index(inplace=True)
            signature_df.rename(columns={'index': '0', 'log2FoldChange': '1'}, inplace=True)
    
    return signature_df
    
def run_blitzGSEA(signature, library):
    result = blitz.gsea(signature, library)
    return result
