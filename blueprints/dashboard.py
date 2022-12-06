def dashboard(app, url_prefix='/dashboard', **kwargs):
    import os
    import pandas as pd
    import plotly
    from plotly import tools
    from plotly.subplots import make_subplots
    import plotly.express as px
    import plotly.graph_objs as go
    import numpy as np
    import dash
    from dash import Dash, dash_table, dcc, html
    from dash.dependencies import Input, Output

    import blitzgsea as blitz

    from functools import lru_cache

    external_stylesheets = ["https://fonts.googleapis.com/css2?family=Nunito:wght@300;400&display=swap"]

    dashapp = dash.Dash(
        'dashboard',
        server=app,
        external_stylesheets=external_stylesheets,
        requests_pathname_prefix='/'.join([app.config['PREFIX'].rstrip('/'), url_prefix.strip('/'), '']),
        routes_pathname_prefix='/'.join(['', url_prefix.strip('/'), '']),
     )

    @lru_cache()
    def session_from_href(href): 
        '''retrieve session'''
        return href.replace(url_prefix, '')

    @lru_cache()
    def signature_for_session(session): 
        '''retrieve signature file for associated session'''
        return pd.read_csv(
        os.path.join(session, 'signature.tsv'),
        sep='\t',
        )

    @lru_cache()
    def full_table_for_session(session):
        '''retrieve full table file for session'''
        return pd.read_csv(
        os.path.join(session, 'full_table.tsv'),
        sep='\t'
        ).set_index('Term')

    def reformat_gmt(file, session):
        '''reformat library file for session'''
        geneset_dict = {} 
        path = os.getcwd()
        session = session.split('/')
        path = path + '/data/output/' + session[3] + '/' + file
        with open(path) as gmt:
            for line in gmt: 
                gene = line.split()
                del gene[1] 
                geneset_dict[gene[0].strip('\"')] = [gene[1]] 
                for item in range(2, len(gene)):
                    geneset_dict[gene[0].strip('\"')].append(gene[item])

            return geneset_dict
    
    def get_n_table(file, session):
        '''get n for gene sets requested'''
        n = 0
        path = os.getcwd()
        session = session.split('/')
        path = path + '/data/output/' + session[3] + '/' + file
        with open(path) as file:
            for line in file:
                n = int(line.strip())

        return n

    def top_table_figure_for_session(session):
        '''plot top table'''
        signature = signature_for_session(session)
        library = reformat_gmt('library.tsv', session)
        result = full_table_for_session(session) 
        n = get_n_table('n_value.tsv', session)

        lines = []
        sig = signature.sort_values('1', ascending=False).set_index('0')
        sig = sig[~sig.index.duplicated(keep='first')]
        
        fig = go.Figure(layout_xaxis_range=[-0.1,2],
                        layout_yaxis_range=[-0.2, 1.1],
                        layout = {'xaxis': {'title': 'x-label',
                            'visible': False,
                            'showticklabels': False},
                            'yaxis': {'title': 'y-label',
                            'visible': False,
                            'showticklabels': False}}
                    )
        fig.add_shape(type="line", x0=0.2, y0=-0.1, x1=0.2, y1=1, line_width=2, line_dash="solid")
        fig.add_shape(type="line", x0=0.8, y0=-0.1, x1=0.8, y1=1, line_width=2, line_dash="solid")
        fig.add_shape(type="line", x0=0.9, y0=-0.1, x1=0.9, y1=1, line_width=2, line_dash="solid")
        
        ln = np.linspace(-0.1,1,n+1)[::-1]
        for line in range(0, len(ln)):
            fig.add_shape(type="line", x0=0, y0=ln[line], x1=1.5, y1=ln[line], line_width=2, line_dash="solid")
        fig.add_annotation(x=0.03, 
                        y=1.04,
                        showarrow=False, 
                        xanchor="left",
                        captureevents= False,
                        text = "NES",
                        font=dict(size=16)
                        )
        fig.add_annotation(x=0.80, 
                        y=1.04,
                        showarrow=False, 
                        xanchor="left",
                        captureevents= False,
                        text = "SIZE",
                        font=dict(size=16)
                        )
        fig.add_annotation(x=0.92, 
                        y=1.04,
                        showarrow=False, 
                        xanchor="left",
                        captureevents= False,
                        text="SET",
                        font=dict(size=16)
                        )
        
        for i in range(n):
            fig.add_annotation(x=0.03, 
                            y=(ln[i]+ln[i+1])/2, 
                            showarrow=False,
                            xanchor="left",
                            captureevents= False,
                            text="{:.3f}".format(result.iloc[i, 1])
                            )
            fig.add_annotation(x=0.82, 
                            y=(ln[i]+ln[i+1])/2, 
                            showarrow=False,
                            captureevents= False,
                            xanchor="left",
                            text="{}".format(result.iloc[i, 5])
                            )
            fig.add_annotation(x=0.92, 
                            y=(ln[i]+ln[i+1])/2,
                            showarrow=False,
                            xanchor="left",
                            captureevents= True,
                            text="{}".format(result.index[i]),
                            bgcolor = '#f2f3f4'
                            )
            gs = set(library[result.index[i]])
            hits = np.array([i for i,x in enumerate(sig.index) if x in gs])
            hits = (hits/len(sig.index))*0.6+0.2
            if result.iloc[i, 1] > 0:
                for hit in hits:
                    lines.append(dict(type="line", x0=hit, y0=ln[i+1],
                                    x1=hit, y1=ln[i],
                                    line_width=0.25,
                                    line_dash="solid",
                                    line=dict(color="red"))
                                )
            else:
                for hit in hits:
                    lines.append(dict(type="line", x0=hit, y0=ln[i+1],
                                    x1=hit, y1=ln[i],
                                    line_width=0.25,
                                    line_dash="solid",
                                    line=dict(color="blue"))
                                )
                        
        fig.update_layout(shapes=lines, plot_bgcolor='white', width = 1000, height = (0.55*100*n))
                                
        return fig

    def create_plot(session, data):
        '''plot running sum plot for gene set'''
        signature = signature_for_session(session)
        library = reformat_gmt('library.tsv', session)
        result = full_table_for_session(session)

        if data == None:
            geneset = result.index[0]
        else:
            geneset = data['annotation']['text']
            
        signature_map = {}
        details_dict = {}
        signature = signature.sort_values('1', ascending=False).set_index('0')
        signature = signature[~signature.index.duplicated(keep='first')]
        leading_edges = result.loc[geneset, 'leading_edge']
        leading_edges = leading_edges.split(',')
        
        for i,h in enumerate(signature.index):
            signature_map[h] = i
        gs = set(library[geneset])
        geneset_list = library[geneset]
        hits = [i for i, x in enumerate(signature.index) if x in gs]
        running_sum, es = blitz.enrichment_score(np.array(np.abs(signature.iloc[:,0])), signature_map, gs)
        running_sum = list(running_sum)

        fig = make_subplots(rows=3, 
                            cols=1,
                            vertical_spacing=0,
                            row_heights=[2,0.5,1],
                            subplot_titles=(geneset, "")
                        )

        # FIRST SUBPLOT - RS
        nn = np.where(np.abs(running_sum)==np.max(np.abs(running_sum)))[0][0]
        fig.add_trace(go.Scatter(y=list(running_sum),
                                mode='lines',
                                line=dict(color='lightgreen', width=3),
                                showlegend=False
                                ),
                    row=1,
                    col=1
                    )
        if es > 0:
            fig.add_shape(go.layout.Shape(type="line", x0=nn, y0=np.min(running_sum),
                                        x1=nn, y1=np.max(running_sum), 
                                        line_width=1.5,
                                        line_dash="dash",
                                        line=dict(color="red"),), 
                        row=1, col=1
                        )
            fig.add_annotation(x=len(running_sum)/30, y=0.1,
                            showarrow=False,
                            xanchor="left",
                            text="NES="+"{:.3f}".format(result.loc[geneset,"nes"]),
                            font=dict(size=16, color="black"),
                            bgcolor='white'
                            )
        else:
            fig.add_shape(go.layout.Shape(type="line", x0=nn, y0=np.min(running_sum),
                                        x1=nn, y1=np.max(running_sum), 
                                        line_width=1.5,
                                        line_dash="dash",
                                        line=dict(color="red"),), 
                        row=1, col=1
                        )
            fig.add_annotation(x=len(running_sum)/1.5, y=-0.1,
                            showarrow=False,
                            xanchor="left",
                            text="NES="+"{:.3f}".format(result.loc[geneset,"nes"]),
                            font=dict(size=16, color="black"),
                            bgcolor='white'
                            )
            
        # SECOND SUBPLOT - HITS
        lines = []
        for hit in hits:
            lines.append(dict(type="line", x0=hit, y0=-1, x1=hit, y1=1,
                            line_width=0.50,
                            line_dash="solid", 
                            line=dict(color="black"), 
                            xref = "x2",
                            yref = "y2"
                            )
                        )
        fig.update_layout(shapes=lines)

        # THIRD SUBPLOT - RANK
        rank_vec = signature['1']
        x = np.arange(0.0, len(rank_vec), 20).astype("int")
        x = np.append(x, signature.shape[0]-1)
        fig.add_trace(go.Scatter(x=x, y=np.array(rank_vec)[x],
                                mode='lines',
                                fill='tonexty',
                                fillcolor='lightgrey',
                                line=dict(color='black', width=1.25),
                                showlegend=False
                                ), row=3, col=1
                    )
        fig.add_shape(go.layout.Shape(type="line", x0=0, y0=0, 
                                    x1=len(rank_vec), y1=0,
                                    line_width=1, 
                                    line_dash="solid",
                                    ), row=3, col=1
                    )
        minabs = np.min(np.abs(rank_vec))
        zero_cross = int(np.where(np.abs(rank_vec)==minabs)[0][0])
        fig.add_shape(go.layout.Shape(type="line", x0=zero_cross, y0=np.min(rank_vec), x1=zero_cross, y1=np.max(rank_vec), line_width=1.5, line_dash="dash", line=dict(color="blue"),), row=3, col=1)
        fig.add_annotation(x=zero_cross, y=np.max(rank_vec)/3,
                        xref="x3", yref="y3",
                        showarrow=False,
                        xanchor="center",
                        text="Zero crosses at "+str(zero_cross),
                        font=dict(size=12, color="black"),
                        bgcolor='white'
                        )
        fig.update_xaxes(range=[0, len(running_sum)],
                        showticklabels=False,
                        showline = True,
                        linecolor = 'black',
                        linewidth = 1,
                        row=1, col=1, 
                        mirror=True
                        )
        fig.update_yaxes(title_text="Enrichment Score (ES)",
                        range=[np.min(running_sum), (np.max(running_sum))],
                        showline = True,
                        linecolor = 'black',
                        linewidth = 1,
                        row=1, col=1,
                        mirror=True
                        )
        fig.update_xaxes(range=[0, len(running_sum)],
                        showticklabels=False,
                        row=2,
                        col=1,
                        )
        fig.update_yaxes(range=[-1, 1], 
                        showticklabels=False,
                        showline = True,
                        linecolor = 'black',
                        linewidth = 1,
                        row=2, col=1,
                        mirror=True
                        )
        fig.update_xaxes(title_text="Rank in Ordered Dataset",
                        range=[0, len(running_sum)],
                        showticklabels=True,
                        tickformat='000',
                        minor_ticks='outside',
                        showline = True,
                        linecolor = 'black',
                        linewidth = 1,
                        row=3, col=1,
                        mirror=True
                        )
        fig.update_yaxes(title_text="Ranked list metric",
                        range=[np.min(rank_vec),np.max(rank_vec)],
                        showticklabels=True,
                        showline = True,
                        linecolor = 'black',
                        linewidth = 1,
                        row=3, col=1,
                        mirror=True
                        )
        fig.update_layout(plot_bgcolor='white', height=600, width=750)
        
        for i, x in enumerate(signature.index):
            if x in geneset_list:
                details_dict[x] = [i]
        for key in details_dict.keys():
            value = signature.loc[key]['1'] 
            details_dict[key].append(value)
            details_dict[key].append(running_sum[details_dict[key][0]])
            if key in leading_edges:
                details_dict[key].append("Yes")
            else:
                details_dict[key].append("No")
                
        table = pd.DataFrame.from_dict(details_dict, orient = 'index')
        table.reset_index(inplace=True)
        table.columns = ['gene', 'rank', 'rank_metric', 'running_es', "leading_edge"]

        return fig, table.to_dict('records'), table.to_dict('records'), geneset

    dashapp.layout = html.Div(children=[
        dcc.Location(id='url', refresh=False),
        dcc.Store(id='detailed-df'),
        dcc.Store(id='current-geneset'),

        html.Div(children=[
            html.H1(children='Gene Set Data Dashboard'),
        ], style={'fontFamily': 'Nunito', 'textAlign': 'center'}),

        html.Div(children=[
            html.H2(children='Top Table:'),
            html.P(children='The top table displays the top n enriched gene sets with its normalized enrichment score (NES) and the distribution of hits relative to the gene ranking of the signature. Please click a gene set under \'SET\' in the top table to produce its running sum plot and detailed output table. While the text is enlarged, data is being computed.'),
            dcc.Loading(id='loading-top-table', type='circle', children=[dcc.Graph(id='top-table', config={'editable': False, 'edits': {'annotationPosition': False, 'axisTitleText': False, 'titleText': False}, 'showAxisDragHandles': False, 'toImageButtonOptions': {'filename': 'top_table'}})])
        ], style={'fontFamily': 'Nunito', 'textAlign':'left', 'marginLeft': '235px', 'marginRight': '235px'}),

        html.Div(id='test'),
        
        html.Div(children=[
            html.H2(children='Running Sum Plot:'),
            dcc.Loading(type='circle', children=[dcc.Graph(id='running-plot',  config={'toImageButtonOptions': {'filename': 'running_sum_plot'}})]),
        ], style={'fontFamily': 'Nunito', 'textAlign':'left', 'marginLeft': '235px', 'marginRight': '235px'}),

        html.Div(children=[
            html.H2(children='Detailed Output Table:'),
            dash_table.DataTable(
                id='detailed-table',
                columns=[
                    {'name': 'gene', 'id': 'gene'},
                    {'name': 'rank', 'id': 'rank'},
                    {'name': 'rank_metric', 'id': 'rank_metric'},
                    {'name': 'running_es', 'id': 'running_es'},
                    {'name': 'leading_edge', 'id': 'leading_edge'}
                ],
                editable=False,
                style_cell={'textAlign': 'left'},
                sort_action='native',
                sort_mode='multi',
                row_deletable=False,
                selected_columns=[],
                selected_rows=[],
                page_action='native',
                page_current=0,
                page_size=10),
            html.Button("Download Detailed Output", id="btn-download-detail"),
            dcc.Download(id='download-detailed-output')
        ], style={'fontFamily': 'Nunito', 'textAlign':'left', 'marginLeft': '235px', 'marginRight': '235px'})
    ])

    @dashapp.callback(
        Output('top-table', 'figure'), 
        Input('url', 'href')
    )
    def plot_top_table(href):
        '''plot top table for session'''
        session = session_from_href(href)

        return top_table_figure_for_session(session)

    @dashapp.callback(
        [Output('running-plot', 'figure'),
        Output('detailed-table', 'data'),
        Output('detailed-df', 'data'), # store detailed output for gene set
        Output('current-geneset', 'data')], # store gene set 
        [Input('url', 'href'),
        Input('top-table', 'clickAnnotationData')]
    )

    def update_plot(href, data):
        '''take in data from annotation click on top table and plot running sum plot for gene set'''
        session = session_from_href(href)

        return create_plot(session, data)

    @dashapp.callback(
        Output('download-detailed-output', 'data'), 
        Input('btn-download-detail', 'n_clicks'),
        Input('detailed-df', 'data'),
        Input('current-geneset', 'data')
    )

    def download_data(n_clicks, df, geneset):
        '''create downloadable .csv file for gene set'''
        if n_clicks == None:

            return dash.no_update
        else:
            df = pd.DataFrame(df)

            return dcc.send_data_frame(df.to_csv, '{0}_detailed_output.csv'.format(geneset))