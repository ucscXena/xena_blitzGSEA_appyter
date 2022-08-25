def top_table(app, url_prefix = '/top_table', DATA_DIR=''):
    import dash
    from dash import dash_table, html, dcc, Input, Output
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import pandas as pd
    import numpy as np

    app = dash.Dash(
        'top_table',
        server = app,
        routes_pathname_prefix=url_prefix + '/'
    )
    app.layout = html.Div(children=[html.H1(children='Top Table'),
        html.Div(children='''Produced with Dash: A web application framework for your data.'''), # more info about reading GSEA plots 
        html.Div([
            html.Div(
                dcc.Graph(id='top-table',
                        config={'editable': True, 'edits': {'annotationPosition': False, 'axisTitleText': False, 'titleText': False}},
                        figure=top_table
                        ),
                style={'display': 'inline-block'}
            ),
            html.Div(
                dcc.Graph(id='running-plot'),
                style={'display': 'inline-block'}
            ),
            html.Div(
                dash_table.DataTable(id = 'datatable-interactivity',
                                    columns = [{"name": i, "id": i, "deletable": False, "selectable": False} for i in df.columns],
                                    data = df.to_dict('records'),
                                    editable = True,
                                    sort_action = "native",
                                    sort_mode = "multi",
                                    row_deletable = False,
                                    selected_columns = [],
                                    selected_rows = [],
                                    page_action = "native",
                                    page_current= 0,
                                    page_size= 15,
                                    ),
                style={'display': 'inline-block'}
            ), 
        ]),
    ])
    @app.callback(
        Output('running-plot', 'figure'),
        Output('datatable-interactivity', 'data'),
        [Input('top-table', 'relayoutData')]
    )

    def update_output(data):
        signature = sig
        library = lib
        result = res
        
        if data is None: 
            geneset = result.index[0]
        else:
            key = list(data)
            key = key[0]
            geneset = data.get(key)
            
        signature_map = {}
        details_dict = {}
        signature = signature.sort_values(1, ascending=False).set_index(0)
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
        rank_vec = signature[1]
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
            value = signature.loc[key][1] 
            details_dict[key].append(value)
            details_dict[key].append(running_sum[details_dict[key][0]])
            if key in leading_edges:
                details_dict[key].append("Yes")
            else:
                details_dict[key].append("No")
                
        table = pd.DataFrame.from_dict(details_dict, orient = 'index')
        table.reset_index(inplace=True)
        table.columns = ["Gene", "Rank in Gene List", "Rank Metric Score", "Running ES", "Leading Edge"]  
        
        return fig, table.to_dict('records')