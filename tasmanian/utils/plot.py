import pandas as pd
import numpy as np
import plotly
from plotly.offline import plot
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import sys


def plot_html(table):
    '''
        function splits table (a dictionary of 3 dataframes as values) and plots them in 
        a html report.
        INPUT: table = dictionary with keys: 'intersections','complement','non_intersection'
    '''

    def table_to_json(df, normalize=True):
        #df = table['intersection'] #pd.read_csv(sys.argv[1])
        
        a = ['a_a','a_t','a_c','a_g']
        c = ['c_a','c_t','c_c','c_g']
        g = ['g_a','g_t','g_c','g_g']
        t = ['t_a','t_t','t_c','t_g']
        mutations = set(a+c+g+t).difference(set(['a_a','c_c','t_t','g_g']))

        df.iloc[:,2:] #/= df[eval(column_name[0])].sum(axis=1).values.reshape(-1,1)
        df = df.replace(np.inf,0).fillna(0)


        # normalize the data to later rescale the data to between 0-1
        dfc = df.copy()
        errors = {'a': a, 'c': c, 'g': g, 't': t}
        
        if normalize:
            for i in errors.keys():
                dfc.loc[:, errors[i]] /= (dfc.loc[:, errors[i]].sum(axis=1).values.reshape(-1,1) + 10e-6 )
            
        table_max = dfc[mutations].max().max() * 1.1
        table_min = dfc[mutations].min().min() * 0.5

        df1 = dfc[(dfc.read==1)]
        df2 = dfc[(dfc.read==2)]
        

        # plot the data in horizontal subplots
        fig = make_subplots(rows=1, cols=2,
                            subplot_titles=("Read 1", "Read 2"))
         
        for mut in mutations:

            fig.add_trace(
                go.Scatter(
                    x=np.arange(df1.shape[0]),
                    y=df1[mut],
                    mode='markers',
                    name= mut + ' - R1'
                ),
                row=1, col=1
            )

            fig.add_trace(
                go.Scatter(
                    x = np.arange(df2.shape[0]),
                    y = df2[mut],
                    mode = 'markers',
                    name = mut + ' - R2',
                ),
                row=1, col=2
            )
            
        # Update markers and color selection
        markers = [
                   "circle", "circle-open-dot", "square", "diamond", "cross", "triangle-up", 
                   "triangle-down","hexagon", "hexagram", "star", "octagon", "square-x"
                  ]
        colors = [
                   'rgb(31, 119, 180)', 'rgb(255, 127, 14)', 'rgb(44, 160, 44)', 'rgb(214, 39, 40)',
                  'rgb(148, 103, 189)', 'rgb(140, 86, 75)', 'rgb(227, 119, 194)', 'rgb(127, 127, 127)',
                  'rgb(188, 189, 34)', 'rgb(23, 190, 207)','rgb(188, 189, 150)', 'rgb(23, 190, 150)'
                 ]
        for n in np.arange(0,len(fig['data']),2):
            m = n//2 + n%2
            fig['data'][n]['marker']['symbol']=markers[m] 
            fig['data'][n+1]['marker']['symbol']=markers[m]

            fig['data'][n]['marker']['color']=colors[m]
            fig['data'][n+1]['marker']['color']=colors[m]


        # final aestetic touches
        fig.update_traces(mode='markers', marker_line_width=2, marker_size=10) #, visible="legendonly")
        fig.update_layout(height=800, xaxis_title="Position in read", yaxis_title="Normalized counts")

        fig.update_yaxes(range=[table_min, table_max], row=1, col=1)
        fig.update_yaxes(range=[table_min, table_max], row=1, col=2)

        # convert it to JSON
        fig_json = fig.to_json()

        return fig_json


    fig_json_i = table_to_json(table['intersection'], normalize=False)
    fig_json_c = table_to_json(table['complement'], normalize=False)
    fig_json_n = table_to_json(table['non_intersection'], normalize=False)

    fig_json_i_norm = table_to_json(table['intersection'])
    fig_json_c_norm = table_to_json(table['complement'])
    fig_json_n_norm = table_to_json(table['non_intersection'])


    # HTML template
    template = """<html>
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <title>Tasmanian Report</title>

        <style>
            .button {{
              background-color: white; /* Green */
              border: 2px solid #8bd9ff;
              color: black;
              padding: 16px 32px 0px -10px;
              text-align: center;
              text-decoration: none;
              display: inline-block;
              font-size: 16px;
              margin: 20px 2px 0px 100px;
              transition-duration: 0.4s;
              cursor: pointer;
            }}

            .button:hover {{
              background-color: #8bd9ff;
              color: white;
            }}
            h1{{ 
                background-color: #8bd9ff;
                padding: 30px 0px 30px 0px;
                margin: 10px 120px 0px 80px;
            }}
        </style>

    </head>
    <body>
        <script type="text/javascript">
           function toggle_visibility(id1, id2, id3, id4, id5, id6) {{
               var e = document.getElementById(id1);
               var f = document.getElementById(id2);
               var g = document.getElementById(id3);
               var h = document.getElementById(id4);
               var i = document.getElementById(id5);
               var j = document.getElementById(id6);

               if(e.style.display == 'block') {{
                  e.style.display = 'none';
                  f.style.display = 'block';
                  g.style.display = 'none';
                  h.style.display = 'block';
                  i.style.display = 'none';
                  j.style.display = 'block';
                }}
               else {{
                  e.style.display = 'block';
                  f.style.display = 'none';
                  g.style.display = 'block';
                  h.style.display = 'none';
                  i.style.display = 'block';
                  j.style.display = 'none';
                }}
           }}
        </script>

        <h1 align="center">Tasmanian artifacts metrics results </h1>

        <button class="button button1"; onclick="toggle_visibility('divPlotly1_norm', 'divPlotly1', 'divPlotly2_norm', 'divPlotly2', 'divPlotly3_norm', 'divPlotly3');">Counts/Normalize Counts</button>
    
        <h2 style="padding-left: 40px; padding-top: 90px;">Intersections</h2>
        <h3 style="padding-left: 40px; padding-right: 800; ">Includes all bases that intersect some fragment provided in the bed-file</h3>



        <!-- plot 1 -->

        <div id='divPlotly1_norm'>
        <script>
            var plotly_data = {}
            Plotly.react('divPlotly1_norm', plotly_data.data, plotly_data.layout);
        </script>
        </div>

        <div id='divPlotly1'>
        <script>
            var plotly_data = {}
            Plotly.react('divPlotly1', plotly_data.data, plotly_data.layout);
        </script>
        </div>

        <h2 style="padding-left: 40px;padding-top: 60px;">Complements</h2>
        <h3 style="padding-left: 40px; padding-right: 800;">Includes all bases from that do not intersect a fragment, from reads that intersect a fragment provided in the bed-file</h3>

        <!-- plot 2 -->

        <div id='divPlotly2_norm'>
        <script>
            var plotly_data = {}
            Plotly.react('divPlotly2_norm', plotly_data.data, plotly_data.layout);
        </script>
        </div>

        <div id='divPlotly2'>
        <script>
            var plotly_data = {}
            Plotly.react('divPlotly2', plotly_data.data, plotly_data.layout);
        </script>
        </div>

        <h2 style="padding-left: 40px;padding-top: 60px;">Non-intersections</h2>
        <h3 style="padding-left: 40px; padding-right: 800;">Includes all bases from reads with no intersections with the bed-file</h3>

        <!-- plot 3 -->

        <div id='divPlotly3_norm'>
        <script>
            var plotly_data = {}
            Plotly.react('divPlotly3_norm', plotly_data.data, plotly_data.layout);
        </script>
        </div>

        <div id='divPlotly3'>
        <script>
            var plotly_data = {}
            Plotly.react('divPlotly3', plotly_data.data, plotly_data.layout);
        </script>
        </div>
        
        <script> 
            divPlotly1_norm.style.display = 'none'; 
            divPlotly2_norm.style.display = 'none'; 
            divPlotly3_norm.style.display = 'none'; 
        </script>

    </body>

    </html>"""

    # write the JSON to the HTML template
    #with open('Tasmanian_artifact_metrics_report.html', 'w') as f:
    #   f.write(template.format(fig_json))
    return template.format(fig_json_i_norm, fig_json_i, fig_json_c_norm, fig_json_c, fig_json_n_norm, fig_json_n)



#if __name__=='__main__':
#    
#    normalize = False
#    for n,i in enumerate(sys.argv):
#        if i in ["-normalize","--normalize","-n","--n","--norm","-norm"]:
#            normalize=True
#        if i in ["--table", "-t","--t","-table"]:
#            
