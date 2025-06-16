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
        mutations = list( set(a+c+g+t).difference(set(['a_a','c_c','t_t','g_g'])) )

        df.iloc[:,2:] #/= df[eval(column_name[0])].sum(axis=1).values.reshape(-1,1)
        df = df.replace(np.inf,0).fillna(0).astype(float)


        # normalize the data to later rescale the data to between 0-1
        dfc = df.copy()
        errors = {'a': a, 'c': c, 'g': g, 't': t}
        
        if normalize:
            row_sums = dfc[mutations].sum(axis=1).values.reshape(-1,1) + 10e-6
            dfc[mutations] = dfc[mutations].div(row_sums) 
        
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
                    x=np.arange(1,df1.shape[0]+1),
                    y=df1[mut],
                    mode='markers',
                    name= mut + ' - R1'
                ),
                row=1, col=1
            )

            fig.add_trace(
                go.Scatter(
                    x = np.arange(1,df2.shape[0]+1),
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
    fig_json_ic = table_to_json(table[ 'intersection_C'], normalize=False)
    fig_json_cc = table_to_json(table['complement_C'], normalize=False)

    fig_json_i_norm = table_to_json(table['intersection'])
    fig_json_c_norm = table_to_json(table['complement'])
    fig_json_n_norm = table_to_json(table['non_intersection'])
    fig_json_ic_norm = table_to_json(table[ 'intersection_C'])
    fig_json_cc_norm = table_to_json(table['complement_C'])


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

            .sidebar {{
              height: 100%;
              width: 180px;
              position: fixed;
              z-index: 1;
              top: 0;
              left: 0;
              background-color: #111;
              overflow-x: hidden;
              padding-top: 16px;
            }}

            .sidebar a {{
              padding: 6px 8px 6px 16px;
              text-decoration: none;
              font-size: 20px;
              color: #818181;
              display: block;
            }}

            .sidebar a:hover {{
              color: #f1f1f1;
            }}

            .main {{
              margin-left: 160px; /* Same as the width of the sidenav */
              padding: 0px 10px;
            }}

            @media screen and (max-height: 450px) {{
              .sidebar {{padding-top: 15px;}}
              .sidebar a {{font-size: 18px;}}
            }}
        </style>

    </head>
    <body>
        <script type="text/javascript">
           function toggle_normalize(id1, id2) {{
               var e = document.getElementById(id1);
               var f = document.getElementById(id2);

               if(e.style.display == 'block') {{
                  e.style.display = 'none';
                  f.style.display = 'block';
                }}
               else {{
                  e.style.display = 'block';
                  f.style.display = 'none';
                }}
           }}
        </script>
        
         <div class="sidebar">
          <a href="#" style="font-size:30px; background-color:#8bd9ff;">Raw Counts</a>
          <a href="#section_divPlotly1">Contained</a>
          <a href="#section_divPlotly2">Boundary</a>
          <a href="#section_divPlotly3">non-overlapping</a>
          <a href="#section_divPlotly4">Contained - confidence</a>
          <a href="#section_divPlotly5">Boundary - confidence</a>
            
          <!--<a href="#" style="font-size:30px; background-color:#8bd9ff;">Normalized Counts</a>
          <a href="#divPlotly1_norm">Intersections</a>
          <a href="#divPlotly2_norm">Complementss</a>
          <a href="#divPlotly4_norm">Intersections - confidence</a>
          <a href="#divPlotly5_norm">Complements - confidence</a>-->
        </div>

        <div class="main">

    
            <h1 align="center">Tasmanian artifacts metrics results </h1>
        
            <!-- plot 1 -->

            <div id='section_divPlotly1'>
              <h2 style="padding-left: 40px; padding-top: 90px;">Intersections</h2>
              <h3 style="padding-left: 40px; padding-right: 800; ">Includes all bases that intersect some fragment provided in the bed-file</h3>
              <button class="button button1"; onclick="toggle_normalize('divPlotly1_norm', 'divPlotly1');">Counts/Normalize Counts</button>
              <div id='divPlotly1_norm'>
                <script>
                    var plotly_data = {};
                    Plotly.react('divPlotly1_norm', plotly_data.data, plotly_data.layout);
                </script>
              </div>
              <div id='divPlotly1'>
                <script>
                    var plotly_data2 = {}
                    Plotly.react('divPlotly1', plotly_data2.data, plotly_data2.layout);
                </script>
              </div>
            </div>

            <!-- plot 2 -->

            <div id='section_divPlotly2'>
              <h2 style="padding-left: 40px;padding-top: 60px;">Complements</h2>
              <h3 style="padding-left: 40px; padding-right: 800;">Includes all bases from that do not intersect a fragment, from reads that intersect a fragment provided in the bed-file</h3>
              <button class="button button1"; onclick="toggle_normalize('divPlotly2_norm', 'divPlotly2');">Counts/Normalize Counts</button>
              <div id='divPlotly2_norm'>
                <script>
                    var plotly_data = {};
                    Plotly.react('divPlotly2_norm', plotly_data.data, plotly_data.layout);
                </script>
              </div>
              <div id='divPlotly2'>
                <script>
                    var plotly_data2 = {}
                    Plotly.react('divPlotly2', plotly_data2.data, plotly_data2.layout);
                </script>
              </div>
            </div>

            <!-- plot 3 -->

            <div id='section_divPlotly3'>
              <h2 style="padding-left: 40px;padding-top: 60px;">Non-intersections</h2>
              <h3 style="padding-left: 40px; padding-right: 800;">Includes all bases from reads with no intersections with the bed-file</h3>
              <button class="button button1"; onclick="toggle_normalize('divPlotly3_norm', 'divPlotly3');">Counts/Normalize Counts</button>
              <div id='divPlotly3_norm'>
                <script>
                    var plotly_data = {};
                    Plotly.react('divPlotly3_norm', plotly_data.data, plotly_data.layout);
                </script>
              </div>
              <div id='divPlotly3'>
                <script>
                    var plotly_data2 = {}
                    Plotly.react('divPlotly3', plotly_data2.data, plotly_data2.layout);
                </script>
              </div>
            </div>

            <!-- plot 4 -->

              <div id='section_divPlotly4'>
              <h2 style="padding-left: 40px; padding-top: 90px;">Intersections confidence</h2>
              <h3 style="padding-left: 40px; padding-right: 800; ">Includes all bases that intersect some fragment provided in the bed-file in confidence reads</h3>
              <button class="button button1"; onclick="toggle_normalize('divPlotly4_norm', 'divPlotly4');">Counts/Normalize Counts</button>
              <div id='divPlotly4_norm'>
                <script>
                    var plotly_data = {};
                    Plotly.react('divPlotly4_norm', plotly_data.data, plotly_data.layout);
                </script>
              </div>
              <div id='divPlotly4'>
                <script>
                    var plotly_data2 = {}
                    Plotly.react('divPlotly4', plotly_data2.data, plotly_data2.layout);
                </script>
              </div>
            </div>

            <!-- plot 5 -->

            <div id='section_divPlotly5'>
              <h2 style="padding-left: 40px; padding-top: 90px;">Complement confidence</h2>
              <h3 style="padding-left: 40px; padding-right: 800; ">Includes all complement bases in confidence reads</h3>
              <button class="button button1"; onclick="toggle_normalize('divPlotly5_norm', 'divPlotly5');">Counts/Normalize Counts</button>
              <div id='divPlotly5_norm'>
                <script>
                    var plotly_data = {};
                    Plotly.react('divPlotly5_norm', plotly_data.data, plotly_data.layout);
                </script>
              </div>
              <div id='divPlotly5'>
                <script>
                    var plotly_data2 = {}
                    Plotly.react('divPlotly5', plotly_data2.data, plotly_data2.layout);
                </script>
              </div>
            </div>

            <script> 
                divPlotly1_norm.style.display = 'none'; 
                divPlotly2_norm.style.display = 'none'; 
                divPlotly3_norm.style.display = 'none';
                divPlotly4_norm.style.display = 'none';
                divPlotly5_norm.style.display = 'none';
            </script>
        </div> <!-- finish with class main here -->
    </body>

    </html>"""

    # write the JSON to the HTML template
    #with open('Tasmanian_artifact_metrics_report.html', 'w') as f:
    #   f.write(template.format(fig_json))
    return template.format(fig_json_i_norm, fig_json_i, fig_json_c_norm, fig_json_c, fig_json_n_norm, fig_json_n, fig_json_ic_norm, fig_json_ic, fig_json_cc_norm, fig_json_cc)



#if __name__=='__main__':
#    
#    normalize = False
#    for n,i in enumerate(sys.argv):
#        if i in ["-normalize","--normalize","-n","--n","--norm","-norm"]:
#            normalize=True
#        if i in ["--table", "-t","--t","-table"]:
#            

