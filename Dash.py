import base64
import datetime
import io
import uuid
import os
from hashlib import md5
import binascii

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go

import pandas as pd

import sanger_analysis_backend as sab

################################################################################

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.scripts.config.serve_locally = True

app.title = 'VVSS'

################################################################################

class appState(object):
    def __init__(self):
        self.primer_table  = None
        self.well_table    = None

        self.abi_dict      = { }
        self.chromatogram_figures = { }
        self.chromatogram_blasts  = None
        self.bres_abi             = None

        self.view          = None
        self.view_data     = { }

        self.buttonCounts = { }
    #edef

    def getButtonCount(self, id):
        if id not in self.buttonCounts:
            self.buttonCounts[id] = 0
        #fi
        return self.buttonCounts[id]
    #edef

    def incrementButtonCount(self, id):
        self.buttonCounts[id] = self.getButtonCount(id) + 1
    #edef

#eclass

aState = appState()

################################################################################
# Define the DASH layout

# +--------------------------------------+
# |[title_bar]                           |
# |+--------------------------++--------++
# ||[view_select_bar]         ||[upload]||
# ||                          ||        ||
# ||                          ||        ||
# |+--------------------------++--------+|
# +--------------------------------------+


upload_button_style = {
    'width': '100%',
    'height': '60px',
    'lineHeight': '60px',
    'borderWidth': '1px',
    'borderStyle': 'dashed',
    'borderRadius': '5px',
    'textAlign': 'center',
    'margin': '10px'
}

upload_side = html.Div([
    dcc.Upload(
        id='upload-primer-table',
        children=html.Div(['Drag and Drop or ', html.A('Select Primer Table File')]),
        style=upload_button_style,
        multiple=False
    ),
    html.Div(id='hidden-primer-upload-div', style={'display':'none'}),
    dcc.Upload(
        id='upload-well-table',
        children=html.Div(['Drag and Drop or ', html.A('Select Well Table File')]),
        style=upload_button_style,
        multiple=False
    ),
    dcc.Upload(
        id='upload-abi-files',
        children=html.Div(['Drag and Drop or ', html.A('Select AB1/ABI Files')]),
        style=upload_button_style,
        multiple=True
    ),
    html.Button('Align ABI to variants', id='button-process', style={'width': '100%', 'margin': '10px'}),
    html.Div(id='hidden-process-div', style={'display':'none'}),

])

title_bar = html.Div(children=[html.H1('Variant Validation with Sanger Sequencing')])


view_frame = html.Div(children=[
  dcc.Tabs(id="view-tabs", children=[
        dcc.Tab(label='Primer Table', children=[
                html.Button('Load Primer Table', id='button-load-primers', style={'width': '45%', 'margin': '10px'}),
                html.Button('Check for mismatches', id='button-process-primers', style={'width': '45%', 'margin': '10px'}),
                html.Div(id='primer-table-view')]),
        dcc.Tab(label='Well Table', children=[html.Div(id='well-table-view')]),
        dcc.Tab(label='Chromatograms',
                children=[html.Div(children=[
                            dcc.Dropdown(
                                id='chromatogram-selection',
                                options=[{"label":"No AB1 files uploaded", "value":None}],
                                value=None),
                            ]),
                          html.Div(id='chromatogram-plots')])
        ])],
      id = 'view_frame',
      style={'width': '100%',
          'lineHeight': '60px',
          'borderWidth': '1px',
          'borderStyle': 'dashed',
          'borderRadius': '5px',
          'textAlign': 'center',
          'margin': '10px'})

app.layout = html.Div(children=[
       title_bar,
       html.Div(children=[
          html.Div(children=[view_frame],
            style={'width':'75%',
                   'margin':'1%'}),
          html.Div(children=[upload_side],
            style={'width':'20%',
                   'margin':'1%'})
          ],
          style={'display': 'flex',
                 'align':'center'})
                 ])

################################################################################

def generate_table(df, names=None):
    header = [html.Tr([html.Th(col) for col in (df.columns if names is None else names)])]
    body = [html.Tr([html.Td(str(r[col])) for col in df.columns]) for i,r in df.iterrows() ]
    return html.Table(header + body, style={'margin':'10px'})
#edef

def parse_tsv_contents(contents, filename):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        return pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t')
    except Exception as e:
        print("ERROR")
        print(e)
        return None
    #etry

    return None
#edef

def parse_abi_contents(contents, filenames):
    D = {}
    sab.biu.utils.fs.mkdirp('uploaded_abi')

    for c, f in zip(contents, filenames):
        content_type, content_string = c.split(',')
        decoded = base64.b64decode(content_string)
        m = md5()
        m.update(decoded)
        chksum = m.hexdigest()
        fname = 'uploaded_abi/%s.abi' % chksum
        with open(fname, 'wb') as ofd:
            ofd.write(binascii.a2b_base64(content_string))
        #ewith
        D[f] = sab.biu.formats.ABI(fname)
    #efor
    return D
#edef

################################################################################

def draw_primer_view():
    if aState.primer_table is None:
        return html.Div(children=[ "Please upload a primer table first", html.Br(),
                                   "The table you upload MUST contain columns: 'variantid', 'primer', 'direction', 'sequence'.", html.Br(),
                                   "  variantid: Variant identifier, String with format: 'chr-pos-ref-alt", html.Br(),
                                   "  primer: Primar name, String", html.Br(),
                                   "  direction: Primer direction, String: [Forward|Reverse]", html.Br(),
                                   "  sequence:  Primer sequence, String: [AaGgCcTt]+", html.Br() ])
    #fi
    print('DRAW ATTEMPT', aState.primer_table.columns)
    df = aState.primer_table

    return html.Div(children=[generate_table(aState.primer_table)])
#edef

def draw_well_view():
    if aState.well_table is None:
        return html.Div(children=[ "Please upload a file table first", html.Br(),
                           "The table you upload MUST contain columns: 'filename', 'variantid'.", html.Br(),
                           "  filename: The name of the ABI file", html.Br(),
                           "  variantid: Variant identifier, String with format: 'chr-pos-ref-alt", html.Br() ])
    #fi

    return html.Div(children=generate_table(aState.well_table))
#edef

################################################################################
# Callbacks

@app.callback(Output('primer-table-view', 'children'),
              [Input('button-process-primers', 'n_clicks'),
               Input('button-load-primers', 'n_clicks')])
def load_or_process_primers(process_clicks, load_clicks):

    if (process_clicks is not None) and (aState.getButtonCount('process_primer') < process_clicks):
        aState.incrementButtonCount('process_primer')
        if (aState.primer_table is not None) and 'true_vid' not in aState.primer_table.columns:
            bres, incorrect = sab.determine_primer_alignments(aState.primer_table)
            aState.view_data['alignments'] = bres
            aState.view_data['incorrect_alignments'] = incorrect
            aState.primer_table = aState.primer_table.join(bres['true_vid'], on='primer')
        #fi
    elif (load_clicks is not None) and (aState.getButtonCount('load_primer') < load_clicks):
        aState.incrementButtonCount('load_primer')
    else:
        return html.Div(children=[ "Please upload a primer table first", html.Br(),
                                   "The table you upload MUST contain columns: 'variantid', 'primer', 'direction', 'sequence'.", html.Br(),
                                   "  variantid: Variant identifier, String with format: 'chr-pos-ref-alt", html.Br(),
                                   "  primer: Primar name, String", html.Br(),
                                   "  direction: Primer direction, String: [Forward|Reverse]", html.Br(),
                                   "  sequence:  Primer sequence, String: [AaGgCcTt]+", html.Br() ])
    #fi
    return draw_primer_view()
#edef


@app.callback(Output('hidden-process-div', 'children'),
              [Input('button-process', 'n_clicks')])
def process_inputs(process_clicks):
    if process_clicks is None:
        return None
    #fi
    if aState.getButtonCount('process') >= process_clicks:
        return None
    #fi

    if (aState.chromatogram_blasts is None) and (len(aState.abi_dict) > 0):
        bres = sab.align_abi_files(aState.abi_dict)
        if aState.chromatogram_blasts is not None:
            bres = pd.concat([aState.chromatogram_blasts, bres])
            bres = bres.drop_duplicates()
        #fi
        bres = bres[bres.groupby('qseqid').transform(max).bitscore == bres.bitscore]
        aState.chromatogram_blasts = bres
    #fi

    if ((aState.chromatogram_blasts is not None) and
       (aState.primer_table is not None)):
        print(aState.chromatogram_blasts)
        aState.bres_abi = sab.process_abi_alignments(aState.primer_table,
                                                     aState.chromatogram_blasts.copy(),
                                                     aState.well_table)

    #fi
    return None
#edef

@app.callback(Output('chromatogram-plots', 'children'),
              [Input('chromatogram-selection', 'value')])
def draw_chromatogram(selected_abi):
    if selected_abi is None:
        return None
    #fi

    figures = []
    if selected_abi not in aState.chromatogram_figures:
        aState.chromatogram_figures[selected_abi] = {}
    #fi

    abi = aState.abi_dict[selected_abi]

    if aState.bres_abi is not None:
        print(aState.bres_abi.columns)
        if 'var_location' not in aState.chromatogram_figures[selected_abi]:
            print('making var_location figures')
            aState.chromatogram_figures[selected_abi]['var_location'] = {}
            bres_abi = aState.bres_abi
            for i, row in bres_abi[bres_abi.qseqid == selected_abi].iterrows():
                chromat_pos = row.chromat_pos
                print(chromat_pos)
                print(row.sseqid)
                print(list(aState.bres_abi.columns))
                keys = list(row.chromat_pos.keys())
                if ('variantid' in row.index) and (row.variantid in keys):
                    keys = [row.variantid]
                #fi
                for true_vid in keys:
                    pos = chromat_pos[true_vid]
                    ax = abi.chromatogram(xlim=(pos-10, pos+10), highlight=(pos-.5,pos+.5))
                    primer = row.primer if 'primer' in row.index else "UNKNOWN"
                    variantid = row.variantid if 'variantid' in row.index else "UNKNOWN"
                    VID = row.primer if 'primer' in row.index else 'UNKNOWN'
                    true_vid_primer = aState.primer_table[aState.primer_table.variantid == true_vid].primer.values[0]
                    ax.set_title('%s | Goal variant: %s\nAligned to: %s -> %s' % (selected_abi, row.variantid, true_vid_primer, true_vid), fontsize=5)
                    aState.chromatogram_figures[selected_abi]['var_location'][true_vid] = (sab.fig_to_uri(ax.figure))
                #efor
            #efor
        #fi
    #fi
    if 'full' not in aState.chromatogram_figures[selected_abi]:
        aState.chromatogram_figures[selected_abi]['full'] = sab.fig_to_uri(abi.chromatogram().figure)
    #fi

    figures = {}
    if 'var_location' in aState.chromatogram_figures[selected_abi]:
        figures = aState.chromatogram_figures[selected_abi]['var_location']
    #fi


    images = [ html.Img(src=aState.chromatogram_figures[selected_abi]['full'], style={"width":"100%"}) ]

    for var_loc, source in figures.items():
        
        images.append(html.H2(children=[var_loc]))
        if len(var_loc.split('-')) == 4:
            chrom, pos, ref, alt = var_loc.split('-')
            pos = int(pos)-1
            context = 3
            fwd_seq = sab.hg.genome[chrom][pos-context:pos+context+1]
            rev_seq = fwd_seq.revcomp().seq
            fwd_seq = fwd_seq.seq
            context_seqs = [
                'Forward: ' + fwd_seq[:context] + '[%s/%s]' % (fwd_seq[context], alt) + fwd_seq[context+1:],
                html.Br(),
                html.Span(children='s'*9, style={"color":"white"}), '-'*context + '[R/A]' + '-'*context,
                html.Br(),
                'Reverse: ' + rev_seq[:context] + '[%s/%s]' % (rev_seq[context], sab.biu.formats.Sequence(None, alt).revcomp().seq) \
                            + rev_seq[context+1:]
                ]
            print(context_seqs)
            images.append(html.Div(children=context_seqs,
                style={"line-height" : "1em",
                       "text-align": "left",
                       "font-size" : "small",
                       "width" : "50%",
                       "font-family":"Consolas,Monaco,Lucida Console,Liberation Mono,DejaVu Sans Mono,Bitstream Vera Sans Mono,Courier New, monospace"}))
        #fi
        images.append(html.Img(src=source, style={"width":"100%"})),
        images.append(html.Br())
        images.append(html.Br())
    #efor
    images = html.Div(children=images, style={"width":"100%"})
    info   = [ "Sequence length : %d" % len(aState.abi_dict[selected_abi].seq) ]
    if aState.well_table is not None:
        print("'%s'" % selected_abi)

        rel_rows = aState.well_table[(aState.well_table.filename == selected_abi)]
        info.append(generate_table(rel_rows))
    #fi
    info = html.Div(children=info, style={"width":"100%"})
    return html.Div(children=[info, images])
#edef

@app.callback(Output('hidden-primer-upload-div', 'children'),
              [Input('upload-primer-table', 'contents')],
              [State('upload-primer-table', 'filename'),
               State('upload-primer-table', 'last_modified')])
def upload_primer_table(contents, names, dates):
    if contents is not None:
        data = parse_tsv_contents(contents, names)

        expected_columns = set([ 'variantid', 'primer', 'direction', 'sequence' ])

        # Only accept a primer table with the above columns
        if len(set(data.columns) & expected_columns) != len(expected_columns):
            aState.primer_table = None
            return None
        #fi

        if ((aState.primer_table is None) or
            not(aState.primer_table.equals(data[[ c for c in aState.primer_table.columns if c in data.columns ]]))):
            aState.primer_table = data
            aState.bres_abi = None
        #fi

    #fi

    return None #draw_primer_view()
#edef

@app.callback(Output('well-table-view', 'children'),
              [Input('upload-well-table', 'contents')],
              [State('upload-well-table', 'filename'),
               State('upload-well-table', 'last_modified')])
def upload_well_table(contents, names, dates):
    if contents is not None:
        data = parse_tsv_contents(contents, names)
        print(data)
        expected_columns = set(['filename', 'variantid'])

        # Only accept a primer table with the above columns
        if len(set(data.columns) & expected_columns) != len(expected_columns):
            aState.well_table = None
            data = None
        elif ((aState.well_table is None) or not(aState.well_table.equals(data))):
            aState.well_table = data
            aState.bres_abi   = None
        #fi
    #fi

    return draw_well_view()
#edef

@app.callback(Output('chromatogram-selection', 'options'),
              [Input('upload-abi-files', 'contents')],
              [State('upload-abi-files', 'filename'),
               State('upload-abi-files', 'last_modified')])
def upload_abi_files(contents, names, dates):

    if contents is not None:
        D = parse_abi_contents(contents, names)

        for k, v in D.items():
            aState.abi_dict[k] = v
            aState.chromatogram_blasts = None
        #efor
    #fi

    if len(aState.abi_dict) == 0:
        return [{"label":"No AB1 files uploaded", "value":None}]
    else:
        return sorted([ {"label":k, "value":k} for (k,v) in aState.abi_dict.items() ], key=lambda x: x['label'])
    #fi
#edef

################################################################################

print(sab.biu.config.settings.platform())

if __name__ == '__main__':
    app.run_server(debug=True)
