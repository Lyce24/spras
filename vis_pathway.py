from pathlib import Path

import dash_cytoscape as cyto
import networkx as nx
from dash import Dash, Input, Output, callback, dcc, html

app = Dash(__name__)


algo_list = []

source_gene_file = Path('./input/source_gene.txt')
# Check whether the directory FlyBase exists

cell_target_gene_file = Path('./input/cell_target_gene.txt')
muscle_target_gene_file = Path('./input/muscle_target_gene.txt')

source_gene = set()
cell_target_gene = set()
muscle_target_gene = set()
with open(source_gene_file, 'r') as f:
    for i, line in enumerate(f):
        # if the first line is the title, skip it
        if i == 0:
            continue
        line = line.strip()
        endpoints = line.split("\t")
        source_gene.add(endpoints[0])

with open(cell_target_gene_file, 'r') as f:
    for i, line in enumerate(f):
        # if the first line is the title, skip it
        if i == 0:
            continue
        line = line.strip()
        endpoints = line.split("\t")
        cell_target_gene.add(endpoints[0])

with open(muscle_target_gene_file, 'r') as f:
    for i, line in enumerate(f):
        # if the first line is the title, skip it
        if i == 0:
            continue
        line = line.strip()
        endpoints = line.split("\t")
        muscle_target_gene.add(endpoints[0])



parent_dir = Path("./output")

cell_parent_dir = parent_dir / 'cell_cell_fusion'

pathlinker_runs = [x for x in cell_parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('cell-cell-fusion-pathlinker-params-')]
rwr_runs = [x for x in cell_parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('cell-cell-fusion-rwr-params-')]
# tiedie_runs = [x for x in parent_dir.iterdir() if x.is_dir() and x.name.startswith('flybase-tiedie-params')]

pathlinker_dict = {}

for i, run in enumerate(pathlinker_runs):
    print(run)
    temp = str(run).split('-')
    pathlinker_dict[temp[-1]] = i
    algo_list.append(f'pathlinker {temp[-1]}')
    

print(pathlinker_dict)    


def produce_elements(algo, mode):
    if algo == 'pathlinker':
        run = pathlinker_dict[mode]
        file = pathlinker_runs[run] / 'pathway.txt'
        edges = []
        with open(file, 'r') as f:
            for line in f:
                temp = line.strip().split('\t')
                edges.append((temp[0], temp[1]))
        G = nx.Graph()
        G.add_edges_from(edges)
        print(f"Using RWR algorithm with {mode}...")
        print(f"Number of nodes: {len(G.nodes)}")
        print(f"Number of edges: {len(G.edges)} \n")

    source = 0
    intermediate = 0
    target = 0
    elements = nx.cytoscape_data(G)
    mid = len(elements['elements']['nodes']) // 2
    for i in range(len(elements['elements']['nodes'])):
        elements['elements']['nodes'][i]['data']['label'] = elements['elements']['nodes'][i]['data']['id']
        if elements['elements']['nodes'][i]['data']['id'] in source_gene:
            elements['elements']['nodes'][i]['classes'] = 'red triangle'
            # if i < mid:
            #     x_coordinate = 40 + 50 * (mid - i)
            # else:
            #     x_coordinate = 40 + 50 * (i - mid)
            elements['elements']['nodes'][i]['position'] = {
                'x': 40 * source, 'y': 0}
            source += 1
        elif elements['elements']['nodes'][i]['data']['id'] in cell_target_gene:
            elements['elements']['nodes'][i]['classes'] = 'blue rectangle'
            # if i < mid:
            #     x_coordinate = 40 + 50 * (mid - i)
            # else:
            #     x_coordinate = 40 + 50 * (i - mid)
            elements['elements']['nodes'][i]['position'] = {
                'x': 40 * target, 'y': 1000}    
            target += 1 
        else:
            # if i < mid:
            #     x_coordinate = 40 + 50 * (mid - i)
            # else:
            #     x_coordinate = 40 + 50 * (i - mid)
            elements['elements']['nodes'][i]['position'] = {
                'x': 40 * intermediate, 'y': 500}
            intermediate += 1

    for i in range(len(elements['elements']['edges'])):
        if elements['elements']['edges'][i]['data']['source'] in source_gene and elements['elements']['edges'][i]['data']['target'] in source_gene:
            elements['elements']['edges'][i]['classes'] = 'red'
            
    for i in range(len(elements['elements']['edges'])):
        if elements['elements']['edges'][i]['data']['source'] in cell_target_gene and elements['elements']['edges'][i]['data']['target'] in cell_target_gene:
            elements['elements']['edges'][i]['classes'] = 'blue'

    return elements['elements']


app.layout = html.Div([
    dcc.Dropdown(
        id='dropdown-update-layout',
        value='preset',
        clearable=False,
        options=[
            {'label': name.capitalize(), 'value': name}
            for name in ['preset', 'grid', 'random', 'circle', 'cose', 'concentric']
        ]
    ),

    dcc.Dropdown(
        id='dropdown-update-algo',
        value=f'pathlinker {str(pathlinker_runs[0]).split("-")[-1]}',
        clearable=False,
        options=[
            {'label': name, 'value': name}
            for name in algo_list
        ]
    ),

    cyto.Cytoscape(
        id='flybase_rwr',
        layout={'name': 'preset'},
        style={'width': '100%', 'height': '960px'},
        elements=produce_elements('pathlinker', str(pathlinker_runs[0]).split('-')[-1]),
        stylesheet=[
            # Group selectors
            {
                'selector': 'node',
                'style': {
                    'content': 'data(label)'
                }
            },

            # Class selectors
            {
                'selector': '.red',
                'style': {
                    'background-color': 'red',
                    'line-color': 'red'
                }
            },
            {
                'selector': '.triangle',
                'style': {
                    'shape': 'triangle'
                }
            },
            {
                'selector': '.rectangle',
                'style': {
                    'shape': 'rectangle'
                }
            },
            {
                'selector': '.blue',
                'style': {
                    'background-color': 'blue',
                    'line-color': 'blue'
                }
            }
        ]
    ),
    

    dcc.Markdown(id='cytoscape-selectedNodeData-markdown')
])


@callback(Output('cytoscape-selectedNodeData-markdown', 'children'),
          Input('flybase_rwr', 'selectedNodeData'))
def displaySelectedNodeData(data_list):
    if data_list is None:
        return "No node selected."

    node_list = [data['label'] for data in data_list]
    return "You selected the following nodes: " + "\t ".join(node_list)


@callback(Output('flybase_rwr', 'layout'),
          Input('dropdown-update-layout', 'value'))
def update_layout(layout):
    return {
        'name': layout,
        'animate': True
    }


@callback(Output('flybase_rwr', 'elements'),
          Input('dropdown-update-algo', 'value'))
def update_algo(algo):
    algo_list = algo.strip().split(' ')
    return produce_elements(algo_list[0], algo_list[1])


if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
