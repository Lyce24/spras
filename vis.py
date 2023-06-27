import dash_cytoscape as cyto
import networkx as nx
from dash import Dash, html, dcc, Input, Output, callback
from pathlib import Path

app = Dash(__name__)

source_gene_file = Path('./input/source_gene.txt')
rwr_pathway_file = Path('./FlyBase/rwr-pathway.txt')


source_gene = set()
with open(source_gene_file, 'r') as f:
    for i, line in enumerate(f):
        # if the first line is the title, skip it
        if i == 0:
            continue
        line = line.strip()
        endpoints = line.split("\t")
        if len(endpoints) != 2:
            raise ValueError(f"Edge {line} does not contain 2 nodes separated by ' ' and a weight")
        source_gene.add(endpoints[0])
        
        
algo_list = ['rwr mode1', 'rwr mode2', 'rwr mode3']


# find all runs of oi1 in the parent directory of rwr_pathway_file
parent_dir = rwr_pathway_file.parent

oi1_runs = [x for x in parent_dir.iterdir() if x.is_dir() and x.name.startswith('flybase-omicsintegrator1-params')]
oi2_runs = [x for x in parent_dir.iterdir() if x.is_dir() and x.name.startswith('flybase-omicsintegrator2-params')]

oi1_dict = {}
oi2_dict = {}

for i, run in enumerate(oi1_runs):
    temp = str(run).split('-')
    oi1_dict[temp[-1]] = i
    algo_list.append(f'oi1 {temp[-1]}')


for i, run in enumerate(oi2_runs):
    temp = str(run).split('-')
    oi2_dict[temp[-1]] = i
    algo_list.append(f'oi2 {temp[-1]}')

def produce_elements(algo, mode):
    if algo == 'rwr':
        file = rwr_pathway_file

        edges = []
        linker_nodes = set()
        with open(file, 'r') as f:
            for line in f:
                temp = line.strip().split('\t')
                edges.append((temp[0], temp[1]))
                linker_nodes.add(temp[0])
                linker_nodes.add(temp[1])
                
        selected_edges = []
        if mode == 'mode1':
            selected_edges = edges
        elif mode == 'mode2':
            for edge in edges:
                if edge[0] in source_gene and edge[1] in linker_nodes:
                    selected_edges.append(edge)
        elif mode == 'mode3':
            for edge in edges:
                if edge[0] in source_gene and edge[1] in (linker_nodes - source_gene):
                    selected_edges.append(edge)
        G = nx.Graph()
        G.add_nodes_from(linker_nodes)
        G.add_edges_from(selected_edges)
        print(f"Using RWR algorithm with {mode}...")
        print(f"Number of nodes: {len(G.nodes)}")
        print(f"Number of edges: {len(G.edges)} \n")
    elif algo == 'oi1':
        run = oi1_dict[mode]
        file = oi1_runs[run] / 'pathway.txt'
        edges = []
        with open(file, 'r') as f:
            for line in f:
                temp = line.strip().split('\t')
                edges.append((temp[0], temp[1]))
        G = nx.Graph()
        G.add_edges_from(edges)
        print(f"Using OI1 algorithm with run {run+1}...")
        print(f"Number of nodes: {len(G.nodes)}")
        print(f"Number of edges: {len(G.edges)} \n")
    
    elif algo == 'oi2':
        run = oi2_dict[mode]
        file = oi2_runs[run] / 'pathway.txt'
        edges = []
        with open(file, 'r') as f:
            for line in f:
                temp = line.strip().split('\t')
                edges.append((temp[0], temp[1]))
        G = nx.Graph()
        G.add_edges_from(edges)
        print(f"Using OI2 algorithm with run {run+1}...")
        print(f"Number of nodes: {len(G.nodes)}")
        print(f"Number of edges: {len(G.edges)} \n")

    elements = nx.cytoscape_data(G)
    mid = len(elements['elements']['nodes']) // 2
    for i in range(len(elements['elements']['nodes'])):
        elements['elements']['nodes'][i]['data']['label'] = elements['elements']['nodes'][i]['data']['id']
        if elements['elements']['nodes'][i]['data']['id'] in source_gene:
            elements['elements']['nodes'][i]['classes'] = 'red triangle'
            if i < mid:
                x_coordinate = 40 + 50 * (mid - i)
            else:
                x_coordinate = 40 + 50 * (i - mid)
            elements['elements']['nodes'][i]['position'] = {'x': x_coordinate, 'y': 40 * i}
        else:
            if i < mid:
                x_coordinate = 120 * mid  + 50 * (i - mid)
            else:
                x_coordinate = 120 * mid + 50 * (mid - i)
            elements['elements']['nodes'][i]['position'] = {'x': x_coordinate, 'y': 40 * i}

    for i in range(len(elements['elements']['edges'])):
        if elements['elements']['edges'][i]['data']['source'] in source_gene and elements['elements']['edges'][i]['data']['target'] in source_gene:
            elements['elements']['edges'][i]['classes'] = 'red'

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
        id = 'dropdown-update-algo',
        value = 'rwr mode1',
        clearable = False,
        options = [
            {'label': name, 'value': name}
            for name in algo_list
        ]
    ),
    
    cyto.Cytoscape(
        id='flybase_rwr',
        layout={'name': 'preset'},
        style={'width': '100%', 'height': '960px'},
        elements=produce_elements('rwr', 'mode1'),
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
    app.run_server(debug=True, port = 8050)
