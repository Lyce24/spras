from pathlib import Path

import dash_cytoscape as cyto
import networkx as nx
from dash import Dash, Input, Output, callback, dcc, html

app = Dash(__name__)


algo_list = []

source_gene_file = Path('./input/source_gene.txt')
# Check whether the directory FlyBase exists

if Path('./FlyBase/rwr-pathway.txt').exists():
    rwr_pathway_file = Path('./FlyBase/rwr-pathway.txt')
    algo_list.append('rwr mode1')
    algo_list.append('rwr mode2')
    algo_list.append('rwr mode3')

source_gene = set()
with open(source_gene_file, 'r') as f:
    for i, line in enumerate(f):
        # if the first line is the title, skip it
        if i == 0:
            continue
        line = line.strip()
        endpoints = line.split("\t")
        if len(endpoints) != 2:
            raise ValueError(
                f"Edge {line} does not contain 2 nodes separated by ' ' and a weight")
        source_gene.add(endpoints[0])


# find all runs of oi1 in the parent directory of rwr_pathway_file
parent_dir = Path("./FlyBase")

oi1_runs = [x for x in parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('flybase-omicsintegrator1-params')]
oi2_runs = [x for x in parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('flybase-omicsintegrator2-params')]

oi1_dict = {}
oi2_dict = {}

possible_params_oi1 = {}
possible_params_oi1['b'] = {}
possible_params_oi1['w'] = {}
possible_params_oi1['mu'] = {}

valid_runs = []


with open(parent_dir / 'flybase-pathway-summary.txt', 'r') as f:
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.strip()
        temp = line.split("\t")
        if int(temp[-1]) == 34:
            name = temp[0].split('/')[1].split('-')[-1]
            valid_runs.append(name)
            if int(temp[1]) >= 40:
                print(name)
                with open(parent_dir / 'logs' / f'parameters-omicsintegrator1-params-{name}.yaml', 'r') as f1:
                    params = ''
                    for j, line1 in enumerate(f1):
                        if j == 0:
                            params += f'b: {line1.strip().split(" ")[-1]} '
                        elif j == 3:
                            params += f'mu: {line1.strip().split(" ")[-1]} '
                        elif j == 5:
                            params += f'w: {line1.strip().split(" ")[-1]} '
                        else:
                            continue
                    print(params)

for i, run in enumerate(oi1_runs):
    # check if the lines of the file is between 50 - 100
    temp = str(run).split('-')
    if temp[-1] in valid_runs:
        with open(parent_dir / 'logs' / f'parameters-omicsintegrator1-params-{temp[-1]}.yaml', 'r') as f:
            for j, line in enumerate(f):
                if j == 0:
                    line = line.strip()
                    params = line.split(" ")
                    if params[-1] in possible_params_oi1['b']:
                        possible_params_oi1['b'][params[-1]] += 1
                    else:
                        possible_params_oi1['b'][params[-1]] = 1

                elif j == 3:
                    line = line.strip()
                    params = line.split(" ")
                    if params[-1] in possible_params_oi1['mu']:
                        possible_params_oi1['mu'][params[-1]] += 1
                    else:
                        possible_params_oi1['mu'][params[-1]] = 1
                elif j == 5:
                    line = line.strip()
                    params = line.split(" ")
                    if params[-1] in possible_params_oi1['w']:
                        possible_params_oi1['w'][params[-1]] += 1
                    else:
                        possible_params_oi1['w'][params[-1]] = 1
                else:
                    continue
        oi1_dict[temp[-1]] = i
        algo_list.append(f'oi1 {temp[-1]}')


print(possible_params_oi1)

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
        print(f"Using OI1 algorithm with run {mode}...")
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
        print(f"Using OI2 algorithm with run {mode}...")
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
            elements['elements']['nodes'][i]['position'] = {
                'x': x_coordinate, 'y': 40 * i}
        else:
            if i < mid:
                x_coordinate = 120 * mid + 50 * (i - mid)
            else:
                x_coordinate = 120 * mid + 50 * (mid - i)
            elements['elements']['nodes'][i]['position'] = {
                'x': x_coordinate, 'y': 40 * i}

    for i in range(len(elements['elements']['edges'])):
        if elements['elements']['edges'][i]['data']['source'] in source_gene and elements['elements']['edges'][i]['data']['target'] in source_gene:
            elements['elements']['edges'][i]['classes'] = 'red'

    return elements['elements']


app.layout = html.Div([
    dcc.Dropdown(
        id='dropdown-update-layout',
        value='cose',
        clearable=False,
        options=[
            {'label': name.capitalize(), 'value': name}
            for name in ['preset', 'grid', 'random', 'circle', 'cose', 'concentric']
        ]
    ),

    dcc.Dropdown(
        id='dropdown-update-algo',
        value=f'oi1 {str(oi1_runs[0]).split("-")[-1]}',
        clearable=False,
        options=[
            {'label': name, 'value': name}
            for name in algo_list
        ]
    ),

    cyto.Cytoscape(
        id='flybase_rwr',
        layout={'name': 'cose'},
        style={'width': '100%', 'height': '960px'},
        elements=produce_elements('oi1', str(oi1_runs[0]).split('-')[-1]),
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
    app.run_server(debug=True, port=8050)
