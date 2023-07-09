from pathlib import Path

import dash_cytoscape as cyto
import networkx as nx
from dash import Dash, Input, Output, callback, dcc, html

app = Dash(__name__)

print("Loading data...")
source_gene_file = Path('./input/source_gene.txt')
receptor_gene_file = Path('./input/receptor.txt')
parent_dir = Path("./output/amigo2")
flybase_parent_dir = Path("./output/flybase")

cell_target_gene_file = Path('./input/cell_target_gene.txt')
muscle_target_gene_file = Path('./input/muscle_target_gene.txt')

source_gene = set()
receptor_gene = set()
cell_target_gene = set()
muscle_target_gene = set()

oi1_dict = {}
oi2_dict = {}
cell_pathlinker_dict = {}
cell_rwr_dict = {}
muscle_pathlinker_dict = {}
muscle_rwr_dict = {}

oi1_runs = [x for x in flybase_parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('flybase-omicsintegrator1-params')]
oi2_runs = [x for x in flybase_parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('flybase-omicsintegrator2-params')]

cell_pathlinker_runs = [x for x in parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('cell-cell-fusion-pathlinker-params-')]
cell_rwr_runs = [x for x in parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('cell-cell-fusion-rwr-params-')]

muscle_pathlinker_runs = [x for x in parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('muscle-development-pathlinker-params-')]
muscle_rwr_runs = [x for x in parent_dir.iterdir() if x.is_dir(
) and x.name.startswith('muscle-development-rwr-params-')]

algo_list_flybase = []
algo_list_cell = []
algo_list_muscle = []

if Path('./FlyBase/rwr-pathway.txt').exists():
    rwr_pathway_file = Path('./FlyBase/rwr-pathway.txt')
    algo_list_flybase.append('rwr mode1')
    algo_list_flybase.append('rwr mode2')
    algo_list_flybase.append('rwr mode3')

def generate_nodes(nodes_file: Path, gene_set) -> None:
    """
    This function is for generating the nodes from the path to the source/target file
    """
    with open(nodes_file, 'r') as f:
        for i, line in enumerate(f):
            # if the first line is the title, skip it
            if i == 0:
                continue
            line = line.strip()
            endpoints = line.split("\t")
            gene_set.add(endpoints[0])


def update_algo_list(arr, dict, runs, name):
    for i, run in enumerate(runs):
        temp = str(run).split('-')
        dict[temp[-1]] = i
        arr.append(f'{name} {temp[-1]}')


def grab_data(parent_dir, gene_dict, algo_runs, algo, run_id):
    if run_id not in gene_dict:
        return nx.Graph()
    run = gene_dict[run_id]
    file = algo_runs[run] / 'pathway.txt'
    edges = []
    with open(file, 'r') as f:
        for line in f:
            temp = line.strip().split('\t')
            edges.append((temp[0], temp[1]))
    G = nx.Graph()
    G.add_edges_from(edges)
    params = ''
    with open(parent_dir / 'logs' / f'parameters-{algo}-params-{run_id}.yaml', 'r') as f1:
        for line1 in f1:
            temp1 = line1.strip().split(':')
            params += f'{temp1[0]} - {temp1[1]} \n'
    print(f"Using {algo} algorithm with run {run_id}...")
    print(f"Number of nodes: {len(G.nodes)}")
    print(f"Number of edges: {len(G.edges)}")
    print(f"Parameters : \n{params}")
    return G


def produce_elements(G: nx.Graph, target_gene):
    source = 0
    receptor = 0
    intermediate = 0
    target = 0
    elements = nx.cytoscape_data(G)
    for i in range(len(elements['elements']['nodes'])):
        elements['elements']['nodes'][i]['data']['label'] = elements['elements']['nodes'][i]['data']['id']
        if elements['elements']['nodes'][i]['data']['id'] in source_gene:
            elements['elements']['nodes'][i]['classes'] = 'red triangle'
            elements['elements']['nodes'][i]['position'] = {
                'x': 60 * source, 'y': 0}
            source += 1
        elif target_gene is not None and elements['elements']['nodes'][i]['data']['id'] in target_gene:
            elements['elements']['nodes'][i]['classes'] = 'blue rectangle'
            elements['elements']['nodes'][i]['position'] = {
                'x': 60 * target, 'y': 1600}
            target += 1
        elif elements['elements']['nodes'][i]['data']['id'] in receptor_gene:
            elements['elements']['nodes'][i]['classes'] = 'green rectangle'
            elements['elements']['nodes'][i]['position'] = {
                'x': 60 * receptor, 'y': 500}
            receptor += 1
        else:
            elements['elements']['nodes'][i]['position'] = {
                'x': 60 * intermediate, 'y': 1000}
            intermediate += 1

    for i in range(len(elements['elements']['edges'])):
        if elements['elements']['edges'][i]['data']['source'] in source_gene and elements['elements']['edges'][i]['data']['target'] in source_gene:
            elements['elements']['edges'][i]['classes'] = 'red'
    if target_gene is not None:
        for i in range(len(elements['elements']['edges'])):
            if elements['elements']['edges'][i]['data']['source'] in target_gene and elements['elements']['edges'][i]['data']['target'] in target_gene:
                elements['elements']['edges'][i]['classes'] = 'blue'
    return elements['elements']


def match_algo(term, algo, run_id):
    if term == 'FlyBase Gene Diffusion':
        if algo == 'omicsintegrator1':
            return produce_elements(grab_data(flybase_parent_dir, oi1_dict, oi1_runs, algo, run_id), None)
        elif algo == 'omicsintegrator2':
            return produce_elements(grab_data(flybase_parent_dir, oi2_dict, oi2_runs, algo, run_id), None)
        elif algo == 'rwr':
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
            if run_id == 'mode1':
                selected_edges = edges
            elif run_id == 'mode2':
                for edge in edges:
                    if edge[0] in source_gene and edge[1] in linker_nodes:
                        selected_edges.append(edge)
            elif run_id == 'mode3':
                for edge in edges:
                    if edge[0] in source_gene and edge[1] in (linker_nodes - source_gene):
                        selected_edges.append(edge)
            G = nx.Graph()
            G.add_edges_from(selected_edges)
            print(f"Using RWR algorithm with {run_id}...")
            print(f"Number of nodes: {len(G.nodes)}")
            print(f"Number of edges: {len(G.edges)} \n")
            return produce_elements(G, None)
        else:
            return {}
    elif term == 'cell-cell fusion':
        if algo == 'pathlinker':
            return produce_elements(grab_data(parent_dir,cell_pathlinker_dict, cell_pathlinker_runs, algo, run_id), cell_target_gene)
        elif algo == 'rwr':
            return produce_elements(grab_data(parent_dir, cell_rwr_dict, cell_rwr_runs, algo, run_id), cell_target_gene)
        else:
            return {}
    elif term == 'muscle cell development':
        if algo == 'pathlinker':
            return produce_elements(grab_data(parent_dir,muscle_pathlinker_dict, muscle_pathlinker_runs, algo, run_id), muscle_target_gene)
        elif algo == 'rwr':
            return produce_elements(grab_data(parent_dir,muscle_rwr_dict, muscle_rwr_runs, algo, run_id), muscle_target_gene)
        else:
            return {}
    else:
        return {}

generate_nodes(source_gene_file, source_gene)
generate_nodes(receptor_gene_file, receptor_gene)
generate_nodes(cell_target_gene_file, cell_target_gene)
generate_nodes(muscle_target_gene_file, muscle_target_gene)

update_algo_list(algo_list_flybase, oi1_dict, oi1_runs, 'omicsintegrator1')
update_algo_list(algo_list_flybase, oi2_dict, oi2_runs, 'omicsintegrator2')
update_algo_list(algo_list_cell, cell_pathlinker_dict,
                 cell_pathlinker_runs, 'pathlinker')
update_algo_list(algo_list_cell, cell_rwr_dict, cell_rwr_runs, 'rwr')
update_algo_list(algo_list_muscle, muscle_pathlinker_dict,
                 muscle_pathlinker_runs, 'pathlinker')
update_algo_list(algo_list_muscle, muscle_rwr_dict, muscle_rwr_runs, 'rwr')

app.layout = html.Div([
    dcc.Dropdown(
        id='dropdown-update-term',
        clearable=True,
        options=[
            {'label': name, 'value': name}
            for name in ['FlyBase Gene Diffusion', 'cell-cell fusion', 'muscle cell development']
        ]
    ),


    dcc.Dropdown(
        id='dropdown-update-algo',
        clearable=True,
        options=[
            {'label': name, 'value': name}
            for name in algo_list_cell
        ]
    ),


    dcc.Dropdown(
        id='dropdown-update-layout',
        value='preset',
        clearable=False,
        options=[
            {'label': name.capitalize(), 'value': name}
            for name in ['preset', 'grid', 'random', 'circle', 'cose', 'concentric']
        ]
    ),

    dcc.Markdown(id='cytoscape-selectedNodeData-markdown'),

    cyto.Cytoscape(
        id='flybase-cytoscape',
        layout={'name': 'preset'},
        style={'width': '100%', 'height': '960px'},
        elements=match_algo('cell-cell fusion', 'pathlinker',
                            str(cell_pathlinker_runs[0]).split('-')[-1]),
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
                'selector': '.green',
                'style': {
                    'background-color': 'green',
                    'line-color': 'green'
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
    )
])


@callback(Output('cytoscape-selectedNodeData-markdown', 'children'),
          Input('flybase-cytoscape', 'selectedNodeData'))
def displaySelectedNodeData(data_list):
    if data_list is None:
        return "No node selected."

    node_list = [data['label'] for data in data_list]
    return "You selected the following nodes: " + "\t ".join(node_list)


@callback(Output('flybase-cytoscape', 'layout'),
          Input('dropdown-update-layout', 'value'))
def update_layout(layout):
    return {
        'name': layout,
        'animate': True
    }


@callback(Output('dropdown-update-algo', 'options'),
          Input('dropdown-update-term', 'value'))
def update_algo_options(term):
    if term == 'cell-cell fusion':
        return algo_list_cell
    elif term == 'muscle cell development':
        return algo_list_muscle
    else:
        return algo_list_flybase


@callback(Output('flybase-cytoscape', 'elements'),
          Input('dropdown-update-algo', 'value'),
          Input('dropdown-update-term', 'value'))
def update_algo(algo, term):
    if algo is None or term is None:
        return {}
    algo_list = algo.strip().split(' ')
    return match_algo(term, algo_list[0], algo_list[1])


if __name__ == '__main__':
    app.run_server(debug=True, port=8050)
