from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
from pathlib import Path
import networkx as nx

graphspace = GraphSpace('liue@reed.edu', 'Lyc20010206!!!')

source_gene_file = Path('./input/source_gene.txt')
receptor_gene_file = Path('./input/receptor.txt')
parent_dir = Path("./output/amigo2")
flybase_parent_dir = Path("./output/flybase")
interactome_file = Path('./input/flybase_interactome.txt')

cell_target_gene_file = Path('./input/cell_target_gene.txt')
muscle_target_gene_file = Path('./input/muscle_target_gene.txt')

source_gene = set()
receptor_gene = set()
cell_target_gene = set()
muscle_target_gene = set()

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
            
generate_nodes(source_gene_file, source_gene)
generate_nodes(receptor_gene_file, receptor_gene)
generate_nodes(cell_target_gene_file, cell_target_gene)
generate_nodes(muscle_target_gene_file, muscle_target_gene)

def update_graphspace_graph(process, algo, run_id, source_size : int, intermediate_size : int):
    print("Logged in to GraphSpace.")
    print("Loading data...")
    if algo == 'omicsintegrator1':
        test_file = Path(flybase_parent_dir / f'{process}-{algo}-params-{run_id}')
        analysis_file = Path('./downstream_analysis/omicsintegrator1.txt')
        source_analysis_file = Path('./downstream_analysis/source_omicsintegrator1.txt')
    else:
        test_file = Path(parent_dir / f'{process}-{algo}-params-{run_id}')
        analysis_file = Path(f'./downstream_analysis/{process}-{algo}.txt')
        source_analysis_file = Path(f'./downstream_analysis/{process}-{algo}-source.txt')
    
    node_dict = {}
    with open(analysis_file, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            line = line.strip()
            endpoints = line.split("\t")
            node_dict[endpoints[0]] = [endpoints[2]]
            
    with open(source_analysis_file, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            line = line.strip()
            endpoints = line.split("\t")
            node_dict[endpoints[0]] = [endpoints[2]]
            
    edge_dict = {}
    with open(interactome_file, 'r') as f:
        for line in f:
            line = line.strip()
            endpoints = line.split("\t")
            edge_dict[(endpoints[0], endpoints[1])] = endpoints[2]
            edge_dict[(endpoints[1], endpoints[0])] = endpoints[2]
    
    G = nx.Graph()
    G.add_edges_from(edge_dict.keys())
    for key in node_dict.keys():
        node_dict[key].append(G.degree[key])
    
    target_set = None
    if process == "cell-cell-fusion":
        target_set = cell_target_gene
    elif process == "muscle-development":
        target_set = muscle_target_gene
        
    node_set = set()
    edge_list = []
    with open(test_file / 'pathway.txt') as f:
        for line in f:
            line = line.strip()
            endpoints = line.split("\t")
            node_set.add(endpoints[0])
            node_set.add(endpoints[1])
            edge_list.append((endpoints[0], endpoints[1], {"weight" :float(edge_dict[(endpoints[0], endpoints[1])])}))
            
    G = nx.Graph()
    G.add_edges_from(edge_list)
    T = nx.maximum_spanning_tree(G, weight='weight', algorithm='kruskal')
    edge_list = T.edges()
    
    G = GSGraph()            
    for node in node_set:
        if node in source_gene:    
            G.add_node(node, label=node, popup = 'Gene: ' + node + '; ' + 'Type: ' + 'source node' + '; ' + 'Node Frequency: ' + node_dict[node][0] + '; ' + 'Degree: ' + str(node_dict[node][1]))
            G.add_node_style(node, shape='rectangle', color='red', border_color = 'red', width=source_size * float(node_dict[node][0]), height=source_size * float(node_dict[node][0]))
        elif target_set is not None and node in target_set:
            G.add_node(node, label=node, popup = 'Gene: ' + node + '; ' + 'Type: ' + 'target node' + '; ' + 'Node Frequency: ' + node_dict[node][0] + '; ' + 'Degree: ' + str(node_dict[node][1]))
            G.add_node_style(node, shape='rectangle', color='green', border_color = 'green', width=source_size * float(node_dict[node][0]), height=source_size * float(node_dict[node][0]))
        else:
            if node in receptor_gene:
                G.add_node(node, label=node, popup = 'Gene: ' + node + '; ' + 'Type: ' + 'receptor node' + '; ' + 'Node Frequency: ' + node_dict[node][0] + '; ' + 'Degree: ' + str(node_dict[node][1]))
                G.add_node_style(node, shape='triangle', color='blue', border_color= 'blue', width=intermediate_size * float(node_dict[node][0]), height=intermediate_size * float(node_dict[node][0])) 
            else:
                G.add_node(node, label=node, popup = 'Gene: ' + node + '; ' + 'Type: ' + 'intermediate node' + '; ' + 'Node Frequency: ' + node_dict[node][0] + '; ' + 'Degree: ' + str(node_dict[node][1]))
                G.add_node_style(node, shape='ellipse', color='grey', border_color = 'grey', width=intermediate_size * float(node_dict[node][0]), height=intermediate_size * float(node_dict[node][0]))

    for edge in edge_list:
        G.add_edge(edge[0], edge[1], directed=False, popup='Edge weight: ' + edge_dict[(edge[0], edge[1])])
        G.add_edge_style(edge[0], edge[1], width= 4 * float(edge_dict[(edge[0], edge[1])]))
        if edge[0] in source_gene and edge[1] in source_gene:
            G.add_edge_style(edge[0], edge[1], color='red', width= 4 * float(edge_dict[(edge[0], edge[1])]))
        elif target_set is not None and edge[0] in target_set and edge[1] in target_set:
            G.add_edge_style(edge[0], edge[1], color='green', width= 4 * float(edge_dict[(edge[0], edge[1])]))
        
    G.set_name(f'{process}-{algo}-{run_id}')
    G.set_tags([f'{process}', f'{algo}'])
    G.set_data(data={'description': f'Using {algo} for {process} analysis, run_id: {run_id}'})
    graph = graphspace.update_graph(G)
    print(graph.id)
    print("Done.")
    
# update_graphspace_graph('FlyBase', 'omicsintegrator1', '5E6I7VA', 60, 120)
# update_graphspace_graph('cell-cell-fusion', 'pathlinker', '6SWY7JS', 70, 70)
# update_graphspace_graph('muscle-development', 'pathlinker', '6SWY7JS', 70, 70)
update_graphspace_graph('cell-cell-fusion', 'rwr', 'Y5BUYRK', 70, 70)
update_graphspace_graph('muscle-development', 'rwr', 'PPOMXB3', 70, 70)