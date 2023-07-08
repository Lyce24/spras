from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
from graphspace_python.graphs.classes.gslayout import GSLayout
from pathlib import Path

graphspace = GraphSpace('liue@reed.edu', 'Lyc20010206!!!')
print("Logged in to GraphSpace.")
print("Loading data...")
interactome_file = Path('./input/flybase_interactome.txt')
source_gene_file = Path('./input/source_gene.txt')
receptor_gene_file = Path('./input/receptor.txt')
test_file = Path('./FlyBase/flybase-omicsintegrator1-params-5Q4YIBB')
analysis_file = Path('./downstream_analysis/omicsintegrator1.txt')
source_analysis_file = Path('./downstream_analysis/source_omicsintegrator1.txt')

source_gene = set()
receptor_gene = set()

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

node_dict = {}
with open(analysis_file, 'r') as f:
    for line in f:
        line = line.strip()
        endpoints = line.split("\t")
        node_dict[endpoints[0]] = endpoints[2]
        
with open(source_analysis_file, 'r') as f:
    for line in f:
        line = line.strip()
        endpoints = line.split("\t")
        node_dict[endpoints[0]] = endpoints[2]
        
edge_dict = {}
with open(interactome_file, 'r') as f:
    for line in f:
        line = line.strip()
        endpoints = line.split("\t")
        edge_dict[(endpoints[0], endpoints[1])] = endpoints[2]

G = GSGraph()
node_set = set()
edge_list = []
with open(test_file / 'pathway.txt') as f:
    for line in f:
        line = line.strip()
        endpoints = line.split("\t")
        node_set.add(endpoints[0])
        node_set.add(endpoints[1])
        edge_list.append((endpoints[0], endpoints[1]))
        
for node in node_set:
    if node in source_gene:    
        G.add_node(node, label=node, popup = 'Gene: ' + node + '; ' + 'Type: ' + 'source' + '; ' + 'Node Frequency: ' + node_dict[node])
        G.add_node_style(node, shape='ellipse', color='red', border_color = 'red', width=60 * float(node_dict[node]), height=60 * float(node_dict[node]))
    elif node in receptor_gene:
        G.add_node(node, label=node, popup = 'Gene: ' + node + '; ' + 'Type: ' + 'receptor' + '; ' + 'Node Frequency: ' + node_dict[node])
        G.add_node_style(node, shape='triangle', color='blue', border_color= 'blue', width=160 * float(node_dict[node]), height=160 * float(node_dict[node]))
    else:
        G.add_node(node, label=node, popup = 'Gene: ' + node + '; ' + 'Type: ' + 'intermediate node' + '; ' + 'Node Frequency: ' + node_dict[node])
        G.add_node_style(node, shape='ellipse', color='grey', border_color = 'grey', width=160 * float(node_dict[node]), height=160 * float(node_dict[node]))

for edge in edge_list:
    G.add_edge(edge[0], edge[1], directed=False, popup='Edge weight: ' + edge_dict[(edge[0], edge[1])])
    G.add_edge_style(edge[0], edge[1], width= 4 * float(edge_dict[(edge[0], edge[1])]))
    if edge[0] in source_gene and edge[1] in source_gene:
        G.add_edge_style(edge[0], edge[1], color='red', width= 4 * float(edge_dict[(edge[0], edge[1])]))

G.set_name('FlyBase-5Q4YIBB')
G.set_tags(['FlyBase'])
G.set_data(data={'description': 'Using Omics Integrator 1 to diffuse from the FlyBase Interactome. (Parameters - b: 2.4 / d: 10 / g: 1e-3 / mu: 0.0005 / r: 0 / w: 2)'})
graph = graphspace.update_graph(G)
print(graph.id)
print("Done.")