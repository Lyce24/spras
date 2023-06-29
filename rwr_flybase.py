from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx

flybase_interactome = Path('./input/flybase_interactome.txt')
source_gene_file = Path('./input/source_gene.txt')
output_file = Path('./FlyBase/rwr-pathway.txt')
output_file.parent.mkdir(parents=True, exist_ok=True)

source_gene = set()

nodes = set()
edges = []
with open(flybase_interactome, 'r') as edges_f:
    for line in edges_f:
         # if the first line is the title, skip it
        line = line.strip()
        endpoints = line.split("\t")
        if len(endpoints) != 3:
            raise ValueError(f"Edge {line} does not contain 2 nodes separated by ' ' and a weight")
        nodes.add(endpoints[0])
        nodes.add(endpoints[1])
        edges.append((endpoints[0], endpoints[1], {'weight' : endpoints[2]}))

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


G = nx.Graph()
G.add_nodes_from(nodes)
G.add_edges_from(edges)

personalization_vector = {}
# assigning value 1 to the source and target nodes to pass it to the random walk function
for i in source_gene:
    personalization_vector[i] = 1

pr = nx.pagerank(G, alpha=0.85, personalization=personalization_vector)
threshold = 0.001
linker_nodes = set()
for i in source_gene:
    linker_nodes.add(i)

for key in pr:
    if pr[key] > threshold:
        linker_nodes.add(key)

# '''
# Edge selection:
# 1. ignore the connecttion between linker nodes
# 2. ignore the connection between source node and source node
# '''

selected_edges = []
for edge in edges:
    if edge[0] in linker_nodes and edge[1] in linker_nodes:
        selected_edges.append(edge)

with open(output_file, 'w') as f:
    for i in selected_edges:
        f.write(i[0] + '\t' + i[1] + '\t' + i[2]['weight'] + '\n')

print(len(linker_nodes))
print(len(selected_edges))

P = nx.Graph()
P.add_nodes_from(linker_nodes)
P.add_edges_from(selected_edges)

pos = nx.random_layout(P)
options = {"edgecolors": "tab:gray", "node_size": 300, "alpha": 0.9}
nx.draw_networkx_nodes(P, pos, nodelist= linker_nodes, node_color="tab:gray", **options)
nx.draw_networkx_nodes(P, pos, nodelist= source_gene, node_color="tab:blue", **options)
# # define a layout for drawing
nx.draw_networkx_edges(P, pos, width=1.0, alpha=0.5)

    # # some math labels
labels = {}
for i in linker_nodes:
    labels[i] = i
nx.draw_networkx_labels(P, pos, labels, font_size=8, font_color="black")

plt.tight_layout()
plt.axis("off")
plt.savefig("./FlyBase/myoblast_fusion_components.png") # save as png
