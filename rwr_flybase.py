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

edge_threshold = input("Please input the threshold for edges: ")

with open(flybase_interactome, 'r') as edges_f:
    total_edges_weight = 0
    number_of_edges = 0
    for line in edges_f:
         # if the first line is the title, skip it
        line = line.strip()
        endpoints = line.split("\t")
        if len(endpoints) != 3:
            raise ValueError(f"Edge {line} does not contain 2 nodes separated by ' ' and a weight")
        nodes.add(endpoints[0])
        nodes.add(endpoints[1])
        total_edges_weight += float(endpoints[2])
        number_of_edges += 1
        if float(endpoints[2]) > float(edge_threshold):
            edges.append((endpoints[0], endpoints[1], endpoints[2]))
    avg_weight = total_edges_weight / number_of_edges
    print(f"Average weight of edges: {avg_weight}")

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
threshold = float(input("Please input the threshold for nodes: "))
linker_nodes = set()
for i in source_gene:
    linker_nodes.add(i)

for key in pr:
    if pr[key] > threshold:
        linker_nodes.add(key)

selected_edges = []
for edge in edges:
    if edge[0] in linker_nodes and edge[1] in linker_nodes:
        selected_edges.append(edge)

with open(output_file, 'w') as f:
    for i in selected_edges:
        f.write(i[0] + '\t' + i[1] + '\t' + i[2]['weight'] + '\n')

P = nx.Graph()
P.add_edges_from(selected_edges)
print("Number of nodes in the pathway: ", P.number_of_nodes())
print("Number of edges in the pathway: ", P.number_of_edges())
