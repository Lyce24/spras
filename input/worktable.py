import dash_cytoscape as cyto
import matplotlib.pyplot as plt
import networkx as nx
from dash import Dash, html

app = Dash(__name__)

interactome = 'interactome-flybase-collapsed-weighted.txt'
flybase_interactome = 'flybase_interactome.txt'
source_gene_file = 'source_gene.txt'
flybase_dict = 'flybase_dict.txt'
myoblast_fusion_components = 'myoblast_fusion_components.txt'
gene_dict = {}
selected_edges_file = 'selected_edges.txt'

index = 0
with open(interactome, 'r') as f:
    with open(flybase_interactome, 'w') as f2:
        with open(flybase_dict, 'w') as f3:
            for line in f:
                if index == 0:
                    index += 1
                    continue
                # split the line by tabs
                temp = line.strip().split('\t')
                if temp[0] == 'nan':
                    print(temp)
                    temp[0] = 'TRPV'
                if temp[1] == 'nan':
                    print(temp)
                    temp[1] = 'TRPV'
                if temp[4] not in gene_dict:
                    gene_dict[temp[4]] = temp[0]
                    f3.write(temp[4] + '\t' + temp[0] + '\n')
                if temp[5] not in gene_dict:
                    gene_dict[temp[5]] = temp[1]
                    f3.write(temp[5] + '\t' + temp[1] + '\n')
                f2.write(temp[0] + '\t' + temp[1] + '\t' + temp[2] + '\n')

source_gene = set()

with open(myoblast_fusion_components, 'r') as f:
    with open(source_gene_file, 'w') as f2:
        f2.write('NODEID\tprize\n')
        for line in f:
            # split the line by tabs
            temp = line.strip().split('\t')
            source_gene.add(gene_dict[temp[1]])
            f2.write(gene_dict[temp[1]] + '\t' + '1' + '\n')

nodes = set()
edges = []
with open(flybase_interactome, 'r') as edges_f:
    for i, line in enumerate(edges_f):
         # if the first line is the title, skip it
        if i == 0:
            continue
        line = line.strip()
        endpoints = line.split("\t")
        if len(endpoints) != 3:
            raise ValueError(f"Edge {line} does not contain 2 nodes separated by ' ' and a weight")
        nodes.add(endpoints[0])
        nodes.add(endpoints[1])
        edges.append((endpoints[0], endpoints[1], {'weight' : endpoints[2]}))

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

selected_edges = []
for edge in edges:
    if edge[0] in linker_nodes and edge[1] in linker_nodes:
        selected_edges.append(edge)

with open(selected_edges_file, 'w') as f:
    for i in selected_edges:
        f.write(i[0] + '\t' + i[1] + '\t' + i[2]['weight'] + '\n')

print(len(linker_nodes))
print(len(selected_edges))

P = nx.Graph()
P.add_nodes_from(linker_nodes)
P.add_edges_from(selected_edges)

pos = nx.spring_layout(P)
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
plt.savefig("myoblast_fusion_components.png") # save as png



elements = nx.cytoscape_data(P)
for i in range(len(elements['elements']['nodes'])):
    elements['elements']['nodes'][i]['data']['label'] = elements['elements']['nodes'][i]['data']['id']
    if elements['elements']['nodes'][i]['data']['id'] in source_gene:
        elements['elements']['nodes'][i]['classes'] = 'red'
        elements['elements']['nodes'][i]['position'] = {'x': 400, 'y': 40 * i}
    else:
        elements['elements']['nodes'][i]['position'] = {'x': 3000, 'y': 40 * i}


app.layout = html.Div([
    cyto.Cytoscape(
        id='cytoscape-styling-1',
        layout={'name': 'preset'},
        style={'width': '100%', 'height': '960px'},
        elements=elements['elements'],
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
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)
