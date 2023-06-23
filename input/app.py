from dash import Dash, html
import dash_cytoscape as cyto
import networkx as nx


app = Dash(__name__)

G = nx.Graph()
G.add_edge('a', 'b')
G.add_edge('b', 'c')
G.add_edge('c', 'd')

elements = nx.cytoscape_data(G)
for i in range(len(elements['elements']['nodes'])):
    elements['elements']['nodes'][i]['data']['label'] = elements['elements']['nodes'][i]['data']['id']
    elements['elements']['nodes'][i]['position'] = {'x': 40 * i, 'y': 40 * i}
    
# nodes = [
#     {
#         'data': {'id': short, 'label': label},
#         'position': {'x': 20 * lat, 'y': -20 * long}
#     }
#     for short, label, long, lat in (
#         ('la', 'Los Angeles', 34.03, -118.25),
#         ('nyc', 'New York', 40.71, -74),
#         ('to', 'Toronto', 43.65, -79.38),
#         ('mtl', 'Montreal', 45.50, -73.57),
#         ('van', 'Vancouver', 49.28, -123.12),
#         ('chi', 'Chicago', 41.88, -87.63),
#         ('bos', 'Boston', 42.36, -71.06),
#         ('hou', 'Houston', 29.76, -95.37)
#     )
# ]

# edges = [
#     {'data': {'source': source, 'target': target}}
#     for source, target in (
#         ('van', 'la'),
#         ('la', 'chi'),
#         ('hou', 'chi'),
#         ('to', 'mtl'),
#         ('mtl', 'bos'),
#         ('nyc', 'bos'),
#         ('to', 'hou'),
#         ('to', 'nyc'),
#         ('la', 'nyc'),
#         ('nyc', 'bos')
#     )
# ]

# elements = nodes + edges

app.layout = html.Div([
    cyto.Cytoscape(
        id='cytoscape-layout-1',
        elements=elements['elements'],
        style={'width': '100%', 'height': '960px'},
        layout={
            'name': 'preset'
        }
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)
