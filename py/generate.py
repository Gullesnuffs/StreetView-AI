"""
  Script that generates a test case from a given place.
  Uses Open Street Map data through the OSMnx python API.

  Takes as input in seperat lines the place to generate the
  graph from the total time in seconds and the number of cars
  to be included in the test case.

  problem:
    - If a street is directed we don't know the direction
    relating to the nodes it is going from and to. So we 
    make all edges bi-directional.

  improvements:
    - If street has max speed data use it to calculate the
    time to drive the street.
"""
from pprint import pprint as pp
from math import ceil
import osmnx as ox

place = input()
T = int(input())
C = int(input())

# Generate the graph using OSMnx
city = ox.graph_from_place(place, network_type='drive')
nodes = {}
node_mapping = {}
edges = {}

for e in city.edges:
    e0 = e[0]
    e1 = e[1]
    if (e0, e1) not in edges and (e1, e0) not in edges:
        nodes[e0] = (city.node[e0]['y'], city.node[e0]['x'])
        nodes[e1] = (city.node[e1]['y'], city.node[e1]['x'])
        edges[(e0,e1)] = city[e0][e1][0]

N = len(nodes) # number of junctions
M = len(edges) # number of streets
S = 0 # starting node for the cars
# First part of the test case
print(N, M, T, C, S)

# map nodes to ids range(N)
for i, n in enumerate(nodes):
    node_mapping[n] = i 

# print nodes
for n in nodes:
    print(nodes[n][0], nodes[n][1])

for e in edges:
    edge = edges[e]
    length = ceil(edge['length'])
    # divide by 1.25 to get 45 km/h speed on roads
    duration = ceil(length / 1.25)
    direction = 2 # 1 if edge['oneway'] else 2
    print(node_mapping[e[0]], node_mapping[e[1]], direction, duration, length)
