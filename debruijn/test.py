import networkx as nx

G = nx.complete_graph(4)
for path in nx.all_simple_paths(G, source=0, target=3):
    print(path)

paths = nx.all_simple_paths(G, source=0, target=3, cutoff=2)
print(list(paths))

