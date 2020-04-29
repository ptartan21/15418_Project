import numpy as np
import networkx as nx

def get_stats(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
        n = int(lines[0].split(" ")[0])
        m = int(lines[0].split(" ")[1])
        print("n:", n)
        print("m:", m)
        max_degree = 0
        for i in range(1, len(lines)):
            line = lines[i]
            line_split = line.split(" ")
            out_degree = len(line_split)-1
            max_degree = max(max_degree, out_degree)
        print("Average degree:", m/n)
        print("Max degree:", max_degree)
    G = construct_graph(filepath)
    print("Diameter:", nx.diameter(G))

def construct_graph(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
        G = nx.Graph()
        for i in range(1, len(lines)):
            line = lines[i]
            line_split = line.split(" ")
            vid = int(line_split[0])
            for j in range(1, len(line_split)):
                nid = int(line_split[j])
                G.add_edge(vid, nid)
        return G


get_stats("powerlaw_cluster/20000_5_0.25/random_powerlaw_cluster_20000_5_0.25_0.txt")
construct_graph("powerlaw_cluster/20000_5_0.25/random_powerlaw_cluster_20000_5_0.25_0.txt")