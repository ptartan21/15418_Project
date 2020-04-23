import networkx as nx
import matplotlib.pyplot as plt
import sys

def export_graph(G, filename):
    adj_list = nx.generate_adjlist(G)
    n = G.number_of_nodes()
    m = G.number_of_edges()
    with open(filename, "w") as f:
        f.write("%d %d" % (n, m))
        f.write("\n")
        for row in adj_list:
            f.write(row)
            f.write("\n")
    print("Exported to %s" % filename)

# Generates G_{n,p} (Erdos-Renyi graph/binomial graph).
#     n - number of vertices
#     p - probability of edge creation
#     directed - True for directed; False for undirected
# Each of the n(n-1)/2 (undirected) or n(n-1) (directed) edges is included
# with probability p.
def export_erdos_renyi_graph(n, p, directed=False):
    G = nx.erdos_renyi_graph(n, p, directed)
    gname = "erdos_renyi_%d_%s" % (n, p)
    export_graph(G, gname + ".txt")
    return G, gname + ".png"

# Generates random G_{n,m}.
#     n - number of vertices
#     m - number of edges
#     directed - True for directed; False for undirected
# Graph is chosen uniform randomly from set of all graphs with n vertices and
# m edges.
def export_random_graph(n, m, directed=False):
    G = nx.gnm_random_graph(n, m, directed)
    gname = "random_graph_%d_%s" % (n, m)
    export_graph(G, gname + ".txt")
    return G, gname + ".png"

# Generates random undirected graph with power law degree distribution and
# approximate average clustering.
#     n - number of vertices
#     m - number of edges
#     p - probability of adding a triangle after adding a random edge
def export_power_law_graph(n, m, p):
    G = nx.powerlaw_cluster_graph(n, m, p)
    gname = "powerlaw_%d_%s" % (n, m) 
    export_graph(G, gname + ".txt")
    return G, gname + ".png" 

def export_complete_graph(n):
    G = nx.complete_graph(n)
    m = (n * (n - 1) // 2) 
    gname = "complete_%d_%s" % (n, m)
    export_graph(G, gname + ".txt")
    return G, gname + ".png"

def mapping(x):
    return x

if __name__ == "__main__":
    if (len(sys.argv) == 3):
        n = int(sys.argv[1])
        m = int(sys.argv[2])
        G, gname = export_power_law_graph(n, m, 0.3)
        nx.draw(G, with_labels = True)
        plt.savefig(gname)
        plt.show()
    else:
        print("Usage: python3 graph_gen.py [n] [p]")
