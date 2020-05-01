import networkx as nx
import matplotlib.pyplot as plt
import sys
import random

def export_graph(G, filename):
    # adj_list = G.adjacency()
    adj_list = nx.generate_adjlist(G)
    n = G.number_of_nodes()
    m = G.number_of_edges()
    with open(filename, "w+") as f:
        f.write("%d %d" % (n, m))
        f.write("\n")
        '''
        for row in adj_list:
            out = [row[0]] + list(row[1].keys())
            out = map(lambda x: str(x), out)
            f.write(" ".join(out))
            f.write("\n")
        '''
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

def export_random_graph_batch(n, m, batch_size, prefix, directed=False):
    for i in range(batch_size):
        print("Generating graph %d" % i)
        G = nx.gnm_random_graph(n, m, directed)
        gname = prefix + "random_graph_%d_%d_%d" % (n, m, i)
        export_graph(G, gname + ".txt")

# Generates random undirected graph with power law degree distribution and
# approximate average clustering.
#     n - number of vertices
#     m - number of edges
#     p - probability of adding a triangle after adding a random edge
def export_powerlaw_graph(n, m, p):
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
    
def export_random_clustered_graph(joint_deg_seq):
    G = nx.random_clustered_graph(joint_deg_seq)
    G = nx.Graph(G) # remove parallel edges
    G.remove_edges_from(nx.selfloop_edges(G)) # remove self-loops
    gname = "random_clustered_%d" % (len(joint_deg_seq))
    export_graph(G, gname + ".txt")
    return G, gname + ".png"

def export_random_tree(n):
    G = nx.random_tree(n)
    gname = "random_tree_%d" % n
    export_graph(G, gname + ".txt")
    return G, gname + ".png"

def export_internet_graph(n):
    G = nx.random_internet_as_graph(n)
    gname = "random_internet_%d" % n
    export_graph(G, gname + ".txt")
    return G, gname + ".png"

def export_internet_graph_batch(n, batch_size, prefix):
    for i in range(batch_size):
        print("Generating graph %d" % i)
        G = nx.random_internet_as_graph(n)
        gname = prefix + "random_internet_%d_%d" % (n, i)
        export_graph(G, gname + ".txt")

# 
#     n - number of vertices
#     m - number of random edges to add per vertex
#     p - probability of adding a triangle after adding a random edge
def export_powerlaw_cluster_graph(n, m, p):
    G = nx.powerlaw_cluster_graph(n, m, p)
    gname = "random_powerlaw_cluster_%d_%d_%d" % (n, m, p)
    export_graph(G, gname + ".txt")
    return G, gname + ".png"

def export_powerlaw_cluster_graph_batch(n, m, p, batch_size, prefix):
    for i in range(batch_size):
        print("Generating graph %d" % i)
        G = nx.powerlaw_cluster_graph(n, m, p)
        gname = prefix + "random_powerlaw_cluster_%d_%d_%s_%d" % (n, m, str(p), i)
        export_graph(G, gname + ".txt")

def export_watts_strogatz_graph(n, k, p, prefix):
    G = nx.watts_strogatz_graph(n, k, p)
    gname = prefix + "watts_strogatz_%d_%d_%s" % (n, k, str(p))
    export_graph(G, gname + ".txt")
    return G, gname + ".png"

def export_barabasi_albert_graph(n, m, prefix):
    G = nx.barabasi_albert_graph(n, m)
    gname = prefix + "barabasi_albert_%d_%d" % (n, m)
    export_graph(G, gname + ".txt")
    return G, gname + ".png"

if __name__ == "__main__":
    if (len(sys.argv) == 3):
        n = int(sys.argv[1])
        m = int(sys.argv[2])
        G, gname = export_random_graph(n, m)
        # nx.draw(G, with_labels = True)
        # plt.savefig(gname)
        # plt.show()
    else:
        # print("Usage: python3 graph_gen.py [n] [p]")
        # joint_deg_seq = list()
        # n = 100000
        # max_degree = 30
        # max_triangle_degree = 30
        # for i in range(n):
        #     deg = random.randrange(0, max_degree)
        #     triangle_deg = random.randrange(0, max_triangle_degree)
        #     joint_deg_seq.append((deg, triangle_deg))
        # G, gname = export_random_clustered_graph(joint_deg_seq)
        # export_random_tree(100000)
        # export_internet_graph_batch(10000, 50, "../graphs/random_internet/")
        # export_random_graph(50, 150)
        # export_powerlaw_cluster_graph_batch(20000, 5, 0.25, 50, "../graphs/powerlaw_cluster/20000_5_0.25/")
        # export_powerlaw_cluster_graph(20000, 5, 0.25)
        # export_erdos_renyi_graph(20000, 0.25)
        # export_random_graph(20000, 40000)
        # export_random_graph_batch(20000, 100000, 50, "../graphs/random_graph/20000_40000/")
        export_watts_strogatz_graph(20000, 20, 0.1, "../graphs/watts_strogatz/20000/")
        export_barabasi_albert_graph(20000, 20, "../graphs/barabasi_albert/20000/")
