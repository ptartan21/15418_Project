// Sequential Implementation for Finding Strongly Connected Components
#include <chrono>
#include <vector>
#include <unordered_set>

#include "../util/graph.h"
#include "../bfs/bfs.h"

// computes the set of forward reachable vertices from vi
std::unordered_set<int> forward_reachability(Graph &g, std::unordered_set<int> &S, int vi, int method) {
    std::unordered_set<int> fr_vertices;
    if (method == 0) {
        fr_vertices = bfs_bottom_up_seq(g, S, vi);
    } else if (method == 1) {
        fr_vertices = bfs_top_down_seq(g, S, vi);
    }
    return fr_vertices;
}

// flip the edges in the graph
std::unordered_set<int> backward_reachability(Graph &g, std::unordered_set<int> &S, int vi, int method) {
    std::unordered_set<int> bw_vertices;
    if (method == 0) {
        bw_vertices = bfs_bottom_up_seq(g, S, vi);
    } else if (method == 1) {
        bw_vertices = bfs_top_down_seq(g, S, vi);
    }
    return bw_vertices;
}

/***** SEQUENTIAL *****/
// stores all strongly connected components into all_scc
void compute_scc_seq(std::vector<std::unordered_set<int>> &all_scc, Graph &g, int method) {
    auto start_time = std::chrono::steady_clock::now();
    Graph rev_g = reverse_graph(g);
    std::vector<std::unordered_set<int>> V;
    std::unordered_set<int> initial_v;
    std::unordered_set<int> S;
    for (int vid = 0; vid < g->n; ++vid) {
        initial_v.insert(vid);
    }
    V.push_back(initial_v);
    for (int vid = 0; vid < g->n; ++vid) {
        // std::cout << "VID: " << vid << std::endl;
        for (auto &s : V) {
            if (s.count(vid)) {
                S = s;
                std::unordered_set<int> r_plus  = forward_reachability(g, S, vid, method);
                std::unordered_set<int> r_minus = backward_reachability(rev_g, S, vid, method);
                std::unordered_set<int> v_scc   = set_i(r_plus, r_minus);
                V.erase(std::remove(V.begin(), V.end(), S), V.end());
                V.push_back(set_d(r_plus , v_scc));
                V.push_back(set_d(r_minus, v_scc));
                std::unordered_set<int> both = set_u(r_plus, r_minus);
                V.push_back(set_d(S, both));
                all_scc.push_back(v_scc);
                break;
            }
        }
    }
    auto end_time = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    std::cout << "Number of Strongly Connected Components: " << all_scc.size() << std::endl;
}
