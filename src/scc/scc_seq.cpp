// Sequential Implementation for Finding Strongly Connected Components
#include "scc_util.cpp"

std::unordered_set<int> set_u(std::unordered_set<int> &S1, std::unordered_set<int> &S2) {
    std::unordered_set<int> res;
    for (auto &v : S1) { res.insert(v); }
    for (auto &v : S2) { res.insert(v); }
    return res;
}

// computes the set intersection of S1 and S2
std::unordered_set<int> set_i(std::unordered_set<int> &S1, std::unordered_set<int> &S2) {
    if (S1.size() > S2.size()) {
        return set_i(S2, S1);
    }
    std::unordered_set intersect;
    for (auto &v : S1) {
        if (S2.count(v)) {
            intersect.insert(v);
        }
    }
    return intersect;
}

// computes the set difference S1 - S2
std::unordered_set<int> set_d(std::unordered_set<int> &S1, std::unordered_set<int> &S2) {
    std::unordered_set<int> diff;
    for (auto &v : S1) {
        if (!S2.count(v)) {
            diff.insert(v);
        }
    }
    return diff;
}

// computes the set of forward reachable vertices from vi
std::unordered_set<int> forward_reachability(Graph &g, std::unordered_set<int> &S, int vi) {
    std::unordered_set<int> fr_vertices = bfs_bottom_up_seq(g, S, vi);
    return fr_vertices;
}

// flip the edges in the graph
std::unordered_set<int> backward_reachability(Graph &g, std::unordered_set<int> &S, int vi) {
    std::unordered_set<int> bw_vertices = bfs_bottom_up_seq(g, S, vi);
    return bw_vertices;
}

/***** SEQUENTIAL *****/
// stores all strongly connected components into all_scc
void compute_scc_seq(std::vector<std::unordered_set<int>> &all_scc, Graph &g) {
    Graph rev_g = reverse_graph(g);
    std::unordered_set<unordered_set<int>> V;
    std::unordered_set<int> initial_v;
    std::unordered_set<int> S;
    for (int vid = 0; vid < g->n; ++vid) {
        initial_v.insert(vid);
    }
    V.insert(initial_v);
    for (int vid = 0; vid < g->n; ++vid) {
        for (auto &s : V) {
            if (s.count(vid)) {
                S = s;
                break;
            }
        }
        std::unordered_set<int> r_plus  = forward_reachability(g, S, vid);
        std::unordered_set<int> r_minus = backward_reachability(rev_g, S, vid);
        std::unordered_set<int> v_scc   = set_i(r_plus, r_minus);
        V.erase(S);
        V.insert(set_d(r_plus , v_scc));
        V.insert(set_d(r_minus, v_scc));
        V.insert(set_d(S, set_u(r_plus, r_minus)));
        all_scc.push_back(v_scc);
    }
}
