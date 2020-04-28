// Sequential Implementation for Finding Strongly Connected Components
#include <unordered_set>
#include <algorithm>
#include "../util/graph.h"

std::unordered_set<int> intersection(std::unordered_set<int> &S1, std::unordered_set<int> &S2) {
    // implement set intersection
}

std::unordered_set<int> difference(std::unordered_set<int> &S1, std::unordered_set<int> &S2) {
    // implement set difference, S1 - S2
}

// computes the set of forward reachable vertices from vi
std::unordered_set<int> forward_reachability(/* */, int vi) {
    std::unordered_set<int> fr_vertices;
    // implement a breadth-first search to find all forward reachable vertices in S from vi
    return fr_vertices;
}

std::unordered_set<int> backward_reachability(/* */, int vi) {
    std::unordered_set<int> bw_vertices;
    // implement a symmetrical BFS to find all backward reachable vertices in S from vi
    return bw_vertices;
}

void compute_scc(std::vector<std::unordered_set<int>> &all_scc, Graph &g) {
    // initialize V
    for (int vi = 0; vi < g->n; ++vi) {
        // construct S as the subgraph of V that contains vi
        auto r_plus  = forward_reachability(g, S, vi);
        auto r_minus = backward_reachability(g, S, vi);
        std::unordered_set<int> v_scc = intersection(r_plus, r_minus);
        // update V
        all_scc.push_back(v_scc);
    }
}
