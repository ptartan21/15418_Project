// Hybrid Implementation for Finding Strongly Connected Components
#include <chrono>
#include <vector>
#include <unordered_set>

#include "../util/graph.h"
#include "../bfs/bfs.h"

// computes the set of forward reachable vertices from vi
void forward_reachability_hybrid(Graph &g, unsigned char *S, unsigned char *reach, int vi) {
    bfs_hybrid(g, S, reach, vi);
}

// flip the edges in the graph
void backward_reachability_hybrid(Graph &g, unsigned char *S, unsigned char *reach, int vi) {
    bfs_hybrid(g, S, reach, vi);
}

/***** SEQUENTIAL *****/
// stores all strongly connected components into all_scc
void compute_scc_hybrid(std::vector<unsigned char *> &all_scc, Graph &g) {
    auto start_time = std::chrono::steady_clock::now();
    int n = g->n;
    Graph rev_g = reverse_graph(g);
    std::vector<unsigned char *> V;
    unsigned char *initial_v = (unsigned char *) calloc(n, sizeof(unsigned char));
    memset(initial_v, 1, g->n * sizeof(unsigned char));
    V.push_back(initial_v);
    unsigned char *S;
    for (int vid = 0; vid < n; ++vid) {
        // std::cout << "VID: " << vid << std::endl;
        for (auto &s : V) {
            if (s[vid]) {
                S = s;
                break;
            }
        }
        if (S) {
            unsigned char *r_plus    = (unsigned char *) calloc(n, sizeof(unsigned char));
            unsigned char *r_minus   = (unsigned char *) calloc(n, sizeof(unsigned char));
            unsigned char *v_scc     = (unsigned char *) calloc(n, sizeof(unsigned char));
            unsigned char *s1        = (unsigned char *) calloc(n, sizeof(unsigned char));
            unsigned char *s2        = (unsigned char *) calloc(n, sizeof(unsigned char));
            unsigned char *s3        = (unsigned char *) calloc(n, sizeof(unsigned char));
            unsigned char *tmp       = (unsigned char *) calloc(n, sizeof(unsigned char)); 

            forward_reachability_hybrid(g, S, r_plus, vid);             // populate r_plus with all forward reachable vertices
            backward_reachability_hybrid(rev_g, S, r_minus, vid);       // populate r_minus with all backward reachable vertices
            set_i(r_plus , r_minus, v_scc, n);                          // v_scc = r_plus AND r_minus
            set_d(r_plus , v_scc  , s1   , n);                          // s1 = r_plus  - v_scc
            set_d(r_minus, v_scc  , s2   , n);                          // s2 = r_minus - v_scc
            set_u(r_plus , r_minus, tmp  , n);                          // tmp = r_plus + r_minus 
            set_d(S,       tmp    , s3   , n);                          // s3 = S - (r_plus + r_minus)
            V.erase(std::remove(V.begin(), V.end(), S), V.end());
            V.push_back(s1);
            V.push_back(s2);
            V.push_back(s3);
            all_scc.push_back(v_scc);
            S = NULL;
        }
    }
    auto end_time = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    std::cout << "Number of Strongly Connected Components: " << all_scc.size() << std::endl;
}