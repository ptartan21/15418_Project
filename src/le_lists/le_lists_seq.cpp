
#include "../util/graph.h"
#include "../bfs/bfs.h"

#include <limits>

void inline le_lists_bfs_seq(Graph g, int source, std::vector<int> deltas,
    std::vector<int> &S, std::vector<int> &distances) {
    std::vector<int> frontier;
    std::vector<int> next_frontier;
    distances[source] = 0;
    frontier.push_back(source);
    while (frontier.size() > 0) {
        next_frontier.clear();
        for (auto &vid : frontier) {
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                if (distances[nid] == UNVISITED) {
                    distances[nid] = distances[vid]+1;
                    if (distances[nid] < deltas[nid]) {
                        S.push_back(nid);
                        // Only need to explore S and its outgoing edges
                        next_frontier.push_back(nid);
                    }
                }
            }
        }
        std::swap(frontier, next_frontier);
    }
}

void le_lists_seq(Graph g, std::vector<std::vector<int>> &L_v, std::vector<std::vector<int>> &L_d,
    std::unordered_map<std::string, double> &metrics) {
    auto start_time = std::chrono::steady_clock::now();
    std::vector<int> deltas = std::vector<int>(g->n,std::numeric_limits<int>::max());
    L_v = std::vector<std::vector<int>>(g->n, std::vector<int>()); // v_i's
    L_d = std::vector<std::vector<int>>(g->n, std::vector<int>()); // d(v_i, u)'s
    for (int vid = 0; vid < g->n; ++vid) {
        std::vector<int> S;
        std::vector<int> distances = std::vector<int>(g->n, UNVISITED);
        // S = {u in V | d(v_i, u) < delta(u)}
        le_lists_bfs_seq(g, vid, deltas, S, distances);
        for (auto &u : S) {
            deltas[u] = distances[u];
            L_v[u].push_back(vid);
            L_d[u].push_back(distances[u]);
        }
    }
    auto end_time = std::chrono::steady_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    metrics.insert(std::make_pair("runtime", runtime));
}