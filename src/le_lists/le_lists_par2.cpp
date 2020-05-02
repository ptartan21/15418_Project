#include "../util/graph.h"
#include "../bfs/bfs.h"

#include <limits>

void inline le_lists_bfs(Graph g, int source, std::vector<int> deltas,
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

bool comparator(const std::pair<int, int> &a, const std::pair<int, int> &b) {
    // Sort by increasing vertex ID
    return a.first < b.first;
}

void le_lists_par2(Graph g, std::vector<std::vector<std::pair<int, int>>> &L,
    std::unordered_map<std::string, double> &metrics) {

    auto start_time = std::chrono::steady_clock::now();
    L = std::vector<std::vector<std::pair<int, int>>>(g->n, std::vector<std::pair<int, int>>());
    std::vector<int> deltas = std::vector<int>(g->n,std::numeric_limits<int>::max());
    
    #pragma omp parallel shared(L)
    {
        #pragma omp for schedule(dynamic)
        for (int vid = 0; vid < g->n; ++vid) {
            std::vector<int> S;
            std::vector<int> distances = std::vector<int>(g->n, UNVISITED);
            // S = {u in V | d(v_i, u) < delta(u)}
            le_lists_bfs(g, vid, deltas, S, distances);
            for (auto &u : S) {
                #pragma omp critical
                {
                    deltas[u] = distances[u];
                    std::pair<int, int> vid_dist = std::pair<int, int>(vid, distances[u]);
                    L[u].push_back(vid_dist);
                }
            }
        }
        #pragma omp for schedule(static)
        for (int vid = 0; vid < g->n; ++vid) {
            if (L[vid].size() > 0) {
                std::sort(L[vid].begin(), L[vid].end(), comparator);
                int dist_threshold = L[vid][0].second;
                std::vector<std::pair<int, int>> new_L_vid;
                
                new_L_vid.push_back(L[vid][0]);
                for (int i = 1; i < L[vid].size(); ++i) {
                    if (L[vid][i].second < dist_threshold) {
                        new_L_vid.push_back(L[vid][i]);
                        dist_threshold = L[vid][i].second;
                    }
                }
                L[vid] = new_L_vid;
            }
        }
    }
    
    auto end_time = std::chrono::steady_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    metrics.insert(std::make_pair("runtime", runtime));
}