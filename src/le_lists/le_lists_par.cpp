#include "../util/graph.h"
#include "../bfs/bfs.h"

#include <limits>



void inline le_lists_bfs_top_down_step(Graph g, vertex_set *frontier, vertex_set *next_frontier, 
    int *distances, int &num_frontier_edges, int &num_edges_checked) {
    
    int local_num_frontier_edges = 0;
    int local_num_edges_checked = 0;
    #pragma omp parallel private(local_num_frontier_edges, local_num_edges_checked)
    {
        vertex_set *local_frontier = (vertex_set *) malloc(sizeof(vertex_set));
        init_vertex_set(local_frontier, g->n);
        // Iterate over frontier
        #pragma omp for schedule(static)
        for (int i = 0; i < frontier->num_vertices; ++i) {
            int vid = frontier->vertices[i];
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                // Use CAS to guarantee nid is only added to the frontier once
                if (__sync_bool_compare_and_swap(&distances[nid], UNVISITED, distances[vid]+1)) {
                    local_frontier->vertices[local_frontier->num_vertices++] = nid;
                    local_num_frontier_edges += g->out_offsets[nid+1] - g->out_offsets[nid];
                }
                local_num_edges_checked++;
            }
        }
        // Copy local next frontier to real next frontier
        #pragma omp critical
        {
            std::memcpy(next_frontier->vertices+next_frontier->num_vertices, local_frontier->vertices, local_frontier->num_vertices*sizeof(int));
            next_frontier->num_vertices += local_frontier->num_vertices;
            num_frontier_edges += local_num_frontier_edges;
            num_edges_checked += local_num_edges_checked;
        }
        free_vertex_set(local_frontier);
    }
}

void inline le_lists_bfs_bottom_up_step(Graph g, int &frontier_size, int iter, int *distances) {
    int shared_frontier_size = 0;
    #pragma omp parallel
    {
        // Iterate over vertices
        #pragma omp for reduction(+:shared_frontier_size) schedule(static)
        for (int vid = 0; vid < g->n; ++vid) {
            if (is_unvisited(vid, distances)) {
                for (int eid = g->in_offsets[vid]; eid < g->in_offsets[vid+1]; ++eid) {
                    int nid = g->in_edge_list[eid];
                    if (distances[nid] == iter-1) {
                        distances[vid] = iter;
                        shared_frontier_size++;
                        break;
                    }
                }
            }
        }
    }
    frontier_size = shared_frontier_size;
}

void inline le_lists_bfs_hybrid(Graph g, int source, int *deltas, vertex_set* S, int *distances) {

}

void le_lists_par(Graph g, std::vector<std::vector<int>> &L_v, std::vector<std::vector<int>> & L_d,
    std::unordered_map<std::string, double> &metrics) {
    auto start_time = std::chrono::steady_clock::now();
    // std::vector<int> deltas = std::vector<int>(g->n,std::numeric_limits<int>::max());

    int *deltas = (int *) malloc(g->n * sizeof(int));
    memset(deltas, std::numeric_limits<int>::max(), g->n*sizeof(int));

    L_v = std::vector<std::vector<int>>(g->n, std::vector<int>()); // v_i's
    L_d = std::vector<std::vector<int>>(g->n, std::vector<int>()); // d(v_i, u)'s



    for (int vid = 0; vid < g->n; ++vid) {
        // std::vector<int> S;
        // std::vector<int> distances = std::vector<int>(g->n, UNVISITED);
        
        vertex_set *S = (vertex_set *) malloc(sizeof(vertex_set));
        init_vertex_set(S, g->n);
        int *distances = (int *) malloc(g->n*sizeof(int));
        memset(distances, UNVISITED, g->n*sizeof(int));
        
        le_lists_bfs_hybrid(g, vid, deltas, S, distances);
        for (int i = 0; i < S->num_vertices; ++i) {
            
        }

        // S = {u in V | d(v_i, u) < delta(u)}
        le_lists_bfs_seq(g, vid, deltas, S, distances);
        for (auto &u : S) {
            deltas[u] = distances[u];
            L_v[u].push_back(vid);
            L_d[u].push_back(distances[u]);
        }

        free_vertex_set(S);
    }
    auto end_time = std::chrono::steady_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    metrics.insert(std::make_pair("runtime", runtime));

    free(deltas);
}