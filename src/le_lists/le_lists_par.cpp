#include "../util/graph.h"
#include "../bfs/bfs.h"

#include <limits>



void inline le_lists_bfs_top_down_step(Graph g, vertex_set *frontier, vertex_set *next_frontier, 
    int *distances, int &num_frontier_edges, int &num_edges_checked, int *deltas, vertex_set *S, int iter) {
            
    int local_num_frontier_edges = 0;
    int local_num_edges_checked = 0;
    #pragma omp parallel private(local_num_frontier_edges, local_num_edges_checked)
    {
        vertex_set *local_frontier = (vertex_set *) malloc(sizeof(vertex_set));
        init_vertex_set(local_frontier, g->n);
        vertex_set *local_S = (vertex_set *) malloc(sizeof(vertex_set));
        init_vertex_set(local_S, g->n);
        // Iterate over frontier
        #pragma omp for schedule(static)
        for (int i = 0; i < frontier->num_vertices; ++i) {
            int vid = frontier->vertices[i];
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                // Use CAS to guarantee nid is only added to the frontier once
                if (__sync_bool_compare_and_swap(&distances[nid], UNVISITED, distances[vid]+1)) {
                    if (iter < deltas[nid]) {
                        local_frontier->vertices[local_frontier->num_vertices++] = nid;
                        local_num_frontier_edges += g->out_offsets[nid+1] - g->out_offsets[nid];
                        local_S->vertices[local_S->num_vertices++] = nid;
                    }

                }
                local_num_edges_checked++;
            }
        }
        // Copy local next frontier to real next frontier
        #pragma omp critical
        {
            std::memcpy(next_frontier->vertices+next_frontier->num_vertices, local_frontier->vertices, local_frontier->num_vertices*sizeof(int));
            next_frontier->num_vertices += local_frontier->num_vertices;
            std::memcpy(S->vertices+S->num_vertices, local_S->vertices, local_S->num_vertices*sizeof(int));
            S->num_vertices += local_S->num_vertices;
            num_frontier_edges += local_num_frontier_edges;
            num_edges_checked += local_num_edges_checked;
        }
        free_vertex_set(local_frontier);
        free_vertex_set(local_S);
    }
}

void inline le_lists_bfs_bottom_up_step(Graph g, int &frontier_size, int iter, int *distances, int* deltas, vertex_set *S) {
    int shared_frontier_size = 0;

    #pragma omp parallel
    {
        // Iterate over vertices
        #pragma omp for reduction(+:shared_frontier_size) schedule(static)
        for (int vid = 0; vid < g->n; ++vid) {
            if (is_unvisited(vid, distances)) {
                for (int eid = g->in_offsets[vid]; eid < g->in_offsets[vid+1]; ++eid) {
                    int nid = g->in_edge_list[eid];
                    // if nid is in the frontier
                    if (distances[nid] == iter-1) {
                        distances[vid] = iter;
                        if (distances[vid] < deltas[vid]) {
                            shared_frontier_size++;
                            #pragma omp critical
                            {
                                S->vertices[S->num_vertices++] = vid;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    frontier_size = shared_frontier_size;
}

void inline le_lists_bfs_hybrid(Graph g, int source, int *deltas, vertex_set* S, int *distances) {
    double alpha = 20.0;
    double gamma = 5.0;
    
    S->num_vertices = 0;
    memset(distances, UNVISITED, g->n*sizeof(int));

    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);

    int frontier_size = 0; // invariant: frontier_size = frontier->num_vertices
    int iter = 1;

    int num_frontier_edges = 0; // number of edges to check from frontier
    int num_unvisited_edges = 2*g->m; // number of edges to check from unvisited vertices
    unsigned char last_step = TOP_DOWN;
    // unsigned char last_step = BOTTOM_UP;

    frontier->vertices[frontier->num_vertices++] = source;
    frontier_size++;
    distances[source] = 0;
    num_frontier_edges += g->out_edge_list[source+1] - g->out_edge_list[source];

    while (frontier_size > 0) {
        if (last_step == TOP_DOWN) {
            bool should_switch = ((double) num_frontier_edges) > (((double) num_unvisited_edges) / alpha);
            should_switch = false;
            if (should_switch) {
                last_step = BOTTOM_UP;
                le_lists_bfs_bottom_up_step(g, frontier_size, iter, distances, deltas, S);
            } else {
                reset_frontier(next_frontier);
                int num_edges_checked = 0;
                le_lists_bfs_top_down_step(g, frontier, next_frontier, distances, num_frontier_edges, num_edges_checked, deltas, S, iter);
                num_unvisited_edges -= num_edges_checked;
                advance_frontier(frontier, next_frontier);
                frontier_size = frontier->num_vertices;
            }
        } else {
            bool should_switch = ((double) frontier_size) < (((double) g->n) / gamma);
            // should_switch = false;
            if (should_switch) {
                last_step = TOP_DOWN;
                // Reconstruct the current frontier
                reset_frontier(frontier);
                num_unvisited_edges = 0;
                for (int vid = 0; vid < g->n; ++vid) {
                    if (distances[vid] == iter-1) {
                        frontier->vertices[frontier->num_vertices++] = vid;
                    } else if (is_unvisited(vid, distances)) {
                        num_unvisited_edges += g->out_offsets[vid+1] - g->out_offsets[vid];
                    }
                }
                reset_frontier(next_frontier);
                int num_edges_checked = 0;
                le_lists_bfs_top_down_step(g, frontier, next_frontier, distances, num_frontier_edges, num_edges_checked, deltas, S, iter);
                num_unvisited_edges -= num_edges_checked;
                advance_frontier(frontier, next_frontier);
                frontier_size = frontier->num_vertices;
            } else {
                le_lists_bfs_bottom_up_step(g, frontier_size, iter, distances, deltas, S);
            }
        }
        iter++;
    }

    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
}

void le_lists_par(Graph g, std::vector<std::vector<int>> &L_v, std::vector<std::vector<int>> &L_d,
    std::unordered_map<std::string, double> &metrics) {
    auto start_time = std::chrono::steady_clock::now();

    int *deltas = (int *) malloc(g->n*sizeof(int));
    for (int vid = 0; vid < g->n; ++vid) {
        deltas[vid] = std::numeric_limits<int>::max();
    }

    L_v = std::vector<std::vector<int>>(g->n, std::vector<int>()); // v_i's
    L_d = std::vector<std::vector<int>>(g->n, std::vector<int>()); // d(v_i, u)'s

    vertex_set *S = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(S, g->n);
    int *distances = (int *) malloc(g->n*sizeof(int));

    for (int vid = 0; vid < g->n; ++vid) {
        le_lists_bfs_hybrid(g, vid, deltas, S, distances);
        for (int i = 0; i < S->num_vertices; ++i) {            
            int u = S->vertices[i];
            deltas[u] = distances[u];
            L_v[u].push_back(vid);
            L_d[u].push_back(distances[u]);
        }

    }
    auto end_time = std::chrono::steady_clock::now();
    double runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    metrics.insert(std::make_pair("runtime", runtime));

    free(deltas);
    free_vertex_set(S);
    free(distances);
}