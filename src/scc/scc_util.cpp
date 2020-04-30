#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_set>
#include <omp.h>

#include "../util/graph.h"
#include "../bfs/bfs.h"

// Copying the elements from src to dest
void deep_copy(int *src, int *dest, int n) {
    for (int i = 0; i < n; ++i) {
        dest[i] = src[i];
    }
}

// computes the set union of S1 and S2
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
    std::unordered_set<int> intersect;
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

// Constructing the reverse graph (flipping the direction of each edge) of g
Graph reverse_graph(Graph g) {
    Graph rev_g = (graph_t *) malloc(sizeof(graph_t));
    rev_g->out_offsets   = (int *) calloc(g->n + 1, sizeof(int));
    rev_g->out_edge_list = (int *) calloc(2 * g->m, sizeof(int));
    rev_g->in_offsets    = (int *) calloc(g->n + 1, sizeof(int));
    rev_g->in_edge_list  = (int *) calloc(2 * g->m, sizeof(int));
    rev_g->n = g->n;
    rev_g->m = g->m;
    deep_copy(g->out_offsets, rev_g->in_offsets, g->n + 1);
    deep_copy(g->in_offsets, rev_g->out_offsets, g->n + 1);
    deep_copy(g->out_edge_list, rev_g->in_edge_list, 2 * g->m);
    deep_copy(g->in_edge_list, rev_g->out_edge_list, 2 * g->m);
    return rev_g;
}

/***** SEQUENTIAL *****/
/*
 * (SEQUENTIAL)
 * Constructs the next frontier using a top-down approach
 *     g             - graph
 *     S             - subgraph
 *     frontier      - current frontier
 *     next_frontier - next frontier
 *     distances     - distance from source
 */
void construct_frontier_top_down(Graph &g, std::unordered_set<int> &S, vertex_set *frontier, vertex_set *next_frontier, std::vector<int> &distances) {
    for (int i = 0; i < frontier->num_vertices; ++i) {
        int vid = frontier->vertices[i];
        for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid + 1]; ++eid) {
            int nid = g->out_edge_list[eid];
            if (S.count(nid) && distances[nid] == UNVISITED) {
                distances[nid] = distances[vid] + 1;
                next_frontier->vertices[next_frontier->num_vertices++] = nid;
            }
        }
    }
}

/*
 * (SEQUENTIAL)
 * Constructs the next frontier using a bottom-up approach
 *     g             - graph
 *     S             - subgraph
 *     frontier_size - number of vertices in the current frontier
 *     iter          - distance of vertices added to frontier from source
 *     distances     - distance from source
 */
void construct_frontier_bottom_up(Graph &g, std::unordered_set<int> &S, int &frontier_size, int iter, std::vector<int> &distances) {
    frontier_size = 0;
    for (int vid = 0; vid < g->n; ++vid) {
        if (S.count(vid) && distances[vid] == UNVISITED) {
            for (int eid = g->in_offsets[vid]; eid <g->in_offsets[vid + 1]; ++eid) {
                int nid = g->in_edge_list[eid];
                if (distances[nid] == iter - 1) {
                    distances[vid] = iter;
                    frontier_size += 1;
                    break;
                }
            }
        }
    }
}

/*
 * (SEQUENTIAL)
 * Top-down BFS for SCC
 *     g      - graph
 *     S      - subgraph
 *     source - starting point for BFS
 * Output
 *     reach  - set of vertices visited
 */
std::unordered_set<int> bfs_top_down_seq(Graph &g, std::unordered_set<int> &S, int source) {
    vertex_set *frontier      = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);
    std::vector<int> distances(g->n, UNVISITED);
    frontier->vertices[frontier->num_vertices++] = source;
    distances[source] = 0;
    while (frontier->num_vertices > 0) {
        reset_frontier(next_frontier);
        construct_frontier_top_down(g, S, frontier, next_frontier, distances);
        advance_frontier(frontier, next_frontier);
    }
    std::unordered_set<int> reach;
    for (int vid = 0; vid < g->n; ++vid) {
        if (distances[vid] != UNVISITED) {
            reach.insert(vid);
        }
    }
    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
    return reach;
}

/*
 * (SEQUENTIAL)
 * Bottom-up BFS for SCC
 *     g      - graph
 *     S      - subgraph
 *     source - starting point for BFS
 * Output
 *     reach  - set of vertices visited
 */
std::unordered_set<int> bfs_bottom_up_seq(Graph &g, std::unordered_set<int> &S, int source) {
    int iter          = 1;
    int frontier_size = 1;
    std::vector<int> distances(g->n, UNVISITED);
    distances[source] = 0;
    while (frontier_size > 0) {
        construct_frontier_bottom_up(g, S, frontier_size, iter, distances);
        ++iter;
    }
    std::unordered_set<int> reach;
    for (int vid = 0; vid < g->n; ++vid) {
        if (distances[vid] != UNVISITED) {
            reach.insert(vid);
        }
    }
    return reach;
}

/***** PARALLEL *****/

/*
 * (PARALLEL)
 * Constructs the next frontier using a top-down approach
 *     g             - graph
 *     S             - subgraph
 *     frontier      - current frontier
 *     next_frontier - next frontier
 *     distances     - distance from source
 */
void construct_frontier_top_down_par(Graph &g, std::unordered_set<int> &S, vertex_set *frontier, vertex_set *next_frontier, std::vector<int> &distances) {
    #pragma omp parallel
    {
        vertex_set *local_frontier = (vertex_set *) malloc(sizeof(vertex_set));
        init_vertex_set(local_frontier, g->n);
        #pragma omp for schedule(static)
        for (int i = 0; i < frontier->num_vertices; ++i) {
            int vid = frontier->vertices[i];
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid + 1]; ++eid) {
                int nid = g->out_edge_list[eid];
                if (S.count(nid) && __sync_bool_compare_and_swap(&distances[nid], UNVISITED, distances[vid] + 1)) {
                    local_frontier->vertices[local_frontier->num_vertices++] = nid;
                }
            }
        }
        #pragma omp critical
        {
            std:memcpy(next_frontier->vertices + next_frontier->num_vertices, local_frontier->vertices, local_frontier->num_vertices * sizeof(int));
            next_frontier->num_vertices += local_frontier->num_vertices;
        }
        free_vertex_set(local_frontier);
    }
}

/*
 * (PARALLEL)
 * Constructs the next frontier using a bottom-up approach
 *     g             - graph
 *     S             - subgraph
 *     frontier_size - number of vertices in the current frontier
 *     iter          - distance of vertices added to frontier from source
 *     distances     - distance from source
 */
void construct_frontier_bottom_up_par(Graph &g, std::unordered_set<int> &S, int &frontier_size, int iter, std::vector<int> &distances) {
    int shared_frontier_size = 0;
    #pragma omp parallel
    {
        #pragma omp for reduction(+:shared_frontier_size) schedule(static)
        for (int vid = 0; vid < g->n; ++vid) {
            if (S.count(vid) && distances[vid] == UNVISITED) {
                for (int eid = g->in_offsets[vid]; eid <g->in_offsets[vid + 1]; ++eid) {
                    int nid = g->in_edge_list[eid];
                    if (distances[nid] == iter - 1) {
                        distances[vid] = iter;
                        shared_frontier_size += 1;
                        break;
                    }
                }
            }
        }
    }
    frontier_size = shared_frontier_size;
}

/*
 * (PARALLEL)
 * Top-down BFS for SCC
 *     g      - graph
 *     S      - subgraph
 *     source - starting point for BFS
 * Output
 *     reach  - set of vertices visited
 */
std::unordered_set<int> bfs_top_down_par(Graph &g, std::unordered_set<int> &S, int source) {
    vertex_set *frontier      = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);
    std::vector<int> distances(g->n, UNVISITED);
    frontier->vertices[frontier->num_vertices++] = source;
    distances[source] = 0;
    while (frontier->num_vertices > 0) {
        reset_frontier(next_frontier);
        construct_frontier_top_down_par(g, S, frontier, next_frontier, distances);
        advance_frontier(frontier, next_frontier);
    }
    std::unordered_set<int> reach;
    for (int vid = 0; vid < g->n; ++vid) {
        if (distances[vid] != UNVISITED) {
            reach.insert(vid);
        }
    }
    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
    return reach;
}

/*
 * (PARALLEL)
 * Bottom-up BFS for SCC
 *     g      - graph
 *     S      - subgraph
 *     source - starting point for BFS
 * Output
 *     reach  - set of vertices visited
 */
std::unordered_set<int> bfs_bottom_up_par(Graph &g, std::unordered_set<int> &S, int source) {
    int iter          = 1;
    int frontier_size = 1;
    std::vector<int> distances(g->n, UNVISITED);
    distances[source] = 0;
    while (frontier_size > 0) {
        construct_frontier_bottom_up_par(g, S, frontier_size, iter, distances);
        ++iter;
    }
    std::unordered_set<int> reach;
    for (int vid = 0; vid < g->n; ++vid) {
        if (distances[vid] != UNVISITED) {
            reach.insert(vid);
        }
    }
    return reach;
}


