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
void set_u(unsigned char *S1, unsigned char *S2, unsigned char *res, int n) {
    for (int i = 0; i < n; ++i) {
        res[i] = S1[i] | S2[i];
    }
}

// computes the set intersection of S1 and S2
void set_i(unsigned char *S1, unsigned char *S2, unsigned char *res, int n) {
    for (int i = 0; i < n; ++i) {
        res[i] = S1[i] & S2[i];
    }
}

// computes the set difference S1 - S2
void set_d(unsigned char *S1, unsigned char *S2, unsigned char *res, int n) {
    for (int i = 0; i < n; ++i) {
        if (!S1[i] || (S1[i] && S2[i])) {
            res[i] = 0;
        } else {
            res[i] = 1;
        }
    }
}

// // computes the set union of S1 and S2
// void set_u_par(unsigned char *S1, unsigned char *S2, unsigned char *res, int n) {
//     for (int i = 0; i < n; ++i) {
//         res[i] = S1[i] | S2[i];
//     }
// }

// // computes the set intersection of S1 and S2
// void set_i_par(unsigned char *S1, unsigned char *S2, unsigned char *res) {
//     for (int i = 0; i < n; ++i) {
//         res[i] = S1[i] & S2[i];
//     }
// }

// // computes the set difference S1 - S2
// void set_d_par(unsigned char *S1, unsigned char *S2, unsigned char *res) {
//     for (int i = 0; i < n; ++i) {
//         if (!S1[i] || (S1[i] && S2[i])) {
//             res[i] = 0;
//         } else {
//             res[i] = 1;
//         }
//     }
// }

// Constructing the reverse graph (flipping the direction of each edge) of g
Graph reverse_graph(Graph g) {
    Graph rev_g = alloc_graph(g->n, g->m);
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
void construct_frontier_top_down(Graph &g, unsigned char *S, unsigned char *reach, vertex_set *frontier, vertex_set *next_frontier, std::vector<int> &distances) {
    for (int i = 0; i < frontier->num_vertices; ++i) {
        int vid = frontier->vertices[i];
        for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid + 1]; ++eid) {
            int nid = g->out_edge_list[eid];
            if (S[nid] && distances[nid] == UNVISITED) {
                reach[nid] = 1;
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
void construct_frontier_bottom_up(Graph &g, unsigned char *S, unsigned char *reach, int &frontier_size, int iter, std::vector<int> &distances) {
    frontier_size = 0;
    for (int vid = 0; vid < g->n; ++vid) {
        if (S[vid] && distances[vid] == UNVISITED) {
            for (int eid = g->in_offsets[vid]; eid <g->in_offsets[vid + 1]; ++eid) {
                int nid = g->in_edge_list[eid];
                if (distances[nid] == iter - 1) {
                    reach[vid] = 1;
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
void bfs_top_down_seq(Graph &g, unsigned char *S, unsigned char *reach, int source) {
    vertex_set *frontier      = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);
    std::vector<int> distances(g->n, UNVISITED);
    frontier->vertices[frontier->num_vertices++] = source;
    distances[source] = 0;
    while (frontier->num_vertices > 0) {
        reset_frontier(next_frontier);
        construct_frontier_top_down(g, S, reach, frontier, next_frontier, distances);
        advance_frontier(frontier, next_frontier);
    }
    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
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
void bfs_bottom_up_seq(Graph &g, unsigned char *S, unsigned char *reach, int source) {
    int iter          = 1;
    int frontier_size = 1;
    std::vector<int> distances(g->n, UNVISITED);
    distances[source] = 0;
    while (frontier_size > 0) {
        construct_frontier_bottom_up(g, S, reach, frontier_size, iter, distances);
        ++iter;
    }
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
            std::memcpy(next_frontier->vertices + next_frontier->num_vertices, local_frontier->vertices, local_frontier->num_vertices * sizeof(int));
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

/***** HYBRID *****/

void inline bfs_top_down_step(Graph &g, std::unordered_set<int> &S, vertex_set *frontier, vertex_set *next_frontier, std::vector<int> &distances, int &num_frontier_edges, int &num_edges_checked) {
    int local_num_frontier_edges = 0;
    int local_num_edges_checked = 0;
    #pragma omp parallel private(local_num_frontier_edges, local_num_edges_checked)
    {
        vertex_set *local_frontier = (vertex_set *) malloc(sizeof(vertex_set));
        init_vertex_set(local_frontier, g->n);
        #pragma omp for schedule(static)
        for (int i = 0; i < frontier->num_vertices; ++i) {
            int vid = frontier->vertices[i];
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                if (S.count(nid) && __sync_bool_compare_and_swap(&distances[nid], UNVISITED, distances[vid] + 1)) {
                    local_frontier->vertices[local_frontier->num_vertices++] = nid;
                    local_num_frontier_edges += g->out_offsets[nid + 1] - g->out_offsets[nid];
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

void inline bfs_bottom_up_step(Graph &g, std::unordered_set<int> &S, int &frontier_size, int iter, std::vector<int> &distances) {
    int shared_frontier_size = 0;
    #pragma omp parallel
    {
        // Iterate over vertices
        #pragma omp for reduction(+:shared_frontier_size) schedule(static)
        for (int vid = 0; vid < g->n; ++vid) {
            if (S.count(vid) && distances[vid] == UNVISITED) {
                for (int eid = g->in_offsets[vid]; eid < g->in_offsets[vid + 1]; ++eid) {
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

std::unordered_set<int> bfs_hybrid(Graph g, std::unordered_set<int> &S, int source) {
    double alpha = 10.0;
    double beta  = 10.0;

    std::vector<int> distances(g->n, UNVISITED);

    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);

    int frontier_size = 0;
    int iter          = 1;

    int num_frontier_edges  = 0;
    int num_unvisited_edges = 2 * g->m;
    unsigned char last_step = TOP_DOWN;

    frontier->vertices[frontier->num_vertices++] = source;
    frontier_size++;
    distances[source] = 0;
    num_frontier_edges += g->out_edge_list[source + 1] - g->out_edge_list[source];

    while (frontier_size > 0) {
        if (last_step == TOP_DOWN) {
            bool should_switch = ((double) num_frontier_edges) > (((double) num_unvisited_edges) / alpha);
            if (should_switch) {
                last_step = BOTTOM_UP;
                bfs_bottom_up_step(g, S, frontier_size, iter, distances);
            } else {
                reset_frontier(next_frontier);
                int num_edges_checked = 0;
                bfs_top_down_step(g, S, frontier, next_frontier, distances, num_frontier_edges, num_edges_checked);
                num_unvisited_edges -= num_edges_checked;
                advance_frontier(frontier, next_frontier);
            }
        } else {
            bool should_switch = ((double) frontier_size) < (((double) g->n) / beta);
            if (should_switch) {
                last_step = TOP_DOWN;
                reset_frontier(frontier);
                num_unvisited_edges = 0;
                for (int vid = 0; vid < g->n; ++vid) {
                    if (distances[vid] == iter - 1) {
                        frontier->vertices[frontier->num_vertices++] = vid;
                    } else if (distances[vid] == UNVISITED) {
                        num_unvisited_edges += g->out_offsets[vid + 1] - g->out_offsets[vid];
                    }
                }
                reset_frontier(next_frontier);
                int num_edges_checked = 0;
                bfs_top_down_step(g, S, frontier, next_frontier, distances, num_frontier_edges, num_edges_checked);
                num_unvisited_edges -= num_edges_checked;
                advance_frontier(frontier, next_frontier);
            } else {
                bfs_bottom_up_step(g, S, frontier_size, iter, distances);
            }
        }
        /*
        std::cout << "Iteration: " << iter << std::endl;
        if (last_step == TOP_DOWN) {
            std::cout << "Top Down" << std::endl;
        } else {
            std::cout << "Bottom Up" << std::endl;
        }
        */
        iter++;
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
