#include "../bfs/bfs.h"
#include <cstdlib>
#include <vector>
#include <unordered_set>

// Constructing the reverse graph (flipping the direction of each edge) of g
Graph reverse_graph(const Graph &g) {
    Graph rev_g = (graph_t *) malloc(sizeof(graph_t));
    rev_g->out_offsets    = (int *) calloc(g->n + 1, sizeof(int));
    rev_g->out_edge_lists = (int *) calloc(2 * g->m, sizeof(int));
    rev_g->in_offsets     = (int *) calloc(g->n + 1, sizeof(int));
    rev_g->in_edge_lists  = (int *) calloc(2 * g->m, sizeof(int));
    rev_g->n = g->n;
    rev_g->m = g->m;
    std::copy(g->out_offsets.begin(), g->out_offsets.end(), rev->g->in_offsets.begin());
    std::copy(g->in_offsets.begin() , g->in_offsets.end() , rev->g->out_offsets.begin());
    std::copy(g->out_edge_lists.begin(), g->out_edge_lists.end(), rev->g->in_edge_lists.begin());
    std::copy(g->in_edge_lists.begin() , g->in_edge_lists.end() , rev->g->out_edge_lists.begin());
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
void construct_frontier_top_down(Graph &g, std::unordered_set<int> &S, vertex_set *frontier, vertex_set *next_frontier, vector<int> &distances) {
    for (int i = 0; i < frontier->num_vertices; ++i) {
        int vid = frontier->vertices[i];
        for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid + 1]; ++eid) {
            int nid = g->out_edge_list[eid];
            if (S.count(nid) && distances[nid] == UNVISITED) {
                distances[nid] = distances[vid] + 1;
                next_frontier->vertices[next_frontier->vertices++] = nid;
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
void construct_frontier_bottom_up(Graph &g, std::unordered_set<int> &S, int &frontier_size, int iter, vector<int> &distances) {
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
    vector<int> distances(g->n, UNVISITED);
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
    vector<int> distances(g->n, UNVISITED);
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
