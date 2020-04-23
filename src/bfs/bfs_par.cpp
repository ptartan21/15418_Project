#include "../util/graph.h"
#include "bfs.h"

#include <cstring>
#include <stdbool.h>
#include <stdlib.h>
#include <omp.h>

#define UNVISITED (-1)

/*
 * Initializes the vertex set.
 */
void inline init_vertex_set(vertex_set *v_set, int max_vertices) {
    v_set->num_vertices = 0;
    v_set->vertices = (int *) calloc(max_vertices, sizeof(int));
}

/*
 * Resets the frontier by setting its number of vertices to zero.
 */
void inline reset_frontier(vertex_set *frontier) {
    frontier->num_vertices = 0;
}

/*
 * Swap the frontier pointers.
 */
void inline advance_frontier(vertex_set **frontier_ptr, vertex_set **next_frontier_ptr) {
    vertex_set *temp = *frontier_ptr;
    *frontier_ptr = *next_frontier_ptr;
    *next_frontier_ptr = temp;
}

bool inline is_unvisited(int vid, int *distances) {
    return distances[vid] == UNVISITED;
}

void inline mark_unvisited(int vid, int*distances) {
    distances[vid] = UNVISITED;
}

/*
 * (PARALLEL)
 * Constructs the next frontier using a top-down approach.
 *     g - graph
 *     frontier - current frontier
 *     next_frontier - next frontier
 *     distances - distances from source
 */
void inline construct_frontier_top_down_par(Graph g, vertex_set *frontier,
    vertex_set *next_frontier, int *distances) {
    // vid - vertex ID
    // eid - edge ID
    // nid - neighbor ID
    int local_num_vertices;
    #pragma omp parallel private(local_num_vertices)
    {
        local_num_vertices = 0;
        int num_threads = omp_get_num_threads();
        int *local_frontier = (int *) calloc(g->n/num_threads+1, sizeof(int));
        // Iterate over frontier
        #pragma omp for schedule(static)
        for (int i = 0; i < frontier->num_vertices; ++i) {
            int vid = frontier->vertices[i];
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                // Use CAS to guarantee nid is only added to the frontier once
                if (__sync_bool_compare_and_swap(&distances[nid], UNVISITED, distances[vid]+1)) {
                    local_frontier[local_num_vertices++] = nid;
                }
            }
        }
        // Copy local next frontier to real next frontier
        #pragma omp critical
        {
            std::memcpy(next_frontier->vertices+next_frontier->num_vertices, local_frontier, local_num_vertices*sizeof(int));
            next_frontier->num_vertices += local_num_vertices;
        }
        free(local_frontier);
    }
}

/*
 * (PARALLEL)
 *     g - graph
 *     source - starting point for the BFS
 *     distances - output; distances from source
 */
void bfs_top_down_par(Graph g, int source, int *distances) {
    // Initialize frontiers
    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);
    #pragma omp parallel for schedule(static)
    for (int vid = 0; vid < g->n; ++vid) {
        mark_unvisited(vid, distances);
    }
    frontier->vertices[frontier->num_vertices++] = source;
    distances[source] = 0;
    // Main loop
    while (frontier->num_vertices > 0) {
        // Reset the next frontier
        reset_frontier(next_frontier);
        // Construct the next frontier
        construct_frontier_top_down_par(g, frontier, next_frontier, distances);
        // Advance to the next frontier
        advance_frontier(&frontier, &next_frontier);
    }
    free(frontier);
    free(next_frontier);
}
