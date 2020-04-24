#include "bfs.h"

#include <chrono>

/*
 * (SEQUENTIAL)
 * Constructs the next frontier using a top-down approach.
 *     g - graph
 *     frontier - current frontier
 *     next_frontier - next frontier
 *     distances - distances from source
 */
void inline construct_frontier_top_down_seq(Graph g, vertex_set *frontier,
    vertex_set *next_frontier, int *distances) {
    // vid - vertex ID
    // eid - edge ID
    // nid - neighbor ID

    // Iterate over frontier
    for (int i = 0; i < frontier->num_vertices; ++i) {
        int vid = frontier->vertices[i];
        // Iterate over out-neighbors
        for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
            int nid = g->out_edge_list[eid];
            if (is_unvisited(nid, distances)) {
                // Store distance
                distances[nid] = distances[vid]+1;
                // Add the neighbor to the next frontier
                next_frontier->vertices[next_frontier->num_vertices++] = nid;
            }
        }
    }
}

/*
 * (SEQUENTIAL)
 * Top-down BFS.
 *     g - graph
 *     source - starting point for the BFS
 *     distances - output; distances from source
 */
void bfs_top_down_seq(Graph g, int source, int *distances) {
    auto start_time = std::chrono::steady_clock::now();
    // Initialize frontiers
    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);
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
        construct_frontier_top_down_seq(g, frontier, next_frontier, distances);
        // Advance to the next frontier
        advance_frontier(&frontier, &next_frontier);
    }
    auto end_time = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
}

/*
 * (SEQUENTIAL)
 * Constructs the next frontier using a bottom-up approach.
 *     g - graph
 *     frontier_size - number of vertices in the frontier
 *     iter - current iteration; distance of vertices added to frontier from source
 *     distances - distances from source
 */
void inline construct_frontier_bottom_up_seq(Graph g, int *frontier_size,
    int iter, int *distances) {
    // Do not store an explicit frontier; use distances instead
    // Only store the size of the frontier
    *frontier_size = 0;
    // Iterate over vertices
    for (int vid = 0; vid < g->n; ++vid) {
        // Add vertex to implicit frontier if any in-neighbors were visited in
        // the previous iteration (i.e. are in the current frontier)
        if (is_unvisited(vid, distances)) {
            for (int eid = g->in_offsets[vid]; eid < g->in_offsets[vid+1]; ++eid) {
                int nid = g->in_edge_list[eid];
                if (distances[nid] == iter-1) {
                    distances[vid] = iter;
                    *frontier_size = *frontier_size + 1;
                    break;
                }
            }
        }
    }
}

/*
 * (SEQUENTIAL)
 * Bottom-up BFS.
 *     g - graph
 *     source - starting point for the BFS
 *     distances - output; distances from source
 */
void bfs_bottom_up_seq(Graph g, int source, int *distances) {
    auto start_time = std::chrono::steady_clock::now();
    int iter = 1;
    int *frontier_size = (int *) calloc(1, sizeof(int));
    for (int vid = 0; vid < g->n; ++vid) {
        mark_unvisited(vid, distances);
    }
    *frontier_size = 1;
    distances[source] = 0;
    while (*frontier_size > 0) {
        construct_frontier_bottom_up_seq(g, frontier_size, iter, distances);
        iter++;
    }
    auto end_time = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    free(frontier_size);
}
