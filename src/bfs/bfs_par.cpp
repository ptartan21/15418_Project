#include "bfs.h"

#include <omp.h>
#include <chrono>

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
    #pragma omp parallel
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
                }
            }
        }
        // Copy local next frontier to real next frontier
        #pragma omp critical
        {
            std::memcpy(next_frontier->vertices+next_frontier->num_vertices, local_frontier->vertices, local_frontier->num_vertices*sizeof(int));
            next_frontier->num_vertices += local_frontier->num_vertices;
        }
        free_vertex_set(local_frontier);
    }
}

/*
 * (PARALLEL)
 *     g - graph
 *     source - starting point for the BFS
 *     distances - output; distances from source
 */
void bfs_top_down_par(Graph g, int source, int *distances, 
    std::unordered_map<std::string, double> &metrics) {
    auto start_time = std::chrono::steady_clock::now();
    // Initialize frontiers
    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);
    memset(distances, UNVISITED, g->n*sizeof(int));
    frontier->vertices[frontier->num_vertices++] = source;
    distances[source] = 0;
    // Main loop
    while (frontier->num_vertices > 0) {
        // Reset the next frontier
        reset_frontier(next_frontier);
        // Construct the next frontier
        construct_frontier_top_down_par(g, frontier, next_frontier, distances);
        // Advance to the next frontier
        advance_frontier(frontier, next_frontier);
    }
    auto end_time = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    // std::cout << "\n";

    double runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    metrics.insert(std::make_pair("runtime", runtime));

    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
}

/*
 * (PARALLEL)
 * Constructs the next frontier using a bottom-up approach.
 *     g - graph
 *     frontier_size - number of vertices in the frontier
 *     iter - current iteration; distance of vertices added to frontier from source
 *     distances - distances from source
 */
 void inline construct_frontier_bottom_up_par(Graph g, int &frontier_size,
    int iter, int *distances) {
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

/*
 * (PARALLEL)
 * Bottom-up BFS.
 *     g - graph
 *     source - starting point for the BFS
 *     distances - output; distances from source
 */
void bfs_bottom_up_par(Graph g, int source, int *distances,
    std::unordered_map<std::string, double> &metrics) {
    auto start_time = std::chrono::steady_clock::now();
    int iter = 1;
    memset(distances, UNVISITED, g->n*sizeof(int));
    int frontier_size = 1;
    distances[source] = 0;
    while (frontier_size > 0) {
        construct_frontier_bottom_up_par(g, frontier_size, iter, distances);
        iter++;
    }
    auto end_time = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    // std::cout << "\n";

    double runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    metrics.insert(std::make_pair("runtime", runtime));
}
