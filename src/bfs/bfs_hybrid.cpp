#include "bfs.h"

#include <omp.h>
#include <chrono>

void inline top_down_step() {
    int local_num_vertices;
    #pragma omp parallel private(local_num_vertices)
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

void inline bottom_up_step() {
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
    *frontier_size = shared_frontier_size;
}

void hybrid_bfs(Graph g, int source, int *distances) {
    // Reset distances
    memset(distances, UNVISITED, g->n*sizeof(int));

    // Top-down data structures
    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);

    // Bottom-up data structures
    int frontier_size = 0; // invariant: frontier_size = frontier->num_vertices

    // Insert the source vertex into the initial frontier
    frontier->vertices[frontier->num_vertices++] = source;
    frontier_size++;
    distances[source] = 0;

    while (frontier_size > 0) {

    }
}