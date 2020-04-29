#include "bfs.h"

#include <omp.h>
#include <chrono>

void inline bfs_top_down_step(Graph g, vertex_set *frontier, vertex_set *next_frontier, 
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

void inline bfs_bottom_up_step(Graph g, int &frontier_size, int iter, int *distances) {
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

// Assume undirected graph
void bfs_hybrid(Graph g, int source, int *distances,
    std::unordered_map<std::string, double> &metrics) {
    double alpha = 10.0;
    double beta = 10.0;
    auto start_time = std::chrono::steady_clock::now();

    // Reset distances
    memset(distances, UNVISITED, g->n*sizeof(int));

    // Top-down data structures
    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);

    // Bottom-up data structures
    int frontier_size = 0; // invariant: frontier_size = frontier->num_vertices
    int iter = 1;

    // Hybrid data structures
    int num_frontier_edges = 0; // number of edges to check from frontier
    int num_unvisited_edges = 2*g->m; // number of edges to check from unvisited vertices
    unsigned char last_step = TOP_DOWN;

    // Insert the source vertex into the initial frontier
    frontier->vertices[frontier->num_vertices++] = source;
    frontier_size++;
    distances[source] = 0;
    num_frontier_edges += g->out_edge_list[source+1] - g->out_edge_list[source];

    while (frontier_size > 0) {
        if (last_step == TOP_DOWN) {
            bool should_switch = ((double) num_frontier_edges) > (((double) num_unvisited_edges) / alpha);
            if (should_switch) {
                last_step = BOTTOM_UP;
                bfs_bottom_up_step(g, frontier_size, iter, distances);
            } else {
                reset_frontier(next_frontier);
                int num_edges_checked = 0;
                bfs_top_down_step(g, frontier, next_frontier, distances, num_frontier_edges, num_edges_checked);
                num_unvisited_edges -= num_edges_checked;
                advance_frontier(frontier, next_frontier);
            }
        } else {
            bool should_switch = ((double) frontier_size) < (((double) g->n) / beta);
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
                bfs_top_down_step(g, frontier, next_frontier, distances, num_frontier_edges, num_edges_checked);
                num_unvisited_edges -= num_edges_checked;
                advance_frontier(frontier, next_frontier);
            } else {
                bfs_bottom_up_step(g, frontier_size, iter, distances);
            }
        }
        std::string step;
        if (last_step == TOP_DOWN) {
            step = "top down";
        } else {
            step = "bottom up";
        }
        std::cout << "Iteration " << iter << ": " << step << std::endl;
        iter++;
    }
    auto end_time = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    std::cout << "\n";

    double runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    metrics.insert(std::make_pair("runtime", runtime));

    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
}