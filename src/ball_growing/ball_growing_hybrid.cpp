#include "../util/graph.h"
#include "../bfs/bfs.h"

#include <chrono>
#include <random>
#include <math.h>
#include <omp.h>


void inline bg_top_down_step(Graph &g, vertex_set *frontier, vertex_set *next_frontier, 
    std::unordered_set<int> &unvisited, double *deltas, int *ball_ids, std::vector<int> &owner,
    int *distances, int &num_unvisited, int &num_frontier_edges, int &num_edges_checked, int iter) {
   
    reset_frontier(next_frontier);

    // Check for vertices starting their own BFS
    for (auto &vid : unvisited) {
        if (deltas[vid] < iter) {
            next_frontier->vertices[next_frontier->num_vertices++] = vid;
            owner[vid] = vid;
            ball_ids[vid] = vid;
            distances[vid] = iter;
            num_frontier_edges += g->out_offsets[vid+1] - g->out_offsets[vid];
        }
    }
    // Mark those vertices as unvisited (must do after to not modify set while iterating)
    for (int i = 0; i < next_frontier->num_vertices; ++i) {
        int vid = next_frontier->vertices[i];
        unvisited.erase(vid);
    }

    #pragma omp parallel
    {
        vertex_set *local_frontier = (vertex_set *) malloc(sizeof(vertex_set));
        init_vertex_set(local_frontier, g->n);
        int local_num_frontier_edges = 0;
        int local_num_edges_checked = 0;
        // Iterate over frontier
        #pragma omp for schedule(static)
        for (int i = 0; i < frontier->num_vertices; ++i) {
            int vid = frontier->vertices[i];
            // Iterate over neighbors
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                #pragma omp critical
                {
                    if (unvisited.find(nid) != unvisited.end()) {
                        int ovid = owner[nid];
                        // vid is first to reach nid
                        if (ovid == NOT_OWNED) {
                            owner[nid] = vid;
                            local_frontier->vertices[local_frontier->num_vertices++] = nid;
                            deltas[nid] = deltas[vid];
                            local_num_frontier_edges += g->out_offsets[nid+1] - g->out_offsets[nid];
                        }
                        // nid has been reached by some other vertex in frontier
                        else {
                            double vid_delta_int, ovid_delta_int;
                            double vid_delta_frac = modf(deltas[vid], &vid_delta_int);
                            double ovid_delta_frac = modf(deltas[ovid], &ovid_delta_int);
                            if (vid_delta_frac < ovid_delta_frac) {
                                owner[nid] = vid;
                                deltas[nid] = deltas[vid];
                            }
                        }
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
            num_frontier_edges += local_num_frontier_edges;
            num_edges_checked += local_num_edges_checked;
        }
        free_vertex_set(local_frontier);
    }
    // Mark all vertices in the next frontier as visited
    // Update ball ID of every vertex in the next frontier
    for (int i = 0; i < next_frontier->num_vertices; ++i) {
        int vid = next_frontier->vertices[i];
        unvisited.erase(vid);
        ball_ids[vid] = ball_ids[owner[vid]];
        distances[vid] = iter;
    }
    num_unvisited = unvisited.size();
    advance_frontier(frontier, next_frontier);
}

void inline bg_bottom_up_step(Graph &g, int *distances, double *deltas, 
    int *ball_ids, int &frontier_size, std::vector<int> &owner, int &num_unvisited, int iter) {
    int shared_frontier_size = 0;
    #pragma omp parallel shared(distances, deltas, ball_ids, shared_frontier_size)
    {
        #pragma omp for reduction(+:shared_frontier_size) schedule(static)
        for (int vid = 0; vid < g->n; ++vid) {
            if (is_unvisited(vid, distances)) {
                // Check if vertex should start its own BFS (i.e. its delay is over)
                if (deltas[vid] < iter) {
                    distances[vid] = iter;
                    ball_ids[vid] = vid;
                    shared_frontier_size++;
                } else {
                    // Check if any neighbors are in the frontier
                    double min_delta_frac = 1.0;
                    int oid = -1;
                    double delta_int = 0.0;
                    bool visited = false;
                    for (int eid = g->in_offsets[vid]; eid < g->in_offsets[vid+1]; ++eid) {
                        int nid = g->in_edge_list[eid];
                        if (distances[nid] == iter-1) {
                            visited = true;
                            double delta_frac = modf(deltas[nid], &delta_int);
                            // Check, of the neighbors to reach vid, which is first
                            if (delta_frac < min_delta_frac) {
                                min_delta_frac = delta_frac;
                                oid = nid;
                            }
                            // oid = nid;
                            // ball_ids[vid] = ball_ids[oid];
                            // deltas[vid] = deltas[oid];
                            // distances[vid] = iter;
                            // shared_frontier_size++;
                            // owner[vid] = oid;
                            // break;
                        }
                    }
                    if (visited) {
                        ball_ids[vid] = ball_ids[oid];
                        deltas[vid] = deltas[oid];
                        distances[vid] = iter;
                        shared_frontier_size++;
                        owner[vid] = oid;
                    }
                }
            }
        }
    }
    frontier_size = shared_frontier_size;
    num_unvisited -= frontier_size;
}

void ball_decomp_hybrid(Graph g, double beta, 
    std::vector<std::unordered_set<int>> &collection, std::vector<int> &radii,
    std::unordered_map<std::string, double> &metrics) {
    double alpha = 10.0;
    double gamma = 10.0;

    auto start_time = std::chrono::steady_clock::now();

    // Compute deltas
    double *deltas = (double *) calloc(g->n, sizeof(double));
    compute_deltas(deltas, beta, g->n);

    // Common data structures
    int *distances = (int *) malloc(g->n * sizeof(int));
    memset(distances, UNVISITED, g->n * sizeof(int));
    int *ball_ids = (int *) calloc(g->n, sizeof(int));
    int num_unvisited = g->n;

    // Top-down data structures
    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);
    std::vector<int> owner = std::vector<int>(g->n, NOT_OWNED);
    std::unordered_set<int> unvisited;
    for (int vid = 0; vid < g->n; ++vid) {
        unvisited.insert(vid);
    }

    // Bottom-up data structures
    int frontier_size = 0;
    int iter = 1;
    
    // Hybrid data structures
    int num_frontier_edges = 0; // number of edges to check from frontier
    int num_unvisited_edges = 2*g->m; // number of edges to check from unvisited vertices
    unsigned char last_step = TOP_DOWN;
    // unsigned char last_step = BOTTOM_UP;

    while (num_unvisited > 0) {
        if (last_step == TOP_DOWN) {
            bool should_switch = ((double) num_frontier_edges) > (((double) num_unvisited_edges) / alpha);
            // bool should_switch = false;
            if (should_switch) {
                last_step = BOTTOM_UP;
                bg_bottom_up_step(g, distances, deltas, ball_ids, frontier_size, owner, num_unvisited, iter);
            } else {
                int num_edges_checked = 0;
                bg_top_down_step(g, frontier, next_frontier, unvisited, deltas, ball_ids, owner, distances, num_unvisited, num_frontier_edges, num_edges_checked, iter);
                num_unvisited_edges -= num_edges_checked;
            }
        } else {
            bool should_switch = ((double) frontier_size) < (((double) g->n) / gamma);
            if (should_switch) {
                last_step = TOP_DOWN;
                // Reconstruct the current frontier
                reset_frontier(frontier);
                num_unvisited_edges = 0;
                for (int vid = 0; vid < g->n; ++vid) {
                    if (is_unvisited(vid, distances)) {
                        num_unvisited_edges += g->out_offsets[vid+1] - g->out_offsets[vid];
                    } else {
                        unvisited.erase(vid);
                        if (distances[vid] == iter-1) {
                            frontier->vertices[frontier->num_vertices++] = vid;
                        }
                    }
                }
                int num_edges_checked = 0;
                bg_top_down_step(g, frontier, next_frontier, unvisited, deltas, ball_ids, owner, distances, num_unvisited, num_frontier_edges, num_edges_checked, iter);
                num_unvisited_edges -= num_edges_checked;
            } else {
                bg_bottom_up_step(g, distances, deltas, ball_ids, frontier_size, owner, num_unvisited, iter);
            }
        }
        iter++;
    }

    auto end_time = std::chrono::steady_clock::now();
    std::cout << "Iterations: " << iter << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    std::cout << "Fraction of Intercluster Edges: " << get_frac_intercluster_edges(g, ball_ids) << std::endl;

    double runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    metrics.insert(std::make_pair("runtime", runtime));
    metrics.insert(std::make_pair("iterations", (double) (iter-1)));

    free(deltas);
    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
    free(distances);
    free(ball_ids);
}

