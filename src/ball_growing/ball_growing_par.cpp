#include "../util/graph.h"
#include "../bfs/bfs.h"

#include <chrono>
#include <random>
#include <math.h>
#include <omp.h>
#include <limits>

#define NOT_OWNED (-1)

/*
 *     ball_ids - maps vertex i to its ball ID
 */
double get_avg_isoperimetric_num(Graph &g, int *ball_ids) {
    // balls - maps owner vertex i to its constituents
    std::vector<std::vector<int>> balls(g->n, std::vector<int>());
    std::unordered_set<int> ball_owners;
    for (int vid = 0; vid < g->n; ++vid) {
        balls[ball_ids[vid]].push_back(vid);
        if (vid == ball_ids[vid]) {
            ball_owners.insert(ball_ids[vid]);
        }
    }
    int num_balls = ball_owners.size();
    std::cout << "Num Balls: " << num_balls << std::endl;
    double total_isoperimetric_num = 0.0;
    
    // Iterate over ball owners; find isoperimetric number of each ball
    for (auto &owner : ball_owners) {
        int degree_sum = 0;
        int boundary_size = 0;
        // Iterate over vertices in ball
        for (auto &vid : balls[owner]) {
            degree_sum += g->out_offsets[vid+1] - g->out_offsets[vid];
            // Iterate over neighbors to count boundary size
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                // Edge spanning this ball and other ball
                if (ball_ids[nid] != ball_ids[vid]) {
                    boundary_size++;
                }
            }
        }
        double isoperimetric_num;
        if (degree_sum == 0) {
            isoperimetric_num = 0;
        } else {
            isoperimetric_num = ((double) boundary_size) / ((double) degree_sum);
        }
        // std::cout << "Isoperimetric Number of Ball " << isoperimetric_num << std::endl;
        total_isoperimetric_num += isoperimetric_num;
    }
    return total_isoperimetric_num / ((double) num_balls);
}

/*
 * Computes S_i = delta_max - delta_i for all i where delta_i ~ Exp(beta).
 *     deltas - output
 *     beta - exponential distribution parameter
 *     n - length of deltas
 */
void inline compute_deltas(double *&deltas, float beta, int n) {
    std::default_random_engine generator;
    // std::random_device rd;
    // std::mt19937 generator(rd());
    std::exponential_distribution<double> distribution(beta);
    double max_delta = -1.0;
    #pragma omp parallel
    {
        #pragma omp for reduction(max:max_delta) schedule(static)
        for (int i = 0; i < n; ++i) {
            deltas[i] = distribution(generator);
            if (deltas[i] > max_delta) {
                max_delta = deltas[i];
            }
        }
        #pragma omp for schedule(static)
        for (int i = 0; i < n; ++i) {
            deltas[i] = max_delta - deltas[i];
        }
    }
    // double min_delta = std::numeric_limits<double>::max();
    // #pragma omp parallel
    // {
    //     #pragma omp for reduction(min:min_delta) schedule(static)
    //     for (int i = 0; i < n; ++i) {
    //         deltas[i] = -distribution(generator);
    //         if (deltas[i] < min_delta) {
    //             min_delta = deltas[i];
    //         }
    //     }
    //     #pragma omp for schedule(static)
    //     for (int i = 0; i < n; ++i) {
    //         deltas[i] = deltas[i] - min_delta;
    //     }
    // }
}

/*
 *
 * Computes a (beta, O(log(m)/beta)) decomposition in expectation in parallel
 * using a bottom-up approach.
 *     g - graph
 *     beta - max isoperimetric number
 *     collection - list of balls
 *     radii - list of radii
 */
void ball_decomp_bottom_up_par(Graph g, float beta, std::vector<std::unordered_set<int>> &collection,
    std::vector<int> &radii) {
    auto start_time = std::chrono::steady_clock::now();

    // Sample deltas from Exp(beta)
    double *deltas = (double *) calloc(g->n, sizeof(double));
    compute_deltas(deltas, beta, g->n);
    // for (int i = 0; i < g->n; i++) {
    //     fprintf(stderr, "vertex %d: %f\n", i, deltas[i]);
    // }

    int iter = 1;
    int frontier_size = 0;
    int *distances = (int *) malloc(g->n * sizeof(int));
    memset(distances, UNVISITED, g->n * sizeof(int)); // careful about UNVISITED being -1
    int *ball_ids = (int *) calloc(g->n, sizeof(int));
    int num_unvisited = g->n;

    // Construct initial frontier
    int initial_frontier_size = 0;
    #pragma omp parallel
    {
        #pragma omp for reduction(+:initial_frontier_size) schedule(static)
        for (int vid = 0; vid < g->n; ++vid) {
            if (deltas[vid] < iter) {
                // Implicitly add to the frontier
                // fprintf(stderr, "starting bfs from vertex %d\n", vid);
                distances[vid] = iter;
                ball_ids[vid] = vid;
                initial_frontier_size++;
            }
        }
    }
    num_unvisited -= initial_frontier_size;

    // BFS
    while (num_unvisited > 0) {
        int shared_frontier_size = 0;
        #pragma omp parallel
        {
            #pragma omp for reduction(+:shared_frontier_size) schedule(static)
            for (int vid = 0; vid < g->n; ++vid) {
                if (is_unvisited(vid, distances)) {
                    // Check if vertex should start its own BFS (i.e. its delay is over)
                    if (deltas[vid] < iter) {
                        // fprintf(stderr, "starting bfs from vertex %d\n", vid);
                        distances[vid] = iter;
                        ball_ids[vid] = vid;
                        shared_frontier_size++;
                    } else {
                        // Check if any neighbors are in the frontier
                        double min_delta_frac = 1.0;
                        int owner = -1;
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
                                    owner = nid;
                                }
                            }
                        }
                        if (visited) {
                            ball_ids[vid] = ball_ids[owner];
                            distances[vid] = iter;
                            shared_frontier_size++;
                        }
                    }
                }
            }
        }
        num_unvisited -= shared_frontier_size;
        iter++;
    }

    // for (int i = 0; i < g->n; ++i) {
    //     fprintf(stderr, "vertex %d is in vertex %d's ball\n", i, ball_ids[i]);
    // }

    auto end_time = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    std::cout << "Avg Isoperimetric Number: " << get_avg_isoperimetric_num(g, ball_ids) << std::endl;

    free(distances);
    free(ball_ids);
}

/*
 * Computes a (beta, O(log(m)/beta)) decomposition in expectation in parallel
 * using a top-down approach.
 *     g - graph
 *     beta - max isoperimetric number
 *     collection - list of balls
 *     radii - list of radii
 */
void ball_decomp_top_down_par(Graph g, float beta, std::vector<std::unordered_set<int>> &collection,
    std::vector<int> &radii) {
    auto start_time = std::chrono::steady_clock::now();

    // Sample deltas from Exp(beta)
    double *deltas = (double *) calloc(g->n, sizeof(double));
    compute_deltas(deltas, beta, g->n);

    // BFS
    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);

    int num_threads = omp_get_num_threads();

    int *ball_ids = (int *) calloc(g->n, sizeof(int));
    int iter = 1;
    std::unordered_set<int> unvisited;
    // Mark all vertices as unvisited
    for (int vid = 0; vid < g->n; ++vid) {
        unvisited.insert(vid);
    }
    std::vector<int> owner = std::vector<int>(g->n, NOT_OWNED);
    while (unvisited.size() > 0) {

        // Add all unvisited vertices with delta < iter into frontier
        for (auto &vid : unvisited) {
            if (deltas[vid] < iter) {
                frontier->vertices[frontier->num_vertices++] = vid;
                ball_ids[vid] = vid;
            }
        }
        // Mark the frontier vertices as visited
        for (int i = 0; i < frontier->num_vertices; ++i) {
            int vid = frontier->vertices[i];
            unvisited.erase(vid);
        }
        // Reset the next frontier
        reset_frontier(next_frontier);
        int thread_id;
        #pragma omp parallel private(thread_id)
        {
            vertex_set *local_frontier = (vertex_set *) malloc(sizeof(vertex_set));
            init_vertex_set(local_frontier, g->n);
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
                            }
                            // nid has been reached by some other vertex in frontier
                            else if (deltas[vid] < deltas[ovid]) {
                                owner[nid] = vid;
                            }
                        }
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
        iter++;
        // Mark all vertices in the next frontier as visited
        // Update ball ID of every vertex in the next frontier
        for (int i = 0; i < next_frontier->num_vertices; ++i) {
            int vid = next_frontier->vertices[i];
            unvisited.erase(vid);
            ball_ids[vid] = ball_ids[owner[vid]];
        }
        advance_frontier(&frontier, &next_frontier);
    }

    auto end_time = std::chrono::steady_clock::now();
    // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    std::cout << "Avg Isoperimetric Number: " << get_avg_isoperimetric_num(g, ball_ids) << std::endl;

    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
    free(deltas);
    free(ball_ids);
}
