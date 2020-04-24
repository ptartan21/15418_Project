#include "../util/graph.h"
#include "../bfs/bfs.h"

#include <random>
#include <omp.h>

#define NOT_OWNED (-1)

/*
 * Computes a (beta, O(log(m)/beta)) decomposition in expectation in parallel.
 *     g - graph
 *     beta - max isoperimetric number
 *     collection - list of balls
 *     radii - list of radii
 */
void ball_decomp_par(Graph g, float beta, std::vector<std::unordered_set<int>> &collection,
    std::vector<int> &radii) {
    double *deltas = (double *) malloc(g->n * sizeof(double));
    std::default_random_engine generator;
    std::exponential_distribution<double> distribution(beta);
    double max_delta = -1.0;
    for (int i = 0; i < g->n; ++i) {
        deltas[i] = distribution(generator);
        if (deltas[i] > max_delta) {
            max_delta = deltas[i];
        }
    }
    // Make every delta positive
    for (int i = 0; i < g->n; ++i) {
        deltas[i] = max_delta - deltas[i];
    }
    // BFS
    vertex_set *frontier = (vertex_set *) malloc(sizeof(vertex_set));
    vertex_set *next_frontier = (vertex_set *) malloc(sizeof(vertex_set));
    init_vertex_set(frontier, g->n);
    init_vertex_set(next_frontier, g->n);

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
        #pragma omp parallel
        {
            int num_threads = omp_get_num_threads();
            vertex_set *local_frontier = (vertex_set *) malloc(sizeof(vertex_set));
            init_vertex_set(local_frontier, g->n/num_threads+1);
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
            iter++;
        }
        // Mark all vertices in the next frontier as visited
        // Update ball ID of every vertex in the next frontier
        for (int i = 0; i < next_frontier->num_vertices; ++i) {
            int vid = next_frontier->vertices[i];
            unvisited.erase(vid);
            ball_ids[vid] = ball_ids[owner[vid]];
        }
        advance_frontier(&frontier, &next_frontier);

        fprintf(stderr, "Unvisited:\n");
        for (auto &vid : unvisited) {
            fprintf(stderr, "%d ", vid);
        }
        fprintf(stderr, "\n\n");
    }

    std::unordered_set<int> ball;
    for (int i = 0; i < g->n; ++i) {
        for (int v = 0; v < g->n; ++v) {
            if (ball_ids[v] == i) {
                ball.insert(v);
            }
        }
        collection.push_back(ball);
        ball.clear();
    }

    free_vertex_set(frontier);
    free_vertex_set(next_frontier);
    free(deltas);
    free(ball_ids);
}
