#include "../util/graph.h"
#include "../bfs/bfs.h"

#include <omp.h>
#include <unordered_set>

/*
 * Computes S_i = delta_max - delta_i for all i where delta_i ~ Exp(beta).
 *     deltas - output
 *     beta - exponential distribution parameter
 *     n - length of deltas
 */
void inline compute_deltas(double *&deltas, double beta, int n) {
    std::default_random_engine generator(418);
    std::exponential_distribution<double> distribution(beta);
    double max_delta = -1.0;
    // omp_set_num_threads(1);
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
    // omp_set_num_threads(8);
}

/*
 * Returns the fraction of edges that are intercluster.
 *     g - graph
 *     ball_ids - maps vertex i to its ball ID
 */
double get_frac_intercluster_edges(Graph &g, int *ball_ids) {
    // balls - maps owner vertex i to its constituents
    std::vector<std::vector<int>> balls(g->n, std::vector<int>());
    std::unordered_set<int> ball_owners;
    for (int vid = 0; vid < g->n; ++vid) {
        balls[ball_ids[vid]].push_back(vid);
        if (vid == ball_ids[vid]) {
            ball_owners.insert(ball_ids[vid]);
        }
    }
    // int num_balls = ball_owners.size();
    // std::cout << "Num Balls: " << num_balls << std::endl;
    int total_intercluster_edges = 0;
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
        total_intercluster_edges += boundary_size;
    }
    // Divide by 2 because each edge is double-counted
    return ((double) total_intercluster_edges) / ((double) g->m) / 2.0;
}