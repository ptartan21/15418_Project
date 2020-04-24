#include "../util/graph.h"

#include <chrono>
#include <vector>
#include <unordered_set>

#define PRESENT 1
#define NOT_PRESENT 0

double get_avg_isoperimetric_num(Graph &g, std::vector<std::unordered_set<int>> &collection) {
    double total_isoperimetric_num = 0.0;
    int num_balls = collection.size();
    for (auto &ball : collection) {
        int degree_sum = 0;
        int boundary_size = 0;
        for (auto &vid : ball) {
            degree_sum += g->out_offsets[vid+1] - g->out_offsets[vid];
            // Iterate over neighbors to count boundary size
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                // Edge spanning this ball and other ball
                if (ball.find(nid) == ball.end()) {
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
 * Grows a ball from the given source until its isoperimetric number is at 
 * most beta.
 *     g - graph
 *     source - source vertex to start growing from
 *     beta - max isoperimetric number
 *     present - indicates if each vertex is considered present in the graph still
 *     ball - collection of vertices representing the ball
 *     R - radius of the final ball
 */
void grow_ball(Graph g, int source, float beta, std::vector<bool> &present, std::unordered_set<int> &ball, int &R) {
    // Initialize metrics
    int r = 1;
    float isoperimetric_num = (float) g->m; // max possible
    int degree_sum = g->out_offsets[source+1] - g->out_offsets[source];
    // Initialize boundary data structures
    std::unordered_set<int> boundary_vertices;
    std::unordered_set<int> next_boundary_vertices;
    boundary_vertices.insert(source);
    // Initialize ball
    ball.clear();
    ball.insert(source);
    // Iterate until satisfying isoperimetric number
    while (isoperimetric_num > beta) {
        int next_boundary_size = 0;
        // Iterate over boundary vertices
        for (const auto &vid : boundary_vertices) {
            // Iterate over neighbors
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                // Ignore if the vertex has already been viisted (no longer present in graph)
                if (present[nid] == NOT_PRESENT) {
                    continue;
                }
                // With respect to vertices
                if (next_boundary_vertices.find(nid) == next_boundary_vertices.end()) {
                    next_boundary_vertices.insert(nid);
                    degree_sum += g->out_offsets[nid+1] - g->out_offsets[nid];
                }
                // With respect to edges
                next_boundary_size++;
            }
        }
        // Compute new isoperimetric number
        if (degree_sum > 0) {
            isoperimetric_num = ((float) next_boundary_size) / ((float) degree_sum);
        } else {
            // Clamp the isoperimetric number in case of division-by-zero
            isoperimetric_num = 0;
        }
        // Mark newly visited vertices as no longer present in the graph and
        // add them to the current ball
        for (const auto &vid : next_boundary_vertices) {
            present[vid] = NOT_PRESENT;
            ball.insert(vid);
        }
        // Advance the boundary
        boundary_vertices.swap(next_boundary_vertices);
        next_boundary_vertices.clear();
        r++;
    }
    R = r-1;
}


/*
 * Computes a (beta, O(log(m)/beta)) decomposition.
 *     g - graph
 *     beta - max isoperimetric number
 *     collection - list of balls
 *     radii - list of radii
 */
void ball_decomp_seq(Graph g, float beta, std::vector<std::unordered_set<int>> &collection, std::vector<int> &radii) {
    auto start_time = std::chrono::steady_clock::now();

    // Initially, all vertices present
    std::vector<bool> present(g->n, PRESENT);
    std::unordered_set<int> ball;
    // Iterate over all vertices (in-order)
    for (int vid = 0; vid < g->n; ++vid) {
        int R;
        if (present[vid] == PRESENT) {
            // Grow ball from vertex
            grow_ball(g, vid, beta, present, ball, R);
            // Store ball and radius
            collection.push_back(ball);
            radii.push_back(R);
            ball.clear();
        }
    }

    auto end_time = std::chrono::steady_clock::now();
    // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() << " ns" << std::endl;
    std::cout << "Num Balls: " << collection.size() << std::endl;
    std::cout << "Avg Isoperimetric Number: " << get_avg_isoperimetric_num(g, collection) << std::endl;
}
