#include "../util/graph.h"

#include <chrono>
#include <vector>
#include <unordered_set>

/*
 * Returns the average isoperimetric number.
 *     g - graph
 *     collection - collection of balls
 */
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
        total_isoperimetric_num += isoperimetric_num;
    }
    return total_isoperimetric_num / ((double) num_balls);
}

/*
 * Returns the fraction of edges that are intercluster.
 *     g - graph
 *     collection - collection of balls
 */
double get_frac_intercluster_edges(Graph &g, std::vector<std::unordered_set<int>> &collection) {
    int total_intercluster_edges = 0;
    for (auto &ball : collection) {
        for (auto &vid : ball) {
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                if (ball.find(nid) == ball.end()) {
                    total_intercluster_edges++;
                }
            }
        }
    }
    // Divide by 2 because each edge is double-counted
    return ((double) total_intercluster_edges) / ((double) g->m) / 2.0;
}

/*
 * Returns the degree of the given vertex according to which vertices are
 * still present.
 *     g - graph
 *     vid - vertex ID
 *     present - indicates if each vertex remains in the graph
 */
int inline find_degree(Graph g, int vid, std::vector<bool> &present) {
    int degree = 0;
    for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
        int nid = g->out_edge_list[eid];
        if (present[nid]) {
            degree++;
        }
    }
    return degree;
}

/*
 * Grows a ball from the source until the isoperimetric number is at most beta.
 *     g - graph
 *     source - source vertex to start BFS
 *     beta - desired isoperimetric number
 *     present - indicates if each vertex remains in the graph
 *     ball - output; set of vertices in the ball
 *     R - radius of the outputted ball
 */
void grow_ball(Graph g, int source, float beta, std::vector<bool> &present, std::unordered_set<int> &ball, int &R) {
    int r = 1;
    float isoperimetric_num;
    ball.insert(source);
    std::unordered_set<int> frontier;
    std::unordered_set<int> next_frontier;
    frontier.insert(source);
    int source_degree = find_degree(g, source, present);
    if (source_degree > 0) {
        isoperimetric_num = 1;
    } else {
        isoperimetric_num = 0;
    }
    int boundary_size = source_degree;
    while (isoperimetric_num > beta) {
        for (auto &vid : frontier) {
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                if (present[nid]) {
                    if (ball.find(nid) == ball.end()) {
                        // (vid, nid) - ball to non-ball edge; nid not in next frontier
                        if (next_frontier.find(nid) == next_frontier.end()) {
                            boundary_size--;
                            next_frontier.insert(nid);
                        }
                        // (vid, nid) - ball to non-ball edge; nid already in next frontier
                        else {
                            boundary_size--;
                        }
                    }
                    // (vid, nid) - ball to ball edge
                    else {
                        // Do nothing
                    }
                }
            }
        }
        for (auto &vid : next_frontier) {
            ball.insert(vid);
        }
        int degree_sum = 0;
        for (auto &vid : next_frontier) {
            int degree = 0;
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                if (present[nid]) {
                    degree++;
                    if (ball.find(nid) == ball.end() && next_frontier.find(nid) == next_frontier.end()) {
                        boundary_size++;
                    }
                }
            }
            degree_sum += degree;
        }
        if (degree_sum == 0) {
            isoperimetric_num = 0;
        } else {
            isoperimetric_num = ((float) boundary_size) / ((float) degree_sum);
        }
        frontier.swap(next_frontier);
        next_frontier.clear();
        r++;
    }
    // Remove the vertices in the ball from the graph
    for (auto &vid : ball) {
        present[vid] = false;
    }
    R = r-1;
}

/*
 * Computes a low diameter decomposition of the graph by repeatedly growing a
 * ball and removing it from the graph. 
 *     g - graph
 *     beta - desired isoperimetric number per ball
 *     collection - output; collection of balls
 *     radii - output; radii of the balls
 */
void ball_decomp_seq(Graph g, float beta, std::vector<std::unordered_set<int>> &collection, std::vector<int> &radii, std::unordered_map<std::string, double> &metrics) {
    auto start_time = std::chrono::steady_clock::now();
    std::vector<bool> present(g->n, true);
    std::unordered_set<int> ball;
    for (int vid = 0; vid < g->n; ++vid) {
        if (present[vid]) {
            int R;
            grow_ball(g, vid, beta, present, ball, R);
            collection.push_back(ball);
            radii.push_back(R);
            ball.clear();
        }
    }
    auto end_time = std::chrono::steady_clock::now();
    // std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    // std::cout << "Num Balls: " << collection.size() << std::endl;
    // std::cout << "Fraction of Intercluster Edges: " << get_frac_intercluster_edges(g, collection) << std::endl;

    double runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    metrics.insert(std::make_pair("runtime", runtime));
}
