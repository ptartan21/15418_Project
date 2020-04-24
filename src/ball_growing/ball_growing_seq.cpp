#include "../util/graph.h"

#include <chrono>
#include <vector>
#include <unordered_set>

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
    for (auto &vid : ball) {
        present[vid] = false;
        // std::cout << vid << " ";
    }
    // std::cout << "\n";
    R = r-1;
}

void ball_decomp_seq(Graph g, float beta, std::vector<std::unordered_set<int>> &collection, std::vector<int> &radii) {
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
            // std::cout << "Avg Isoperimetric Number: " << get_avg_isoperimetric_num(g, collection) << std::endl;
        }
    }
    auto end_time = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << " ms" << std::endl;
    std::cout << "Num Balls: " << collection.size() << std::endl;
    std::cout << "Avg Isoperimetric Number: " << get_avg_isoperimetric_num(g, collection) << std::endl;
}
