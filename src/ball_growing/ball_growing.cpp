#include "../util/graph.h"

#include <vector>
#include <unordered_set>

#define PRESENT 1
#define NOT_PRESENT 0

void grow_ball(Graph g, int source, float beta, std::vector<bool> &present, std::unordered_set<int> &ball, int &R) {
    int r = 1;
    float isoperimetric_num = (float) g->m; // max possible
    std::unordered_set<int> boundary_vertices;
    std::unordered_set<int> next_boundary_vertices;
    boundary_vertices.insert(source);
    int degree_sum = g->out_offsets[source+1] - g->out_offsets[source];
    ball.clear();
    ball.insert(source);
    while (isoperimetric_num > beta) {
        int next_boundary_size = 0;
        for (const auto &vid : boundary_vertices) {
            for (int eid = g->out_offsets[vid]; eid < g->out_offsets[vid+1]; ++eid) {
                int nid = g->out_edge_list[eid];
                if (!present[nid]) {
                    continue;
                }
                // Unvisited
                if (next_boundary_vertices.find(nid) != next_boundary_vertices.end()) {
                    next_boundary_vertices.insert(nid);
                    degree_sum += g->out_offsets[nid+1] - g->out_offsets[nid];
                }
                next_boundary_size++;
            }
        }
        // Update present
        // Recompute isoperimetric_num by dividing
        isoperimetric_num = next_boundary_size / degree_sum;
        for (const auto &vid : next_boundary_vertices) {
            present[vid] = NOT_PRESENT;
            ball.insert(vid);
        }
        boundary_vertices.swap(next_boundary_vertices);
        next_boundary_vertices.clear();
        r++;
    }
    R = r;
}



void ball_decomp_seq(Graph g, float beta, std::vector<std::unordered_set<int>> &collection, std::vector<int> &radii) {
    //
    std::vector<bool> present(g->n, PRESENT);
    std::unordered_set<int> ball;
    for (int vid = 0; vid < g->n; ++vid) {
        int R;
        if (present[vid] == PRESENT) {
          grow_ball(g, vid, beta, present, ball, R);
          collection.push_back(ball);
          radii.push_back(R);
          ball.clear();
        }
    }
}
