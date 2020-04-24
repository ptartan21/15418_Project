#include <algorithm>
#include <random>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "util/graph.h"
#include "bfs/bfs_seq.cpp"
#include "bfs/bfs_par.cpp"
#include "ball_growing/ball_growing_seq.cpp"
#include "ball_growing/ball_growing_par.cpp"

/*
 * Populate g
 */
void inline load_graph(std::string graph_in, Graph &g) {
    int n, m;
    std::ifstream source;
    source.open(graph_in);
    std::string line;
    std::getline(source, line);
    std::istringstream iss(line);
    // retrieving the number of vertices and the number of edges
    iss >> n >> m;
    std::vector<std::vector<int>> in_mapper(n + 1, std::vector<int>{});
    // populating the out-neighbor portion of the graph
    g->out_offsets   = (int *) calloc(n + 1, sizeof(int));
    g->out_edge_list = (int *) calloc(m,     sizeof(int));
    g->n = n;
    g->m = m;
    int v   = -1;
    int off = 0;
    int pos = 0;
    int cur = -1;
    while (std::getline(source, line)) {
        std::istringstream iss(line);
        iss >> cur;
        g->out_offsets[off++] = pos;
        while (iss >> v) {
            g->out_edge_list[pos++] = v;
            in_mapper[v].push_back(cur);
        }
    }
    g->out_offsets[n] = m;

    // populating the in-neighbor portion of the graph
    g->in_offsets   = (int *) calloc(n + 1, sizeof(int));
    g->in_edge_list = (int *) calloc(m,     sizeof(int));
    off = 0;
    pos = 0;
    // for each vertex->in_neighbor pair
    for (int vtx = 0; vtx < in_mapper.size(); ++vtx) {
        g->in_offsets[off++] = pos;
        // for each in-neighbor of the current vertex
        for (auto &in_nbor : in_mapper[vtx]) {
            g->in_edge_list[pos++] = in_nbor;
        }
    }
    std::cout << std::endl;
    source.close();
    std::cout << "Succesfully Loaded Graph" << std::endl;

}

void inline bfs_top_down_seq_wrapper(Graph &g) {
    int n = g->n;
    int *distances = (int *) calloc(n, sizeof(int));

    free(distances);
}

void inline bfs_top_down_par_wrapper(Graph &g) {
    int n = g->n;
    int *distances = (int *) calloc(n, sizeof(int));
    bfs_top_down_par(g, 0, distances);
    for (int i = 0; i < g->n; ++i) {
        std::cout << distances[i] << " ";
    }
    std::cout << std::endl;
    free(distances);
}

void inline bfs_bottom_up_seq_wrapper(Graph &g) {
    int n = g->n;
    int *distances = (int *) calloc(n, sizeof(int));

    free(distances);
}

void inline bfs_bottom_up_par_wrapper(Graph &g) {
    int n = g->n;
    int *distances = (int *) calloc(n, sizeof(int));

    free(distances);
}

void inline ball_decomp_seq_wrapper(Graph g, float beta) {
    std::cout << "Ball Decomposition with beta = " << std::setprecision(2) << beta << " (Sequential)" << std::endl;
    std::vector<std::unordered_set<int>> collection;
    std::vector<int> radii;
    ball_decomp_seq(g, beta, collection, radii);
    for (int i = 0; i < radii.size(); ++i) {
        for (auto &vid : collection[i]) {
            std::cout << vid << " ";
        }
        std::cout << "\nRadius: " << radii[i] << "\n\n";
    }
}

void inline ball_decomp_par_wrapper(Graph g, float beta) {
    std::cout << "Ball Decomposition with beta = " << std::setprecision(2) << beta << " (Parallel)" << std::endl;
    std::vector<std::unordered_set<int>> collection;
    std::vector<int> radii;
    ball_decomp_par(g, beta, collection, radii);
}

int main(int argc, char **argv) {
    std::string graph_in(argv[1]);
    Graph g = (graph_t *) malloc(sizeof(graph_t));
    load_graph(graph_in, g);

    // ball_decomp_seq_wrapper(g, 0.25);
    // bfs_top_down_par_wrapper(g);
    ball_decomp_par_wrapper(g, 0.25);

    free(g);
}
