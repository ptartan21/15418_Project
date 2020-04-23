#include <algorithm>
#include <random>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "util/graph.h"
#include "bfs/bfs_seq.cpp"
#include "bfs/bfs_par.cpp"
#include "ball_growing/ball_growing.cpp"

/*
 * Populate g
 */
void inline load_graph(std::string source, Graph g) {

}

int main(int argc, char **argv) {
    int n, m;
    Graph g = (graph_t *) malloc(sizeof(graph_t));
    std::string graph_in(argv[1]);
    std::string line;
    std::ifstream source;
    source.open(graph_in);
    // retrieving the number of vertices and the number of edges
    std::getline(source, line);
    std::istringstream iss(line);
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

    int *distances = (int *) calloc(n, sizeof(int));
    bfs_top_down_seq(g, 0, distances);
    for (int i = 0; i < n; ++i) {
        std::cout << distances[i] << " ";
    }
    std::cout << "\n\n";

    for (int i = 0; i < n; ++i) {
        distances[i] = 0;
    }
    bfs_bottom_up_par(g, 0, distances);
    for (int i = 0; i < n; ++i) {
        std::cout << distances[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Sequential Ball Growing" << std::endl;
    vector<unordered_set<int>> collection;
    vector<int> radii;
    ball_decomp_seq(g, 0.75, collection, radii);
}
