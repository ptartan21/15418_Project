#include <algorithm>
#include <random>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>

#include "util/graph.h"
#include "bfs/bfs_seq.cpp"
#include "bfs/bfs_par.cpp"
#include "bfs/bfs_hybrid.cpp"
#include "ball_growing/ball_growing_seq.cpp"
#include "ball_growing/ball_growing_par.cpp"
#include "ball_growing/ball_growing_hybrid.cpp"
#include "scc/scc_seq.cpp"
#include "scc/scc_par.cpp"
#include "scc/scc_hybrid.cpp"
#include "le_lists/le_lists_seq.cpp"
#include "le_lists/le_lists_par.cpp"

/*
 * Populate g from the given input file.
 *     graph_in - input filename
 *     g - graph
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
    g->out_edge_list = (int *) calloc(2*m,   sizeof(int));
    g->n = n;
    g->m = m;
    int v   = -1;
    int off = 0;
    int pos = 0;
    int cur = -1;
    int num_edges_stored = 0;
    while (std::getline(source, line)) {
        std::istringstream iss(line);
        iss >> cur;
        g->out_offsets[off++] = pos;
        while (iss >> v) {
            g->out_edge_list[pos++] = v;
            in_mapper[v].push_back(cur);
            num_edges_stored++;
        }
    }
    g->out_offsets[n] = num_edges_stored;

    // populating the in-neighbor portion of the graph
    g->in_offsets   = (int *) calloc(n + 1, sizeof(int));
    g->in_edge_list = (int *) calloc(2*m,   sizeof(int));
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
    g->in_offsets[n] = num_edges_stored;
    source.close();
    std::cout << "Succesfully Loaded Graph" << std::endl;

}

void inline bfs_top_down_seq_wrapper(Graph g, std::string out_filename) {
    std::cout << "Top Down BFS (Sequential)" << std::endl;
    int *distances = (int *) calloc(g->n, sizeof(int));
    std::unordered_map<std::string, double> metrics;
    double runtime = std::numeric_limits<double>::max();
    // Run 3 times and take min runtime
    for (int i = 0; i < 3; ++i) {
        bfs_top_down_seq(g, 0, distances, metrics);
        #pragma omp barrier
        runtime = std::min(runtime, metrics.find("runtime")->second);
        metrics.clear();
    }

    // Write metrics to outfile
    std::ofstream outfile;
    outfile.open(out_filename, std::ios_base::app);
    outfile << std::to_string(runtime) << " ";

    free(distances);
}

void inline bfs_top_down_par_wrapper(Graph g, std::string out_filename) {
    std::cout << "Top Down BFS (Parallel)" << std::endl;
    int *distances = (int *) calloc(g->n, sizeof(int));
    std::unordered_map<std::string, double> metrics;
    double runtime = std::numeric_limits<double>::max();
    // Run 3 times and take min runtime
    for (int i = 0; i < 3; ++i) {
        bfs_top_down_par(g, 0, distances, metrics);
        #pragma omp barrier
        runtime = std::min(runtime, metrics.find("runtime")->second);
        metrics.clear();
    }

    // Write metrics to outfile
    std::ofstream outfile;
    outfile.open(out_filename, std::ios_base::app);
    outfile << std::to_string(runtime) << " ";

    free(distances);
}

void inline bfs_bottom_up_seq_wrapper(Graph g, std::string out_filename) {
    std::cout << "Bottom Up BFS (Sequential)" << std::endl;
    int *distances = (int *) calloc(g->n, sizeof(int));
    std::unordered_map<std::string, double> metrics;
    double runtime = std::numeric_limits<double>::max();
    // Run 3 times and take min runtime
    for (int i = 0; i < 3; ++i) {
        bfs_bottom_up_seq(g, 0, distances, metrics);
        #pragma omp barrier
        runtime = std::min(runtime, metrics.find("runtime")->second);
        metrics.clear();
    }

    // Write metrics to outfile
    std::ofstream outfile;
    outfile.open(out_filename, std::ios_base::app);
    outfile << std::to_string(runtime) << " ";

    free(distances);
}

void inline bfs_bottom_up_par_wrapper(Graph g, std::string out_filename) {
    std::cout << "Bottom Up BFS (Parallel)" << std::endl;
    int *distances = (int *) calloc(g->n, sizeof(int));
    std::unordered_map<std::string, double> metrics;
    double runtime = std::numeric_limits<double>::max();
    // Run 3 times and take min runtime
    for (int i = 0; i < 3; ++i) {
        bfs_bottom_up_par(g, 0, distances, metrics);
        #pragma omp barrier
        runtime = std::min(runtime, metrics.find("runtime")->second);
        metrics.clear();
    }

    // Write metrics to outfile
    std::ofstream outfile;
    outfile.open(out_filename, std::ios_base::app);
    outfile << std::to_string(runtime) << " ";

    free(distances);
}

void inline bfs_hybrid_wrapper(Graph g, std::string out_filename) {
    std::cout << "Hybrid BFS (Parallel)" << std::endl;
    int *distances = (int *) calloc(g->n, sizeof(int));
    std::unordered_map<std::string, double> metrics;
    double runtime = std::numeric_limits<double>::max();
    // Run 3 times and take min runtime
    for (int i = 0; i < 3; ++i) {
        bfs_hybrid(g, 0, distances, metrics);
        #pragma omp barrier
        runtime = std::min(runtime, metrics.find("runtime")->second);
        metrics.clear();
    }

    // Write metrics to outfile
    std::ofstream outfile;
    outfile.open(out_filename, std::ios_base::app);
    outfile << std::to_string(runtime) << "\n";

    free(distances);
}

void inline ball_decomp_seq_wrapper(Graph g, float beta) {
    std::cout << "Ball Decomposition with beta = " << std::setprecision(2) << beta << " (Sequential)" << std::endl;
    std::vector<std::unordered_set<int>> collection;
    std::vector<int> radii;
    ball_decomp_seq(g, beta, collection, radii);
}

void inline ball_decomp_bottom_up_par_wrapper(Graph g, float beta, std::string out_filename) {
    std::cout << "Ball Decomposition with beta = " << std::setprecision(2) << beta << " (Bottom-Up, Parallel)" << std::endl;
    std::vector<std::unordered_set<int>> collection;
    std::vector<int> radii;
    std::unordered_map<std::string, double> metrics;
    ball_decomp_bottom_up_par(g, beta, collection, radii, metrics);

    // Write metrics to outfile
    double runtime = metrics.find("runtime")->second;
    std::ofstream outfile;
    outfile.open(out_filename, std::ios_base::app);
    outfile << std::to_string(runtime) << " ";
}

void inline ball_decomp_top_down_par_wrapper(Graph g, float beta, std::string out_filename) {
    std::cout << "Ball Decomposition with beta = " << std::setprecision(2) << beta << " (Top-Down, Parallel)" << std::endl;
    std::vector<std::unordered_set<int>> collection;
    std::vector<int> radii;
    std::unordered_map<std::string, double> metrics;
    ball_decomp_top_down_par(g, beta, collection, radii, metrics);

    // Write metrics to outfile
    double runtime = metrics.find("runtime")->second;
    std::ofstream outfile;
    outfile.open(out_filename, std::ios_base::app);
    outfile << std::to_string(runtime) << " ";
}

void inline ball_decomp_hybrid_wrapper(Graph g, float beta, std::string out_filename) {
    std::cout << "Ball Decomposition with beta = " << std::setprecision(2) << beta << " (Hybrid, Parallel)" << std::endl;
    std::vector<std::unordered_set<int>> collection;
    std::vector<int> radii;
    std::unordered_map<std::string, double> metrics;
    ball_decomp_hybrid(g, beta, collection, radii, metrics);

    // Write metrics to outfile
    double runtime = metrics.find("runtime")->second;
    std::ofstream outfile;
    outfile.open(out_filename, std::ios_base::app);
    outfile << std::to_string(runtime) << "\n";
}

void inline bfs_correctness_wrapper(Graph g) {
    int n = g->n;
    int *distances_ref = (int *) calloc(n, sizeof(int));
    int *distances_test = (int *) calloc(n, sizeof(int));
    std::unordered_map<std::string, double> metrics_ref;
    std::unordered_map<std::string, double> metrics_test;
    bfs_top_down_seq(g, 0, distances_ref, metrics_ref);
    bfs_hybrid(g, 0, distances_test, metrics_test);
    for (int i = 0; i < n; ++i) {
        if (distances_ref[i] != distances_test[i]) {
            fprintf(stderr, "Mismatch at %d\n", i);
        }
    }
    free(distances_ref);
    free(distances_test);
}

// 0 == bottom up seq, 1 == top down seq
void inline scc_seq_wrapper(Graph &g, int method) {
    std::vector<std::unordered_set<int>> all_scc;
    if (method == 0) {
        std::cout << "SCC Bottom Up BFS (Sequential)" << std::endl;
        compute_scc_seq(all_scc, g, 0);
    } else if (method == 1) {
        std::cout << "SCC Top Down BFS (Sequential)" << std::endl;
        compute_scc_seq(all_scc, g, 1);
    }
}

// 0 == bottom up seq, 1 == top down seq
void inline scc_par_wrapper(Graph &g, int method) {
    std::vector<std::unordered_set<int>> all_scc;
    if (method == 0) {
        std::cout << "SCC Bottom Up BFS (Parallel)" << std::endl;
        compute_scc_par(all_scc, g, 0);
    } else if (method == 1) {
        std::cout << "SCC Top Down BFS (Parallel)" << std::endl;
        compute_scc_par(all_scc, g, 1);
    }
}

void inline scc_hybrid_wrapper(Graph &g) {
    std::vector<std::unordered_set<int>> all_scc;
    std::cout << "SCC Hybrid BFS" << std::endl;
    compute_scc_hybrid(all_scc, g);
}

void inline le_lists_seq_wrapper(Graph g) {
    std::cout << "LE-Lists (Seq)" << std::endl;
    std::vector<std::vector<int>> L_v;
    std::vector<std::vector<int>> L_d;
    std::unordered_map<std::string, double> metrics;
    double runtime = std::numeric_limits<double>::max();
    for (int i = 0; i < 3; ++i) {
        le_lists_seq(g, L_v, L_d, metrics);
        #pragma omp barrier
        runtime = std::min(runtime, metrics.find("runtime")->second);
        metrics.clear();
    }



    // for (int vid = 0; vid < g->n; ++vid) {
    //     std::cout << "LE-List for vertex " << vid << std::endl;
    //     std::vector<int> L_vid_v = L_v[vid];
    //     std::vector<int> L_vid_d = L_d[vid];
    //     for (int j = 0; j < L_vid_v.size(); ++j) {
    //         int nid = L_vid_v[j];
    //         std::cout << "vertex: " << nid << ", distance: " << L_vid_d[j] << "    ";
    //     }
    //     std::cout << "\n";
    // }

    std::cout << std::to_string(runtime) << std::endl;

}

void inline le_lists_par_wrapper(Graph g) {
    std::cout << "LE-Lists (Par)" << std::endl;
    std::vector<std::vector<int>> L_v;
    std::vector<std::vector<int>> L_d;
    std::unordered_map<std::string, double> metrics;
    double runtime = std::numeric_limits<double>::max();
    for (int i = 0; i < 3; ++i) {
        le_lists_par(g, L_v, L_d, metrics);
        #pragma omp barrier
        runtime = std::min(runtime, metrics.find("runtime")->second);
        metrics.clear();
    }



    // for (int vid = 0; vid < g->n; ++vid) {
    //     std::cout << "LE-List for vertex " << vid << std::endl;
    //     std::vector<int> L_vid_v = L_v[vid];
    //     std::vector<int> L_vid_d = L_d[vid];
    //     for (int j = 0; j < L_vid_v.size(); ++j) {
    //         int nid = L_vid_v[j];
    //         std::cout << "vertex: " << nid << ", distance: " << L_vid_d[j] << "    ";
    //     }
    //     std::cout << "\n";
    // }

    std::cout << std::to_string(runtime) << std::endl;

}

int main(int argc, char **argv) {
    std::string graph_in(argv[1]);
    Graph g = (graph_t *) malloc(sizeof(graph_t));
    load_graph(graph_in, g);

    // int num_threads = 8;
    // omp_set_num_threads(num_threads);
    // std::cout << "Number of Threads: " << num_threads << std::endl;

    // ball_decomp_seq_wrapper(g, 0.25);
    ball_decomp_top_down_par_wrapper(g, 0.5, "results/ball_growing/bg.txt");
    ball_decomp_bottom_up_par_wrapper(g, 0.5, "results/ball_growing/bg.txt");
    ball_decomp_hybrid_wrapper(g, 0.5, "results/ball_growing/bg.txt");

    // for (int num_threads = 1; num_threads <= 8; ++num_threads) {
    //     std::string num_threads_str = std::to_string(num_threads);
    //     omp_set_num_threads(num_threads);
    //     std::cout << "Number of Threads: " << num_threads << std::endl;
    //     std::string out_filename = "results/bfs/bfs_powerlaw2_" + num_threads_str + ".txt";
    //     bfs_top_down_seq_wrapper(g, out_filename);
    //     #pragma omp barrier
    //     bfs_top_down_par_wrapper(g, out_filename);
    //     #pragma omp barrier
    //     bfs_bottom_up_seq_wrapper(g, out_filename);
    //     #pragma omp barrier
    //     bfs_bottom_up_par_wrapper(g, out_filename);
    //     #pragma omp barrier
    //     bfs_hybrid_wrapper(g, out_filename);
    //     #pragma omp barrier
    // }

    // bfs_correctness_wrapper(g);

    // scc_seq_wrapper(g, 0);
    // scc_seq_wrapper(g, 1);
    // scc_par_wrapper(g, 0);
    // scc_par_wrapper(g, 1);
    // scc_hybrid_wrapper(g);
    // le_lists_seq_wrapper(g);
    // le_lists_par_wrapper(g);

    free(g);

    return 0;
}
