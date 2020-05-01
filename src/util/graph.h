#ifndef __GRAPH_H__
#define __GRAPH_H__

typedef struct {
    int n; // number of vertices
    int m; // number of edges

    // Add dummy element (at the end) mapping n to m
    // This allows for offsets[nid] to offsets[nid+1] indexing without out-of-bounds
    int *out_offsets;
    int *out_edge_list;

    int *in_offsets;
    int *in_edge_list;

} graph_t, *Graph;

struct vertex_set {
    int num_vertices;
    int *vertices;
};

/*
 * Allocates a graph data structure
 */
Graph alloc_graph(int n, int m) {
    Graph g = (graph_t *) malloc(sizeof(graph_t));
    g->n = n;
    g->m = m;
    g->out_offsets   = (int *) calloc(n + 1, sizeof(int));
    g->out_edge_list = (int *) calloc(2 * m, sizeof(int));
    g->in_offsets    = (int *) calloc(n + 1, sizeof(int));
    g->in_edge_list  = (int *) calloc(2 * m, sizeof(int));
    return g;
}

/*
 * Initializes the vertex set.
 */
void inline init_vertex_set(vertex_set *v_set, int max_vertices) {
    v_set->num_vertices = 0;
    v_set->vertices = (int *) calloc(max_vertices, sizeof(int));
}

void inline free_vertex_set(vertex_set *v_set) {
    free(v_set->vertices);
    free(v_set);
}

#endif
