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

#endif
