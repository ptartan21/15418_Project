#ifndef __GRAPH_H__
#define __GRAPH_H__

typedef struct {
    int n; // number of vertices
    int m; // number of edges

    int *vertex_offsets;
    // Consider vertex i. Let offset_i = vertex_offsets[i], and let 
    // offset_{i+1} = vertex_offsets[i+1]
    // edge_list[offset_i], edge_list[offset_i+1], ..., edge_list[offset_{i+1}-1]
    // are the out neighbors of vertex i
    int *edge_list;

} graph_t;

#endif
