#include "bfs.h"

/*
 * Initializes the vertex set.
 */
void inline init_vertex_set(vertex_set *v_set, int max_vertices) {
    v_set->num_vertices = 0;
    v_set->vertices = (int *) calloc(max_vertices, sizeof(int));
}

/*
 * Resets the frontier by setting its number of vertices to zero.
 */
void inline reset_frontier(vertex_set *frontier) {
    frontier->num_vertices = 0;
}

/*
 * Swap the frontier pointers.
 */
void inline advance_frontier(vertex_set **frontier_ptr, vertex_set **next_frontier_ptr) {
    vertex_set *temp = *frontier_ptr;
    *frontier_ptr = *next_frontier_ptr;
    *next_frontier_ptr = temp;
}

bool inline is_unvisited(int vid, int *distances) {
    return distances[vid] == UNVISITED;
}

void inline mark_unvisited(int vid, int*distances) {
    distances[vid] = UNVISITED;
}
