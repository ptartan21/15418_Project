/*
 * Resets the frontier by setting its number of vertices to zero.
 */
void inline reset_frontier(vertex_set *frontier) {
    frontier->num_vertices = 0;
}

/*
 * Swap the frontier pointers.
 */
void inline advance_frontier(vertex_set *frontier, vertex_set *next_frontier) {
    vertex_set temp = *frontier;
    *frontier = *next_frontier;
    *next_frontier = temp;
}


bool inline is_unvisited(int vid, int *distances) {
    return distances[vid] == UNVISITED;
}

void inline mark_unvisited(int vid, int*distances) {
    distances[vid] = UNVISITED;
}

void print_array(int *arr, int n) {
    for (int i = 0; i < n; ++i) {
        std::cout << arr[i] << " ";
    }
    std::cout << "\n";
}