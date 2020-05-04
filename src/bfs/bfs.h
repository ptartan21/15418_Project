#ifndef __BFS_H__
#define __BFS_H__

#include "../util/graph.h"

#include <cstring>
#include <stdbool.h>
#include <stdlib.h>
#include <unordered_map>
#include <string>

#define UNVISITED (-1)
#define TOP_DOWN 1
#define BOTTOM_UP 0
#define NOT_OWNED (-1)

#include "bfs_util.cpp"
#include "../ball_growing/ball_growing_util.cpp"
#include "../scc/scc_util.cpp"

#endif
