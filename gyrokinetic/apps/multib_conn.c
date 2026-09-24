#include <gkyl_multib_conn.h>
#include <gkyl_gyrokinetic_multib.h>
#include <gkyl_gyrokinetic_multib_priv.h>
#include <gkyl_multib_comm_conn.h>

/**
 * Count number of distinct elements in an array of ints
 * .
 * @param a input array of ints
 * @param n length of a
 * return number of unique elements in a
 */
static
int count_distinct(int a[], int n)
{
   int i, j, count = 1;
   for (i = 1; i < n; i++) { // Check if a[i] is a new element
     for (j = 0; j < i; j++) {
       if (a[i] == a[j])    // Check if a[i] has already been found 
          break;            // Break if it is a duplicate
     }
     if (i == j)
       count++;     //increment the number of distinct elements
   }
   return count;
}

/**
 * Populate an output array with the distinct elements in an array of ints
 * .
 * @param a input array of ints
 * @param n length of input array
 * @param unique_array on output contains the unique elements in a
 * return number of unique elements in a
 */
static int
get_unique(int *a, int n, int *unique_array) {
   unique_array[0] = a[0]; // The first element of a is the first unique element
   int i, j, count = 1;
   for (i = 1; i < n; i++) { // Check if a[i] is a new element
     for (j = 0; j < i; j++) {
       if (a[i] == a[j])    // Check if a[i] has already been found 
          break;            // Break if it is a duplicate
     }
     if (i == j) {
       count++;     //increment the number of distinct elements
       unique_array[i] = a[i];
     }
   }
   return count;
}


/** Insert an element at the beginning of an array of ints
 * @param arr of length n+1 (padded with one dummy value at the end)
 * @param n number of values in arr before insertion
 * @param new val value to insert
*/
static void
insert_below(int* arr, int n, int new_val)
{
  int temp_arr[GKYL_MAX_BLOCKS] = {-1};
  for (int i = 0; i<n; i++) {
    temp_arr[i+1] = arr[i];
  }
  temp_arr[0] = new_val;
  n+=1;
  for (int i = 0; i<n; i++) {
    arr[i] = temp_arr[i];
  }

}

/** Insert an element at the end of an array of ints
 * @param arr of length n+1 (padded with one dummy value at the end)
 * @param n number of values in arr before insertion
 * @param new val value to insert
*/
static void
insert_above(int* arr, int n, int new_val)
{
  arr[n] = new_val;
}

/**
 *  Get the block indices of neighbors (adjacent blocks) in a direction.
 * 
 * @param block_topo block topology object
 * @param bidx block index
 * @param dir
 * @param neighbor_idxs on output indices of neighboring blocks
 * return number of neighbors
 */
int
get_neighbors(struct gkyl_block_topo *block_topo, int bidx, int dir, int *neighbor_idxs)
{
  struct gkyl_block_connections conn = block_topo->conn[bidx];
  int neighbor_num = 0;
  for (int e = 0; e < 2; e++) {
    if (conn.connections[dir][e].edge != GKYL_PHYSICAL) {
      neighbor_idxs[neighbor_num] = conn.connections[dir][e].bid;
      neighbor_num += 1;
    }
  }
  return neighbor_num;
}

/**
 *  Get the number of neighbors (adjacent blocks) in a direction.
 * 
 * @param block_topo block topology object
 * @param bidx block index
 * @param direction
 * return number of neighbors
 */
int
get_num_neighbors(struct gkyl_block_topo *block_topo, int bidx, int dir)
{
  int neighbor_idxs[1000] = {-1};
  int neighbor_num = get_neighbors(block_topo, bidx, dir, neighbor_idxs);
  return neighbor_num;
}

/**
 *  Get the block index of adjacent block below in a direction.
 * 
 * @param block_topo block topology object
 * @param bidx block index
 * @param dir
 * @param neighbor_idxs on output index of neighboring block below
 * return number of neighbors below (0 or 1)
*/
int
get_below(struct gkyl_block_topo *block_topo, int bidx, int dir, int *neighbor_idxs)
{
  struct gkyl_block_connections conn = block_topo->conn[bidx];
  int neighbor_num = 0;
  int e = 0;
  if (conn.connections[dir][e].edge != GKYL_PHYSICAL) {
    neighbor_idxs[neighbor_num] = conn.connections[dir][e].bid;
    neighbor_num += 1;
  }
  return neighbor_num;
}

/**
 *  Get the number of adjacent blocks below in a direction.
 * 
 * @param block_topo block topology object
 * @param bidx block index
 * @param direction
 * return number of neighbors
 */
int
get_num_below(struct gkyl_block_topo *block_topo, int bidx, int dir)
{
  int neighbor_idxs[1000] = {-1};
  int neighbor_num = get_neighbors(block_topo, bidx, dir, neighbor_idxs);
  return neighbor_num;
}

/**
 *  Get the block index of adjacent block above in a direction.
 * 
 * @param block_topo block topology object
 * @param bidx block index
 * @param dir
 * @param neighbor_idxs on output index of neighboring block above
 * return number of neighbors above (0 or 1)
*/
int
get_above(struct gkyl_block_topo *block_topo, int bidx, int dir, int *neighbor_idxs)
{
  struct gkyl_block_connections conn = block_topo->conn[bidx];
  int neighbor_num = 0;
  int e = 1;
  if (conn.connections[dir][e].edge != GKYL_PHYSICAL) {
    neighbor_idxs[neighbor_num] = conn.connections[dir][e].bid;
    neighbor_num += 1;
  }
  return neighbor_num;
}

/**
 *  Get the number of adjacent blocks above in a direction.
 * 
 * @param block_topo block topology object
 * @param bidx block index
 * @param direction
 * return number of neighbors
 */
int
get_num_above(struct gkyl_block_topo *block_topo, int bidx, int dir)
{
  int neighbor_idxs[1000] = {-1};
  int neighbor_num = get_neighbors(block_topo, bidx, dir, neighbor_idxs);
  return neighbor_num;
}




/**
 *  Get the indices of connected blocks in a direction.
 * 
 * @param block_topo block topology object
 * @param bidx block index
 * @param direction
 * @param block_list ordered indices of connected blocks including self
 * return number of connected blocks 
 */
int
get_connected(struct gkyl_block_topo *block_topo, int bidx, int dir, int* block_list)
{
  struct gkyl_block_connections conn;
  block_list[0] = bidx;
  int num_blocks = 1;
  int curr_bidx = bidx;
  while(true) {
    conn = block_topo->conn[curr_bidx];
    if (conn.connections[dir][0].edge == GKYL_PHYSICAL) {
      break;
    }
    else if (conn.connections[dir][0].edge == GKYL_UPPER_POSITIVE) { 
      if (conn.connections[dir][0].bid == bidx) return num_blocks;
      insert_below(block_list, num_blocks, conn.connections[dir][0].bid);
      curr_bidx = conn.connections[dir][0].bid;
      num_blocks+=1;
    }
  }

  curr_bidx = bidx;
  while(true) {
    conn = block_topo->conn[curr_bidx];
    if (conn.connections[dir][1].edge == GKYL_PHYSICAL) {
      break;
    }
    else if (conn.connections[dir][1].edge == GKYL_LOWER_POSITIVE) { 
      if (conn.connections[dir][1].bid == bidx) return num_blocks;
      insert_above(block_list, num_blocks, conn.connections[dir][1].bid);
      curr_bidx = conn.connections[dir][1].bid;
      num_blocks+=1;
    }
  }
  
  return num_blocks;
}

/**
 *  Get the number of connected blocks in a direction.
 * 
 * @param block_topo block topology object
 * @param bidx block index
 * @param direction
 * return number of connected blocks 
 */
int
get_num_connected(struct gkyl_block_topo *block_topo, int bidx, int dir)
{
  int block_list[1000] = {-1};
  int num_blocks = get_connected(block_topo, bidx, dir, block_list);
  return num_blocks;
}



/**
 * Check if a block corner is an interior corner
 * @param block_topo block topology object
 * @param bidx block index
 * @param edges list of edges (0 for lower, 1 for upper)
*/
int
check_corner(struct gkyl_block_topo *block_topo, int bidx, int* edges)
{
  struct gkyl_block_connections conn = block_topo->conn[bidx];
  int ndim = block_topo->ndim;
  int interior = 1; // true
  for (int i = 0; i < ndim; i++) {
    if(conn.connections[i][edges[i]].edge == GKYL_PHYSICAL)  {
      interior = 0;
      break;
    }
  }
  return interior;
}

/**
 * Get the list of blocks that touch a specific corner of a block
 * @param block_topo block topology object
 * @param bidx block index
 * @param edges list of edges of length ndim (0 for lower, 1 for upper)
 * @param block list on output a list of block indices that touch the corner
*/
int 
get_corner_connected(struct gkyl_block_topo *block_topo, int bidx, int* edges, int* block_list)
{
  int ndim = block_topo->ndim;
  int num_corner_connected = 0;
  struct gkyl_block_connections conn = block_topo->conn[bidx];
  int interior = check_corner(block_topo, bidx, edges);
  if (interior == 0) return num_corner_connected;
  num_corner_connected+=1;
  block_list[0] = bidx;
  int next_dir = 0;
  int next_edges[ndim];
  for (int i = 0; i < ndim; i++) next_edges[i] = edges[i];

  while(true) {
    int next_bidx = conn.connections[next_dir][next_edges[next_dir]].bid;
    if (next_bidx == bidx) break; // back at original
    block_list[num_corner_connected] = next_bidx;
    num_corner_connected+=1;
    next_edges[next_dir] = !next_edges[next_dir]; // 0 ->1 or 1 ->0
    next_dir = !next_dir;                                         
    interior = check_corner(block_topo, next_bidx, next_edges);
    if(interior == 0) break; // no more corners
    conn = block_topo->conn[next_bidx];
  }

  return num_corner_connected;
}

/**
 * Get the number of blocks that touch a specific corner of a block
 * @param block_topo block topology object
 * @param bidx block index
 * @param edges list of edges of length ndim (0 for lower, 1 for upper)
 * return number of blocks touching this corner
*/
int
get_num_corner_connected(struct gkyl_block_topo *block_topo, int bidx, int* edges)
{
  int block_list[1000] = {-1};
  int num_corner_connected = get_corner_connected(block_topo, bidx, edges, block_list);
  return num_corner_connected;
}


int
gkyl_multib_conn_get_connection(struct gkyl_block_topo *block_topo, int bidx, int dir, int corner_num, enum gkyl_conn_id conn_id, int *block_list)
{
  struct gkyl_block_connections conn = block_topo->conn[bidx];
  int num_connected = 0;

  if (conn_id == GKYL_CONN_NEIGHBOR) {
    num_connected = get_num_neighbors(block_topo, bidx, dir);
    get_neighbors(block_topo, bidx, dir, block_list);
  }
  else if (conn_id == GKYL_CONN_ALL) {
    num_connected = get_num_connected(block_topo, bidx, dir);
    get_connected(block_topo, bidx, dir, block_list);
  }
  else if (conn_id == GKYL_CONN_CORNER) {
    int e0 = corner_num/2;
    int e1 = corner_num%2;
    int edges[2] = {e0,e1};
    num_connected = get_num_corner_connected(block_topo, bidx, edges);
    get_corner_connected(block_topo, bidx, edges, block_list);
  }
  else if (conn_id == GKYL_CONN_BELOW) {
    num_connected = get_num_below(block_topo, bidx, dir);
    get_below(block_topo, bidx, dir, block_list);
  }
  else if (conn_id == GKYL_CONN_ABOVE) {
    num_connected = get_num_above(block_topo, bidx, dir);
    get_above(block_topo, bidx, dir, block_list);
  }

  return num_connected;
}

int
gkyl_multib_conn_get_num_connected(struct gkyl_block_topo *block_topo, int bidx, int dir, int corner_num, enum gkyl_conn_id conn_id)
{
  int block_list[1000] = {-1};
  int num_connected = gkyl_multib_conn_get_connection(block_topo, bidx, dir, corner_num, conn_id, block_list);
  return num_connected;
}

