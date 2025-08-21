/*
//@HEADER
// *****************************************************************************
//
// PULP: Multi-Objective Multi-Constraint Partitioning Using Label Propagation
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?  Contact  George M. Slota (gmslota@sandia.gov)
//                      Siva Rajamanickam (srajama@sandia.gov)
//
// *****************************************************************************
//@HEADER
*/


using namespace std;

extern int seed;

/*
########  ########  ######     #### ##    ## #### ########
##     ## ##       ##    ##     ##  ###   ##  ##     ##
##     ## ##       ##           ##  ####  ##  ##     ##
########  ######    ######      ##  ## ## ##  ##     ##
##     ## ##             ##     ##  ##  ####  ##     ##
##     ## ##       ##    ##     ##  ##   ###  ##     ##
########  ##        ######     #### ##    ## ####    ##
*/

// multi-source bfs with random start points
int* init_nonrandom(pulp_graph_t& g, int num_parts, int* parts)
{
  int  num_verts  = g.n;
  int* queue      = new int[num_verts*QUEUE_MULTIPLIER];
  int* queue_next = new int[num_verts*QUEUE_MULTIPLIER];   // The use of a QUEUE_MULTIPLIER of 2 in the code is likely a heuristic choice
  // to ensure that the queue has enough capacity to handle the vertices that
  // need to be processed during the breadth-first search (BFS) without running
  // out of space.
  int  queue_size = num_parts;
  int  next_size  = 0;

  #pragma omp parallel
  {
    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    xs1024star_t xs;
    xs1024star_seed((unsigned long)(seed + omp_get_thread_num()), &xs);

    // Initialize parts to -1 (not yet assigned to a part)
    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
      parts[i] = -1;

    // Randomly assign start points to parts
    #pragma omp single
    {
      for (int i = 0; i < num_parts; ++i)
      {
        int vert = (int)xs1024star_next(&xs) % num_verts;  // random start point
        while (parts[vert] != -1)
          vert = (int)xs1024star_next(&xs) % num_verts;    // random start point not already taken
        parts[vert] = i;                                   // assign start point to part i
        queue[i]    = vert;                                // add start point to queue
      }
    }

    while (queue_size > 0)
    {
      #pragma omp for schedule(guided) nowait
      for (int i = 0; i < queue_size; ++i)
      {
        int  vert       = queue[i];     // start point
        int  part       = parts[vert];  // part of start point
        long out_degree = out_degree(g, vert);
        int* outs       = out_vertices(g, vert);
        for (long j = 0; j < out_degree; ++j)
        {
          int out = outs[j];
          if (parts[out] == -1)
          {
            parts[out]                        = part; // assign part to out vertex
            thread_queue[thread_queue_size++] = out;

            if (thread_queue_size == THREAD_QUEUE_SIZE)
            {
              #pragma omp atomic capture
              thread_start = next_size += thread_queue_size;

              thread_start -= thread_queue_size;
              for (int l = 0; l < thread_queue_size; ++l)
                queue_next[thread_start+l] = thread_queue[l];
              thread_queue_size = 0;
            }
          }
        }
      }
      #pragma omp atomic capture
      thread_start = next_size += thread_queue_size;

      thread_start -= thread_queue_size;
      for (int l = 0; l < thread_queue_size; ++l)
        queue_next[thread_start+l] = thread_queue[l];
      thread_queue_size = 0;

      #pragma omp barrier

      #pragma omp single
      {
          int* temp  = queue;
          queue      = queue_next;
          queue_next = temp;

          queue_size = next_size;
          next_size  = 0;
      }
    }

    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
      if (parts[i] < 0)
        parts[i] = (unsigned)xs1024star_next(&xs) % num_parts;
  }

  #if OUTPUT_STEP
    evaluate_quality(g, num_parts, parts);
  #endif

  delete [] queue;
  delete [] queue_next;

  return parts;
}


/*
########  ########  ######     #### ##    ## #### ########     ######  #### ######## ########
##     ## ##       ##    ##     ##  ###   ##  ##     ##       ##    ##  ##       ##  ##
##     ## ##       ##           ##  ####  ##  ##     ##       ##        ##      ##   ##
########  ######    ######      ##  ## ## ##  ##     ##        ######   ##     ##    ######
##     ## ##             ##     ##  ##  ####  ##     ##             ##  ##    ##     ##
##     ## ##       ##    ##     ##  ##   ###  ##     ##       ##    ##  ##   ##      ##
########  ##        ######     #### ##    ## ####    ##        ######  #### ######## ########
*/

// multi-source bfs with random start points
// constrain to approximately equal part sizes, randomly assign any remaining.
int* init_nonrandom_constrained(pulp_graph_t& g, int num_parts, int* parts)
{
  bool has_vwgts = (g.vertex_weights != NULL);
  if (!has_vwgts)
  {
    printf("Error: vertex weights required for size-constrained init\n");
    exit(0);
  }

  int num_verts = g.n;

  int* queue      = new int[num_verts*QUEUE_MULTIPLIER];
  int* queue_next = new int[num_verts*QUEUE_MULTIPLIER];
  int* part_sizes = new int[num_parts];
  int  queue_size = num_parts;
  int  next_size  = 0;

  int max_partition_size = num_verts / num_parts * 2; // 2x average part size
  if (g.max_partition_size > 0)
    max_partition_size = g.max_partition_size; // max_partition_size is user-defined

  int rand_count = 0; //debug

  #pragma omp parallel
  {
    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    xs1024star_t xs;
    xs1024star_seed((unsigned long)(seed + omp_get_thread_num()), &xs);

    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
      parts[i] = -1;

    // Randomly assign start points to parts
    #pragma omp single
    {
      for (int i = 0; i < num_parts; ++i)
      {
        int vert = (int)xs1024star_next(&xs) % num_verts;
        while (parts[vert] != -1) {vert = xs1024star_next(&xs) % num_verts;}
        parts[vert]   = i;
        queue[i]      = vert;
        part_sizes[i] = g.vertex_weights[vert];
      }
    }

    // BFS from start points
    while (queue_size > 0)
    {
      #pragma omp for schedule(guided) nowait
      for (int i = 0; i < queue_size; ++i)
      {
        int  vert       = queue[i];               // start point
        int  part       = parts[vert];            // part of start point
        long out_degree = out_degree(g, vert);    // out degree of start point
        int* outs       = out_vertices(g, vert);  // out vertices of start point
        for (long j = 0; j < out_degree; ++j)
        {
          int out = outs[j];
          if (parts[out] == -1)
          {
            // If part size is less than max_partition_size (x2 average size),
            // assign it to the same part as the start point.
            // Otherwise, assign it to a random part
            if ((part_sizes[part] + g.vertex_weights[out]) < max_partition_size)
              parts[out] = part;
            else
            {
              int new_part = -1;
              for (int i = 0; i < num_parts; ++i)
              {
                new_part = (new_part + 1) % num_parts;
                if ((part_sizes[new_part] + g.vertex_weights[out]) < max_partition_size)
                {
                  parts[out] = new_part;
                  break;
                }
              }
              if (parts[out] == -1)
              {
                parts[out] = (int)((unsigned)xs1024star_next(&xs) % (unsigned)num_parts);
                printf("Random assignment (1) for vertex %d of size %d on part %d\n", out, g.vertex_weights[out], parts[out]);
                rand_count++;
              }
            }

            #pragma omp atomic
            part_sizes[parts[out]] += g.vertex_weights[out]; // increment part size

            // thread_queue is a local array within each thread, used to
            // temporarily store vertices that need to be processed in the next
            // iteration of the BFS. Each thread maintains its own thread_queue
            // to avoid race conditions when adding vertices to the queue.
            thread_queue[thread_queue_size++] = out;

            // Once thread_queue reaches its maximum size (THREAD_QUEUE_SIZE),
            // its contents are transferred to the global queue_next array.
            if (thread_queue_size == THREAD_QUEUE_SIZE)
            {
              // Ensures that each thread gets a unique and non-overlapping
              // segment of the queue_next array
              #pragma omp atomic capture
              thread_start = next_size += thread_queue_size;

              thread_start -= thread_queue_size;
              for (int l = 0; l < thread_queue_size; ++l)
                queue_next[thread_start+l] = thread_queue[l];
              thread_queue_size = 0;
            }
          } // end if part not assigned yet
        } // end for loop over out_degree
      } // end for loop over queue_size

      // Transfer remaining vertices from thread_queue to queue_next
      #pragma omp atomic capture
      thread_start = next_size += thread_queue_size;

      thread_start -= thread_queue_size;
      for (int l = 0; l < thread_queue_size; ++l)
        queue_next[thread_start+l] = thread_queue[l];
      thread_queue_size = 0;

      #pragma omp barrier

      // Swap queue and queue_next pointers and sizes after each iteration of
      // the BFS to process the next level
      #pragma omp single
      {
          int* temp  = queue;
          queue      = queue_next;
          queue_next = temp;

          queue_size = next_size;
          next_size  = 0;
      }
    } // end while

    // Assign any remaining vertices to random parts
    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
      if (parts[i] == -1)
        {
          int new_part = -1;
          for (int attempt = 0; attempt < RANDOM_BFS_TRIES; ++attempt)
          {
              new_part = (int)((unsigned)xs1024star_next(&xs) % (unsigned)num_parts);
              if ((part_sizes[new_part] + g.vertex_weights[i]) < max_partition_size)
              {
                  parts[attempt] = new_part;
                  break;
              }
          }

          // If no valid part was found after RANDOM_BFS_TRIES, assign to a random part
          if (parts[i] == -1)
          {
            parts[i] = (int)((unsigned)xs1024star_next(&xs) % (unsigned)num_parts);
            rand_count++;
          }
        }
  } // end parallel

  printf("Random assignments: %d\n", rand_count);

  #if OUTPUT_STEP
    evaluate_quality(g, num_parts, parts);
  #endif

  delete [] queue;
  delete [] queue_next;
  delete [] part_sizes;

  return parts;
}



/*
########  ########  ######      ######     ###    ########     ###     ######  #### ######## ##    ##    #### ##    ## #### ########
##     ## ##       ##    ##    ##    ##   ## ##   ##     ##   ## ##   ##    ##  ##     ##     ##  ##      ##  ###   ##  ##     ##
##     ## ##       ##          ##        ##   ##  ##     ##  ##   ##  ##        ##     ##      ####       ##  ####  ##  ##     ##
########  ######    ######     ##       ##     ## ########  ##     ## ##        ##     ##       ##        ##  ## ## ##  ##     ##
##     ## ##             ##    ##       ######### ##        ######### ##        ##     ##       ##        ##  ##  ####  ##     ##
##     ## ##       ##    ##    ##    ## ##     ## ##        ##     ## ##    ##  ##     ##       ##        ##  ##   ###  ##     ##
########  ##        ######      ######  ##     ## ##        ##     ##  ######  ####    ##       ##       #### ##    ## ####    ##
*/

// multi-source bfs with random start points
// part with heterogeneous capacities, randomly assign any remaining.
int*
init_nonrandom_constrained_capacity(pulp_graph_t& g, int num_parts, int* parts, int vertex_balance)
{
  int rand_count = 0; //debug

  bool has_vwgts        = (g.vertex_weights != NULL);
  bool has_p_capacities = (g.partition_capacities != NULL);

  if (!has_vwgts || !has_p_capacities)
  {
    printf("ERROR: Vertex weights and partition capacities are required for this function\n");
    return NULL;
  }

  int num_verts = g.n;

  int* queue      = new int[num_verts*QUEUE_MULTIPLIER];
  int* queue_next = new int[num_verts*QUEUE_MULTIPLIER];
  int* part_sizes = new int[num_parts];
  int  queue_size = num_parts;
  int  next_size  = 0;

  // Vector of maximum partition sizes
  int* max_partition_size = new int[num_parts];
  if (g.max_partition_size > 0)
    for (int i = 0; i < num_parts; ++i)
    {
      max_partition_size[i] =  vertex_balance * g.max_partition_size * g.partition_capacities[i];
      printf("max_partition_size[%2d] = %d\n", i, max_partition_size[i]);
    }
  else
  {
    int unit_capacity_wgts = (int)g.vertex_weights_sum / (int)g.partition_capacities_sum;
    for (int i = 0; i < num_parts; ++i)
      max_partition_size[i] = vertex_balance * g.partition_capacities[i] * unit_capacity_wgts;
  }

  #pragma omp parallel
  {
    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    xs1024star_t xs;
    xs1024star_seed((unsigned long)(seed + omp_get_thread_num()), &xs);

    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
      parts[i] = -1;

    // Randomly assign start points to parts
    #pragma omp single
    {
      for (int i = 0; i < num_parts; ++i)
      {
        int vert = (int)xs1024star_next(&xs) % num_verts;
        while (parts[vert] != -1) {vert = xs1024star_next(&xs) % num_verts;}
        parts[vert]   = i;
        queue[i]      = vert;
        part_sizes[i] = g.vertex_weights[vert];
        // Partition size is the sum of vertex weights. Partition capacity will be checked later.
      }
    }

    // BFS from start points
    while (queue_size > 0)
    {
      #pragma omp for schedule(guided) nowait
      for (int i = 0; i < queue_size; ++i)
      {
        int  vert       = queue[i];               // start point
        int  part       = parts[vert];            // part of start point
        long out_degree = out_degree(g, vert);    // out degree of start point
        int* outs       = out_vertices(g, vert);  // out vertices of start point
        for (long j = 0; j < out_degree; ++j)
        {
          int out = outs[j];
          if (parts[out] == -1) // If part not assigned yet
          {
            // If part size is less than max_partition_size, assign it to the same part as the start point.
            // Otherwise, assign it to a random part
            if ((part_sizes[part] + g.vertex_weights[out]) < max_partition_size[part])
                parts[out] = part;
            else
            {
              int new_part = -1;
              for (int i = 0; i < num_parts; ++i)
              {
                new_part = (new_part + 1) % num_parts;
                if ((part_sizes[new_part] + g.vertex_weights[out]) < max_partition_size[new_part])
                {
                  parts[out] = new_part;
                  break;
                }
              }
              if (parts[out] == -1)
              {
                parts[out] = (int)((unsigned)xs1024star_next(&xs) % (unsigned)num_parts);

                printf("Random assignment (1) for vertex %d of size %d on part %d\n", out, g.vertex_weights[out], parts[out]);
                #pragma omp atomic
                rand_count++;
              }
            }

            // Update the part size with the vertex weight of the new vertex
            #pragma omp atomic
            part_sizes[parts[out]] += g.vertex_weights[out];

            // thread_queue is a local array within each thread, used to
            // temporarily store vertices that need to be processed in the
            // next iteration of the BFS. Each thread maintains its own thread_queue
            // to avoid race conditions when adding vertices to the queue.
            thread_queue[thread_queue_size++] = out;

            // Once thread_queue reaches its maximum size (THREAD_QUEUE_SIZE),
            // its contents are transferred to the global queue_next array.
            if (thread_queue_size == THREAD_QUEUE_SIZE)
            {
              // Ensures that each thread gets a unique and non-overlapping
              // segment of the queue_next array
              #pragma omp atomic capture
              thread_start = next_size += thread_queue_size;

              thread_start -= thread_queue_size;
              for (int l = 0; l < thread_queue_size; ++l)
                queue_next[thread_start+l] = thread_queue[l];
              thread_queue_size = 0;
            }
          } // end if part not assigned yet
        } // end for loop over out_degree
      } // end for loop over queue_size

      #pragma omp atomic capture
      thread_start = next_size += thread_queue_size;

      thread_start -= thread_queue_size;
      for (int l = 0; l < thread_queue_size; ++l)
        queue_next[thread_start+l] = thread_queue[l];
      thread_queue_size = 0;

      #pragma omp barrier

      #pragma omp single
      {
          int* temp  = queue;
          queue      = queue_next;
          queue_next = temp;

          queue_size = next_size;
          next_size  = 0;
      }
    } // end while loop over queue_size

    // Assign any remaining vertices to random parts
    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
      if (parts[i] == -1)
      {
        int new_part = -1;
        for (int attempt = 0; attempt < RANDOM_BFS_TRIES; ++attempt)
        {
            new_part = (new_part + 1) % num_parts; //(int)((unsigned)xs1024star_next(&xs) % (unsigned)num_parts);
            if ((part_sizes[new_part] + g.vertex_weights[i]) < max_partition_size[new_part])
            {
                parts[attempt] = new_part;
                break;
            }
        }

        // If no valid part was found after RANDOM_BFS_TRIES, assign to a random part
        if (parts[i] == -1)
        {
          parts[i] = (int)((unsigned)xs1024star_next(&xs) % (unsigned)num_parts);
          printf("Random assignment (2) for vertex %d of size %d on part %d\n", i, g.vertex_weights[i], parts[i]);
          #pragma omp atomic
          rand_count++;
        }
      }
  } // end parallel

  printf("Random assignments: %d\n", rand_count);

  #if OUTPUT_STEP
    evaluate_quality(g, num_parts, parts);
  #endif

  delete [] queue;
  delete [] queue_next;
  delete [] part_sizes;

  return parts;
}