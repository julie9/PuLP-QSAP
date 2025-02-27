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

#include <cassert>

using namespace std;

extern int seed;

/*
'########::'########:::'#######::'########::
 ##.... ##: ##.... ##:'##.... ##: ##.... ##:
 ##:::: ##: ##:::: ##: ##:::: ##: ##:::: ##:
 ########:: ########:: ##:::: ##: ########::
 ##.....::: ##.. ##::: ##:::: ##: ##.....:::
 ##:::::::: ##::. ##:: ##:::: ##: ##::::::::
 ##:::::::: ##:::. ##:. #######:: ##::::::::
..:::::::::..:::::..:::.......:::..:::::::::
*/
int* label_prop(pulp_graph_t& g, int num_parts, int* parts,
  int label_prop_iter, double balance_vert_lower)
{
  int  num_verts  = g.n;
  int* part_sizes = new int[num_parts];
  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  int num_changes;
  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int   queue_size    = num_verts;
  int   next_size     = 0;

  double avg_size = (double)num_verts / (double)num_parts;
  double min_size = avg_size * balance_vert_lower;

  #pragma omp parallel
  {
    xs1024star_t xs;
    xs1024star_seed((unsigned long)(seed + omp_get_thread_num()), &xs);

    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
      parts[i] = (int)((unsigned)xs1024star_next(&xs) % (unsigned)num_parts);

    long* part_sizes_thread = new long[num_parts];
    for (int i = 0; i < num_parts; ++i)
      part_sizes_thread[i] = 0;

    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
        ++part_sizes_thread[parts[i]];

    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_sizes[i] += part_sizes_thread[i];

    delete [] part_sizes_thread;


    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      queue[i] = i;
    #pragma omp for schedule(static)
    for (int i = 0; i < num_verts; ++i)
      in_queue_next[i] = false;

    int* part_counts = new int[num_parts];
    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    for (int num_iter = 0; num_iter < label_prop_iter; ++num_iter)
    {
      num_changes = 0;

      #pragma omp for schedule(guided) reduction(+:num_changes)
      for (int i = 0; i < queue_size; ++i)
      {
        int v = queue[i];
        in_queue[v] = false;
        for (int j = 0; j < num_parts; ++j)
          part_counts[j] = 0;

        unsigned out_degree = out_degree(g, v);
        int* outs = out_vertices(g, v);
        for (unsigned j = 0; j < out_degree; ++j)
        {
          int out = outs[j];
          int part = parts[out];
          part_counts[part] += out_degree(g, out);
        }

        int part = parts[v];
        int max_count = -1;
        int max_part = -1;
        int num_max = 0;
        for (int p = 0; p < num_parts; ++p)
        {
          if (part_counts[p] == max_count &&
              (part_sizes[p]-1) > (int)min_size)
          {
            part_counts[num_max++] = p;
          }
          else if (part_counts[p] > max_count &&
                  (part_sizes[p]-1) > (int)min_size)
          {
            max_count = part_counts[p];
            max_part = p;
            num_max = 0;
            part_counts[num_max++] = p;
          }
        }

        if (num_max > 1)
          max_part = part_counts[(int)rand() % num_max];

        if (max_part != part &&
            (part_sizes[part]-1) > (int)min_size)
        {
          #pragma omp atomic
          ++part_sizes[max_part];
          #pragma omp atomic
          --part_sizes[part];

          parts[v] = max_part;
          ++num_changes;

          if (!in_queue_next[v])
          {
            in_queue_next[v] = true;
            thread_queue[thread_queue_size++] = v;

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
          for (unsigned j = 0; j < out_degree; ++j)
          {
            if (!in_queue_next[outs[j]])
            {
              in_queue_next[outs[j]] = true;
              thread_queue[thread_queue_size++] = outs[j];

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
      }
      #pragma omp atomic capture
      thread_start = next_size += thread_queue_size;

      thread_start -= thread_queue_size;
      for (int l = 0; l < thread_queue_size; ++l)
        queue_next[thread_start+l] = thread_queue[l];
      thread_queue_size = 0;

      #pragma omp barrier

      ++num_iter;
      #pragma omp single
      {
        #if VERBOSE
          printf("%d\n", next_size);
        #endif

        int* temp = queue;
        queue = queue_next;
        queue_next = temp;
        bool* temp_b = in_queue;
        in_queue = in_queue_next;
        in_queue_next = temp_b;

        queue_size = next_size;
        next_size = 0;

        #if OUTPUT_STEP
          evaluate_quality(g, num_parts, parts);
        #endif
      } // end single
    } // end while

    delete [] part_counts;
  } // end parallel

  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;

  return parts;
}


/*
'##:::::'##::::'########::'########:::'#######::'########::
 ##:'##: ##:::: ##.... ##: ##.... ##:'##.... ##: ##.... ##:
 ##: ##: ##:::: ##:::: ##: ##:::: ##: ##:::: ##: ##:::: ##:
 ##: ##: ##:::: ########:: ########:: ##:::: ##: ########::
 ##: ##: ##:::: ##.....::: ##.. ##::: ##:::: ##: ##.....:::
 ##: ##: ##:::: ##:::::::: ##::. ##:: ##:::: ##: ##::::::::
. ###. ###::::: ##:::::::: ##:::. ##:. #######:: ##::::::::
:...::...::::::..:::::::::..:::::..:::.......:::..:::::::::
*/
int* label_prop_weighted(pulp_graph_t& g, int num_parts, int* parts,
  int label_prop_iter, double balance_vert_lower)
{
  int  num_verts = g.n;
  bool has_vwgts = (g.vertex_weights != NULL);
  bool has_ewgts = (g.edge_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;

  int* part_sizes = new int[num_parts];
  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  int num_changes;
  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int   queue_size    = num_verts;
  int   next_size     = 0;

  double avg_size = (double)g.vertex_weights_sum / (double)num_parts;
  double min_size = avg_size * balance_vert_lower;

  #pragma omp parallel
  {
    xs1024star_t xs;
    xs1024star_seed((unsigned long)(seed + omp_get_thread_num()), &xs);

    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
      parts[i] = (int)((unsigned)xs1024star_next(&xs) % (unsigned)num_parts);

    long* part_sizes_thread = new long[num_parts];
    for (int i = 0; i < num_parts; ++i)
      part_sizes_thread[i] = 0;

    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      if (has_vwgts)
        part_sizes_thread[parts[i]] += g.vertex_weights[i];
      else
        ++part_sizes_thread[parts[i]];

    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_sizes[i] += part_sizes_thread[i];

    delete [] part_sizes_thread;


    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      queue[i] = i;
    #pragma omp for schedule(static)
    for (int i = 0; i < num_verts; ++i)
      in_queue_next[i] = false;

    int* part_counts = new int[num_parts];
    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    for (int num_iter = 0; num_iter < label_prop_iter; ++num_iter)
    {
      num_changes = 0;

      #pragma omp for schedule(guided) reduction(+:num_changes)
      for (int i = 0; i < queue_size; ++i)
      {
        int v = queue[i];
        int v_weight = 1;
        if (has_vwgts) v_weight = g.vertex_weights[v];

        in_queue[v] = false;
        for (int j = 0; j < num_parts; ++j)
          part_counts[j] = 0;

        unsigned out_degree = out_degree(g, v);
        int* outs = out_vertices(g, v);
        int* weights = out_weights(g, v);
        for (unsigned j = 0; j < out_degree; ++j)
        {
          int out = outs[j];
          int part_out = parts[out];
          double weight_out = 1.0;
          if (has_ewgts) weight_out = (double)weights[j];
          part_counts[part_out] += (double)out_degree(g, out)*weight_out;
        }

        int part = parts[v];
        int max_count = -1;
        int max_part = -1;
        int num_max = 0;
        for (int p = 0; p < num_parts; ++p)
        {
          if (part_counts[p] == max_count)
          {
            part_counts[num_max++] = p;
          }
          else if (part_counts[p] > max_count)
          {
            max_count = part_counts[p];
            max_part = p;
            num_max = 0;
            part_counts[num_max++] = p;
          }
        }

        if (num_max > 1)
          max_part = part_counts[(int)rand() % num_max];

        if (max_part != part &&
            (part_sizes[part] - v_weight > (int)min_size))
        {
          #pragma omp atomic
          part_sizes[max_part] += v_weight;
          #pragma omp atomic
          part_sizes[part] -= v_weight;

          parts[v] = max_part;
          ++num_changes;

          if (!in_queue_next[v])
          {
            in_queue_next[v] = true;
            thread_queue[thread_queue_size++] = v;

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
          for (unsigned j = 0; j < out_degree; ++j)
          {
            if (!in_queue_next[outs[j]])
            {
              in_queue_next[outs[j]] = true;
              thread_queue[thread_queue_size++] = outs[j];

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
      }
      #pragma omp atomic capture
      thread_start = next_size += thread_queue_size;

      thread_start -= thread_queue_size;
      for (int l = 0; l < thread_queue_size; ++l)
        queue_next[thread_start+l] = thread_queue[l];
      thread_queue_size = 0;

      #pragma omp barrier

      ++num_iter;

      #pragma omp single
      {
        #if VERBOSE
          printf("%d\n", next_size);
        #endif

        int*  temp = queue;
        queue      = queue_next;
        queue_next = temp;
        bool* temp_b  = in_queue;
        in_queue      = in_queue_next;
        in_queue_next = temp_b;

        queue_size = next_size;
        next_size  = 0;

        #if OUTPUT_STEP
          evaluate_quality(g, num_parts, parts);
        #endif
      } // end single
    } // end while

    delete [] part_counts;
  } // end parallel

  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;

  return parts;
}




/*
##      ##    #### ##    ## ######## ######## ########     ##          ###    ########  ######## ##          ########  ########   #######  ########
##  ##  ##     ##  ###   ##    ##    ##       ##     ##    ##         ## ##   ##     ## ##       ##          ##     ## ##     ## ##     ## ##     ##
##  ##  ##     ##  ####  ##    ##    ##       ##     ##    ##        ##   ##  ##     ## ##       ##          ##     ## ##     ## ##     ## ##     ##
##  ##  ##     ##  ## ## ##    ##    ######   ########     ##       ##     ## ########  ######   ##          ########  ########  ##     ## ########
##  ##  ##     ##  ##  ####    ##    ##       ##   ##      ##       ######### ##     ## ##       ##          ##        ##   ##   ##     ## ##
##  ##  ##     ##  ##   ###    ##    ##       ##    ##     ##       ##     ## ##     ## ##       ##          ##        ##    ##  ##     ## ##
 ###  ###     #### ##    ##    ##    ######## ##     ##    ######## ##     ## ########  ######## ########    ##        ##     ##  #######  ##
*/
int*
label_prop_weighted_interpart(pulp_graph_t& g, int num_parts, int* parts, int label_prop_iter,
                              double balance_vert_lower)
{
  bool has_ipwgts = (g.interpartition_weights != NULL);
  if (!has_ipwgts)
  {
    printf("Error: interpartition weights required for weighted interpart label prop\n");
    exit(0);
  }

  int  num_verts = g.n;
  bool has_vwgts = (g.vertex_weights != NULL);
  bool has_ewgts = (g.edge_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;

  int* part_sizes = new int[num_parts];
  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  int   num_changes   = 0;
  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int   queue_size    = num_verts;
  int   next_size     = 0;

  double avg_size = (double)g.vertex_weights_sum / (double)num_parts;
  double min_size = avg_size * balance_vert_lower;

  #pragma omp parallel
  {
    // =====================================================
    // Initialize parts randomly
    // =====================================================
    xs1024star_t xs;
    xs1024star_seed((unsigned long)(seed + omp_get_thread_num()), &xs);

    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
      parts[i] = (int)((unsigned)xs1024star_next(&xs) % (unsigned)num_parts);

    // Compute the size of each partition
    long* part_sizes_thread = new long[num_parts];
    for (int i = 0; i < num_parts; ++i)
      part_sizes_thread[i] = 0;

    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      if (has_vwgts)
        part_sizes_thread[parts[i]] += g.vertex_weights[i];
      else
        ++part_sizes_thread[parts[i]];

    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_sizes[i] += part_sizes_thread[i];

    delete [] part_sizes_thread;


    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      queue[i] = i;

    #pragma omp for schedule(static)
    for (int i = 0; i < num_verts; ++i)
      in_queue_next[i] = false;

    int* part_counts = new int[num_parts];
    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    // =====================================================
    // Outer loop for label propagation
    // =====================================================
    for (int num_iter = 0; num_iter < label_prop_iter; ++num_iter)
    {
      num_changes = 0;

      // =====================================================
      // Perform label propagation iterations for all vertices
      // =====================================================
      #pragma omp for schedule(guided) reduction(+:num_changes)
      for (int i = 0; i < queue_size; ++i)
      {
        int v    = queue[i]; // For each vertex in the queue (i.e. all vertices)
        int part = parts[v]; // Get the partition of the vertex

        // -----------------------------------------------------
        // Vertex weight
        // -----------------------------------------------------
        int v_weight = 1; // Get the weight of the vertex
        if (has_vwgts) v_weight = g.vertex_weights[v];

        in_queue[v] = false;
        for (int j = 0; j < num_parts; ++j)
          part_counts[j] = 0;

        unsigned out_degree = out_degree(g, v);   // Get the outgoing degree of the vertex
        int*     outs       = out_vertices(g, v); // Get the outgoing vertices (id of the vertices)
        int*     weights    = out_weights(g, v);  // Get the weights of the outgoing edges
        for (unsigned j = 0; j < out_degree; ++j)
        {
          int out      = outs[j];    // Get the id of one of the outgoing vertices
          int part_out = parts[out]; // Get the partition of the outgoing vertex
          double weight_out = 1.0;
          if (has_ewgts) weight_out = (double)weights[j];
          part_counts[part_out] += (double)out_degree(g, out)*weight_out;
        }

        // -----------------------------------------------------
        // Inter-partition weight
        // -----------------------------------------------------
        int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
        for (int p = 0; p < num_parts; ++p)
            part_counts[p] *= partition_comm_weights[p];

        // -----------------------------------------------------
        // Check which partition has the maximum count of outgoing edges
        // -----------------------------------------------------
        int max_count = -1;
        int max_part  = -1;
        int num_max   =  0; // Number of partitions with the maximum number of edges
        for (int p = 0; p < num_parts; ++p)
        {
          if (part_counts[p] == max_count)
          {
            part_counts[num_max++] = p;
          }
          else if (part_counts[p] > max_count)
          {
            max_count = part_counts[p];
            max_part  = p;
            num_max   = 0;
            part_counts[num_max++] = p;
          }
        }
        // If there are multiple partitions with the maximum count, randomly select one
        if (num_max > 1)
          max_part = part_counts[(int)rand() % num_max];

        // -----------------------------------------------------
        // Swap the vertex to the partition with the maximum count
        // -----------------------------------------------------
        if (max_part != part)
        {

          if (!g.do_bin_packing &&
              (part_sizes[part] - v_weight < (int)min_size))
            continue;

          parts[v] = max_part; // Move vertex v to the partition with the maximum count
          ++num_changes;

          // Partition weight of the partition with the maximum count
          #pragma omp atomic
          part_sizes[max_part] += v_weight;
          #pragma omp atomic
          part_sizes[part]     -= v_weight;

          if (!in_queue_next[v])
          {
            in_queue_next[v] = true;
            thread_queue[thread_queue_size++] = v;

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

        // Add the neighbors of vertex v to the queue
          for (unsigned j = 0; j < out_degree; ++j)
          {
            if (!in_queue_next[outs[j]])
            {
              in_queue_next[outs[j]] = true;
              thread_queue[thread_queue_size++] = outs[j];

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
      }

      #pragma omp atomic capture
      thread_start = next_size += thread_queue_size;

      thread_start -= thread_queue_size;
      for (int l = 0; l < thread_queue_size; ++l)
        queue_next[thread_start+l] = thread_queue[l];
      thread_queue_size = 0;

      #pragma omp barrier

      // ++num_iter; // Notes(julie9): removed, already incremented in the for loop

      #pragma omp single
      {
        #if VERBOSE
          printf("%d\n", next_size);
        #endif

        int*  temp    = queue;
        queue         = queue_next;
        queue_next    = temp;
        bool* temp_b  = in_queue;
        in_queue      = in_queue_next;
        in_queue_next = temp_b;

        queue_size = next_size;
        next_size = 0;

        #if OUTPUT_STEP
          evaluate_quality(g, num_parts, parts);
        #endif
      } // end single
    } // end for num_iter

    delete [] part_counts;
  } // end parallel

  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;

  return parts;
}






/*
##      ##     ######     ###    ########     #### ##    ## ######## ######## ########     ##          ###    ########  ######## ##          ########  ########   #######  ########
##  ##  ##    ##    ##   ## ##   ##     ##     ##  ###   ##    ##    ##       ##     ##    ##         ## ##   ##     ## ##       ##          ##     ## ##     ## ##     ## ##     ##
##  ##  ##    ##        ##   ##  ##     ##     ##  ####  ##    ##    ##       ##     ##    ##        ##   ##  ##     ## ##       ##          ##     ## ##     ## ##     ## ##     ##
##  ##  ##    ##       ##     ## ########      ##  ## ## ##    ##    ######   ########     ##       ##     ## ########  ######   ##          ########  ########  ##     ## ########
##  ##  ##    ##       ######### ##            ##  ##  ####    ##    ##       ##   ##      ##       ######### ##     ## ##       ##          ##        ##   ##   ##     ## ##
##  ##  ##    ##    ## ##     ## ##            ##  ##   ###    ##    ##       ##    ##     ##       ##     ## ##     ## ##       ##          ##        ##    ##  ##     ## ##
 ###  ###      ######  ##     ## ##           #### ##    ##    ##    ######## ##     ##    ######## ##     ## ########  ######## ########    ##        ##     ##  #######  ##
 */
int*
label_prop_weighted_interpart_capacity(pulp_graph_t& g, int num_parts, int* parts, int label_prop_iter,
                                       double balance_vert_lower)
{
  bool has_ipwgts       = (g.interpartition_weights != NULL);
  bool has_p_capacities = (g.partition_capacities != NULL);
  if (!has_p_capacities || !has_ipwgts)
  {
    printf("Error: partition capacities required for weighted interpart label prop\n");
    exit(0);
  }

  int  num_verts = g.n;
  bool has_vwgts = (g.vertex_weights != NULL);
  bool has_ewgts = (g.edge_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;

  int* part_sizes = new int[num_parts];
  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  int   num_changes   = 0;
  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int   queue_size    = num_verts;
  int   next_size     = 0;

  double  unit_avg_size = (double)g.vertex_weights_sum / (double)g.partition_capacities_sum;

  // TODO(julie9): Update based on the version in random init.
  double* min_sizes          = new double[num_parts];
  for (int i = 0; i < num_parts; ++i)
    min_sizes[i] = unit_avg_size * g.partition_capacities[i] * balance_vert_lower;

  // Compute the cumulative distribution of partition capacities
  double* cumulative_distribution = nullptr;
  cumulative_distribution    = new double[num_parts];  // Initialize the cumulative distribution
  cumulative_distribution[0] = (double)g.partition_capacities[0] / g.partition_capacities_sum;
  for (int i = 1; i < num_parts; ++i)
    cumulative_distribution[i] = cumulative_distribution[i - 1]
                                + (double)g.partition_capacities[i] / g.partition_capacities_sum;

  #pragma omp parallel
  {
    // =====================================================
    // Initialize parts randomly
    // =====================================================
    xs1024star_t xs;
    xs1024star_seed((unsigned long)(seed + omp_get_thread_num()), &xs);

    #pragma omp for
    for (int i = 0; i < num_verts; ++i)
    {
      double rand_val = (double)xs1024star_next(&xs) / (double)UINT64_MAX; // Random value between 0 and 1
      for (int j = 0; j < num_parts; ++j)
        if (rand_val <= cumulative_distribution[j])
        {
          parts[i] = j; // Assign the vertex to the partition
          break;        // Break the loop after assigning the vertex to the partition
        }
    }

    // Compute the size of each partition
    long* part_sizes_thread = new long[num_parts];
    for (int i = 0; i < num_parts; ++i)
      part_sizes_thread[i] = 0;

    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      if (has_vwgts)
        part_sizes_thread[parts[i]] += g.vertex_weights[i];
      else
        ++part_sizes_thread[parts[i]];

    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_sizes[i] += part_sizes_thread[i];

    delete [] part_sizes_thread;


    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      queue[i] = i;

    #pragma omp for schedule(static)
    for (int i = 0; i < num_verts; ++i)
      in_queue_next[i] = false;

    int* part_counts = new int[num_parts];
    int* max_parts   = new int[num_parts]; // Array to store the partitions with the maximum number of edges

    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    // =====================================================
    // Outer loop for label propagation
    // =====================================================
    for (int num_iter = 0; num_iter < label_prop_iter; ++num_iter)
    {
      num_changes = 0;

      // =====================================================
      // Perform label propagation iterations for all vertices
      // =====================================================
      #pragma omp for schedule(guided) reduction(+:num_changes)
      for (int i = 0; i < queue_size; ++i)
      {
        int v    = queue[i]; // For each vertex in the queue (i.e. all vertices)
        int part = parts[v]; // Get the partition of the vertex

        // -----------------------------------------------------
        // Vertex weight
        // -----------------------------------------------------
        int v_weight = 1; // Get the weight of the vertex
        if (has_vwgts) v_weight = g.vertex_weights[v];

        in_queue[v] = false;
        for (int j = 0; j < num_parts; ++j)
          part_counts[j] = 0;

        unsigned out_degree = out_degree(g, v);   // Get the outgoing degree of the vertex
        int*     outs       = out_vertices(g, v); // Get the outgoing vertices (id of the vertices)
        int*     weights    = out_weights(g, v);  // Get the weights of the outgoing edges
        for (unsigned j = 0; j < out_degree; ++j)
        {
          int out      = outs[j];    // Get the id of one of the outgoing vertices
          int part_out = parts[out]; // Get the partition of the outgoing vertex
          double weight_out = 1.0;
          if (has_ewgts) weight_out = (double)weights[j];
          part_counts[part_out] += (double)out_degree(g, out)*weight_out;
        }

        // -----------------------------------------------------
        // Inter-partition weight
        // -----------------------------------------------------
        int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
        for (int p = 0; p < num_parts; ++p)
            part_counts[p] *= partition_comm_weights[p];

        // -----------------------------------------------------
        // Check which partition has the maximum count of outgoing edges
        // -----------------------------------------------------
        int max_count = -1;
        int max_part  = -1;
        int num_max   =  0; // Number of partitions with the maximum number of edges
        for (int p = 0; p < num_parts; ++p)
        {
          if (part_counts[p] == max_count)
            max_parts[num_max++] = p;
          else if (part_counts[p] > max_count)
          {
            max_count = part_counts[p];
            max_part  = p;
            num_max   = 0;
            // part_counts[num_max++] = p;
            max_parts[num_max++] = p;
          }
        }
        // If there are multiple partitions with the maximum count, randomly select one
        if (num_max > 1)
          max_part = max_parts[(int)rand() % num_max];

        // -----------------------------------------------------
        // Swap the vertex to the partition with the maximum count
        // -----------------------------------------------------
        if (max_part != part &&
            part_sizes[part] - v_weight > (int)min_sizes[part])
        {
          parts[v] = max_part; // Move vertex v to the partition with the maximum count
          ++num_changes;

          // Partition weight of the partition with the maximum count
          #pragma omp atomic
          part_sizes[max_part] += v_weight;
          #pragma omp atomic
          part_sizes[part]     -= v_weight;

          if (!in_queue_next[v])
          {
            in_queue_next[v] = true;
            thread_queue[thread_queue_size++] = v;

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

          // Add the neighbors of vertex v to the queue
          for (unsigned j = 0; j < out_degree; ++j)
          {
            if (!in_queue_next[outs[j]])
            {
              in_queue_next[outs[j]] = true;
              thread_queue[thread_queue_size++] = outs[j];

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
      }

      #pragma omp atomic capture
      thread_start = next_size += thread_queue_size;

      thread_start -= thread_queue_size;
      for (int l = 0; l < thread_queue_size; ++l)
        queue_next[thread_start+l] = thread_queue[l];
      thread_queue_size = 0;

      #pragma omp barrier

      // ++num_iter; // Notes(julie9): removed, already incremented in the for loop

      #pragma omp single
      {
        #if VERBOSE
          printf("%d\n", next_size);
        #endif

        int*  temp    = queue;
        queue         = queue_next;
        queue_next    = temp;
        bool* temp_b  = in_queue;
        in_queue      = in_queue_next;
        in_queue_next = temp_b;

        queue_size = next_size;
        next_size = 0;

        #if OUTPUT_STEP
          evaluate_quality(g, num_parts, parts);
        #endif
      } // end single
    } // end for num_iter

    delete [] part_counts;
    delete [] max_parts;
  } // end parallel

  delete [] min_sizes;
  delete [] cumulative_distribution;

  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;

  return parts;
}

