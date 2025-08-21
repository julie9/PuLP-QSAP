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

/*
'########:::::'###::::'##::::::::::'##::::'##:'########:'########::'########:
 ##.... ##:::'## ##::: ##:::::::::: ##:::: ##: ##.....:: ##.... ##:... ##..::
 ##:::: ##::'##:. ##:: ##:::::::::: ##:::: ##: ##::::::: ##:::: ##:::: ##::::
 ########::'##:::. ##: ##:::::::::: ##:::: ##: ######::: ########::::: ##::::
 ##.... ##: #########: ##::::::::::. ##:: ##:: ##...:::: ##.. ##:::::: ##::::
 ##:::: ##: ##.... ##: ##:::::::::::. ## ##::: ##::::::: ##::. ##::::: ##::::
 ########:: ##:::: ##: ########::::::. ###:::: ########: ##:::. ##:::: ##::::
........:::..:::::..::........::::::::...:::::........::..:::::..:::::..:::::
*/

void label_balance_verts(pulp_graph_t& g, int num_parts, int* parts,
  int vert_outer_iter, int vert_balance_iter, int vert_refine_iter,
  double vert_balance)
{
  int  num_verts  = g.n;
  int* part_sizes = new int[num_parts];

  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  double avg_size      = num_verts / num_parts;
  int    num_swapped_1 = 0;
  int    num_swapped_2 = 0;
  double max_v;
  double running_max_v = (double)num_verts;

  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int   queue_size    = num_verts;
  int   next_size     = 0;
  int   t             = 0;
  int   num_tries     = 0;

  #pragma omp parallel
  {
    int* part_sizes_thread = new int[num_parts];
    for (int i = 0; i < num_parts; ++i)
      part_sizes_thread[i] = 0;

    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      ++part_sizes_thread[parts[i]];

    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_sizes[i] += part_sizes_thread[i];

    delete [] part_sizes_thread;

    double* part_counts  = new double[num_parts];
    double* part_weights = new double[num_parts];

    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    for (int p = 0; p < num_parts; ++p)
    {
      if (part_sizes[p] == 0)
        part_weights[p] = numeric_limits<double>::max();
      else
      {
        part_weights[p] = (vert_balance * avg_size / (double)part_sizes[p]) - 1.0;
        if (part_weights[p] < 0.0)
          part_weights[p] = 0.0;
      }
    }

    while(t < vert_outer_iter)
    {

      #pragma omp for schedule(static) nowait
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;
      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        in_queue_next[i] = false;

      #pragma omp single
      {
        num_swapped_1 = 0;
        queue_size    = num_verts;
        next_size     = 0;
      }

      int num_iter = 0;
      while (/*swapped &&*/ num_iter < vert_balance_iter)
      {
        #pragma omp for schedule(guided) reduction(+:num_swapped_1) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v    = queue[i];
          int part = parts[v];
          in_queue[v] = false;
          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0.0;

          unsigned out_degree = out_degree(g, v);
          int* outs = out_vertices(g, v);
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out      = outs[j];
            int part_out = parts[out];
            part_counts[part_out] += out_degree(g, out);
            //part_counts[part_out] += 1.0;//out_degree(g, out);
          }

          int max_part = part;
          double max_val = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            part_counts[p] *= part_weights[p];

            if (part_counts[p] > max_val)
            {
              max_val = part_counts[p];
              max_part = p;
            }
          }

          if (max_part != part)
          {
            parts[v] = max_part;
            ++num_swapped_1;
            #pragma omp atomic
            --part_sizes[part];
            #pragma omp atomic
            ++part_sizes[max_part];

            part_weights[part] = vert_balance * avg_size / (double)part_sizes[part] - 1.0;
            part_weights[max_part] = vert_balance * avg_size / (double)part_sizes[max_part]  - 1.0;

            if (part_weights[part] < 0.0)
              part_weights[part] = 0.0;
            if (part_weights[max_part] < 0.0)
              part_weights[max_part] = 0.0;

            if (!in_queue_next[v])
            {
              in_queue_next[v]                  = true;
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
                in_queue_next[outs[j]]            = true;
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
				    printf("num_swapped_1 (vert, balance): %d\n", num_swapped_1);
          #endif
          int*  temp    = queue;
          queue         = queue_next;
          queue_next    = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          num_swapped_1 = 0;

          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        }
      } // end while

      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;

      #pragma omp single
      {
        num_swapped_2 = 0;
        queue_size = num_verts;
        next_size = 0;
      }

      num_iter = 0;
      while (/*swapped &&*/ num_iter < vert_refine_iter)
      {
        #pragma omp for schedule(guided) reduction(+:num_swapped_2) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v = queue[i];
          in_queue[v] = false;
          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0;

          int part = parts[v];
          unsigned out_degree = out_degree(g, v);
          int* outs = out_vertices(g, v);
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out = outs[j];
            int part_out = parts[out];
            part_counts[part_out]++;
          }

          int max_part = -1;
          int max_count = -1;
          for (int p = 0; p < num_parts; ++p)
            if (part_counts[p] > max_count)
            {
              max_count = part_counts[p];
              max_part = p;
            }

          if (max_part != part)
          {
            double new_max_imb = (double)(part_sizes[max_part] + 1) / avg_size;
            if (new_max_imb < vert_balance)
            {
              ++num_swapped_2;
              parts[v] = max_part;
              #pragma omp atomic
              ++part_sizes[max_part];
              #pragma omp atomic
              --part_sizes[part];

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
					  printf("num_swapped_2 (vert, refine): %d\n", num_swapped_2);
          #endif
          int* temp = queue;
          queue = queue_next;
          queue_next = temp;
          bool* temp_b = in_queue;
          in_queue = in_queue_next;
          in_queue_next = temp_b;
          queue_size = next_size;
          next_size = 0;

          num_swapped_2 = 0;

          max_v = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_sizes[p] / avg_size > max_v)
              max_v = (double)part_sizes[p] / avg_size;
          }
          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        }
      } // end while

      #pragma omp single
      {
        if (max_v > vert_balance*1.01 &&
            t == vert_outer_iter-1 &&
            num_tries < MAX_TRIES)
        {
          --t;
          if (max_v < running_max_v)
          {
            running_max_v = max_v;
            printf("Vertex balance missed, attempting further iterations: (%2.3lf)\n", max_v);
          }
          else
            ++num_tries;
        }
        else
          ++t;
      }
    } // end for

    delete [] part_counts;
    delete [] part_weights;

  } // end par

  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}





/*
'##:::::'##::::'########:::::'###::::'##::::::::::'##::::'##:'########:'########::'########::'######::
 ##:'##: ##:::: ##.... ##:::'## ##::: ##:::::::::: ##:::: ##: ##.....:: ##.... ##:... ##..::'##... ##:
 ##: ##: ##:::: ##:::: ##::'##:. ##:: ##:::::::::: ##:::: ##: ##::::::: ##:::: ##:::: ##:::: ##:::..::
 ##: ##: ##:::: ########::'##:::. ##: ##:::::::::: ##:::: ##: ######::: ########::::: ##::::. ######::
 ##: ##: ##:::: ##.... ##: #########: ##::::::::::. ##:: ##:: ##...:::: ##.. ##:::::: ##:::::..... ##:
 ##: ##: ##:::: ##:::: ##: ##.... ##: ##:::::::::::. ## ##::: ##::::::: ##::. ##::::: ##::::'##::: ##:
. ###. ###::::: ########:: ##:::: ##: ########::::::. ###:::: ########: ##:::. ##:::: ##::::. ######::
:...::...::::::........:::..:::::..::........::::::::...:::::........::..:::::..:::::..::::::......:::
*/

void label_balance_verts_weighted(
  pulp_graph_t& g, int num_parts, int* parts,
  int vert_outer_iter, int vert_balance_iter, int vert_refine_iter,
  double vert_balance)
{
  int   num_verts  = g.n;
  long* part_sizes = new long[num_parts];

  bool has_vwgts = (g.vertex_weights != NULL);
  bool has_ewgts = (g.edge_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;

  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  double avg_size      = (double)g.vertex_weights_sum / (double)num_parts;
  int    num_swapped_1 = 0;
  int    num_swapped_2 = 0;
  double max_v;
  double running_max_v = (double)num_verts;

  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int   queue_size    = num_verts;
  int   next_size     = 0;
  int   t             = 0;
  int   num_tries     = 0;

  #pragma omp parallel
  {
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

    #pragma omp barrier

    double* part_counts  = new double[num_parts];
    double* part_weights = new double[num_parts];

    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    for (int p = 0; p < num_parts; ++p)
    {
      if (part_sizes[p] == 0)
        part_weights[p] = 1; //numeric_limits<double>::max();
        //TODO(julie): max here is not aligned with the idea of bin packing
        //where we want to minimize the number of bins
      else
      {
        part_weights[p] = (vert_balance * avg_size / (double)part_sizes[p]) - 1.0;
        if (part_weights[p] < 0.0) part_weights[p] = 0.0;
      }
    }

    while(t < vert_outer_iter)
    {

      #pragma omp for schedule(static) nowait
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;

      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        in_queue_next[i] = false;

      #pragma omp single
      {
        num_swapped_1 = 0;
        queue_size = num_verts;
        next_size = 0;
      }

      int num_iter = 0;
      while (/*swapped &&*/ num_iter < vert_balance_iter)
      {
        #pragma omp for schedule(guided) reduction(+:num_swapped_1) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v        = queue[i];
          in_queue[v]  = false;
          int part     = parts[v];
          int v_weight = 1;
          if (has_vwgts) v_weight = g.vertex_weights[v];

          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0.0;

          unsigned out_degree = out_degree(g, v);
          int*     outs       = out_vertices(g, v);
          int*     weights    = out_weights(g, v);
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int    out        = outs[j];
            int    part_out   = parts[out];
            double weight_out = 1.0;
            if (has_ewgts) weight_out = (double)weights[j];
            part_counts[part_out] += (double)out_degree(g, out)*weight_out;
          }

          int max_part = part;
          double max_val = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            part_counts[p] *= part_weights[p];

            if (part_counts[p] > max_val)
            {
              max_val = part_counts[p];
              max_part = p;
            }
          }

          if (max_part != part)
          {

            if (g.max_partition_size > 0 &&
                part_sizes[max_part] + v_weight > g.max_partition_size)
                continue;

            parts[v] = max_part;
            ++num_swapped_1;

            #pragma omp atomic
            part_sizes[max_part] += v_weight;
            #pragma omp atomic
            part_sizes[part]     -= v_weight;

            part_weights[part]     = vert_balance * avg_size / (double)part_sizes[part] - 1.0;
            part_weights[max_part] = vert_balance * avg_size / (double)part_sizes[max_part]  - 1.0;

            if (part_weights[part] < 0.0)
              part_weights[part] = 0.0;
            if (part_weights[max_part] < 0.0)
              part_weights[max_part] = 0.0;

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
				    printf("num_swapped_1 (vert, balance): %d\n", num_swapped_1);
          #endif

          int*  temp    = queue;
          queue         = queue_next;
          queue_next    = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          num_swapped_1 = 0;

          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        }
      } // end while

      #pragma omp for schedule(static)
        for (int i = 0; i < num_verts; ++i)
          queue[i] = i;

      #pragma omp single
      {
        num_swapped_2 = 0;
        queue_size    = num_verts;
        next_size     = 0;
      }

      num_iter = 0;
      while (/*swapped &&*/ num_iter < vert_refine_iter)
      {
        #pragma omp for schedule(guided) reduction(+:num_swapped_2) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v        = queue[i];
          in_queue[v]  = false;
          int part     = parts[v];
          int v_weight = 1;
          if (has_vwgts) v_weight = g.vertex_weights[v];

          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0;

          unsigned out_degree = out_degree(g, v);
          int*     outs       = out_vertices(g, v);
          int*     weights    = out_weights(g, v);
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out        = outs[j];
            int part_out   = parts[out];
            int out_weight = 1;
            if (has_ewgts) out_weight = weights[j];
            part_counts[part_out] += out_weight;
          }

          int max_part  = -1;
          int max_count = -1;
          for (int p = 0; p < num_parts; ++p)
            if (part_counts[p] > max_count)
            {
              max_count = part_counts[p];
              max_part  = p;
            }

          if (max_part != part)
          {

            if (g.max_partition_size > 0 &&
                part_sizes[max_part] + v_weight > g.max_partition_size)
              continue;

            double new_max_imb = (double)(part_sizes[max_part] + v_weight) / avg_size;
            if (new_max_imb < vert_balance)
            {
              ++num_swapped_2;
              parts[v] = max_part;
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
					  printf("num_swapped_2 (vert, refine): %d\n", num_swapped_2);
          #endif
          int*  temp    = queue;
          queue         = queue_next;
          queue_next    = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          num_swapped_2 = 0;

          max_v = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_sizes[p] / avg_size > max_v)
              max_v = (double)part_sizes[p] / avg_size;
          }
          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        }
      } // end while

      #pragma omp single
      {
        if (max_v > vert_balance*1.01 &&
            t == vert_outer_iter-1 &&
            num_tries < MAX_TRIES)
        {
          --t;
          if (max_v < running_max_v)
          {
            running_max_v = max_v;
            printf("Vertex balance missed, attempting further iterations: (%2.3lf)\n", max_v);
          }
          else
            ++num_tries;
        }
        else
          ++t;
      }
    } // end while t < vert_outer_iter

    delete [] part_counts;
    delete [] part_weights;

  } // end par


  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}





/*
##      ##    ####    ########     ###    ##          ##     ## ######## ########  ########
##  ##  ##     ##     ##     ##   ## ##   ##          ##     ## ##       ##     ##    ##
##  ##  ##     ##     ##     ##  ##   ##  ##          ##     ## ##       ##     ##    ##
##  ##  ##     ##     ########  ##     ## ##          ##     ## ######   ########     ##
##  ##  ##     ##     ##     ## ######### ##           ##   ##  ##       ##   ##      ##
##  ##  ##     ##     ##     ## ##     ## ##            ## ##   ##       ##    ##     ##
 ###  ###     ####    ########  ##     ## ########       ###    ######## ##     ##    ##
*/

void
label_balance_verts_weighted_interpart(pulp_graph_t& g, int num_parts, int* parts,
  int vert_outer_iter, int vert_balance_iter, int vert_refine_iter, double vert_balance)
{
  bool has_ipwgts = (g.interpartition_weights != NULL);
  if (!has_ipwgts)
  {
    printf("Error: interpartition weights required for weighted interpart label prop\n");
    exit(0);
  }


  int  num_verts = g.n;
  bool has_ewgts = (g.edge_weights != NULL);
  bool has_vwgts = (g.vertex_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;

  long*  part_sizes = new long[num_parts];
  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  double avg_size = (double)g.vertex_weights_sum / (double)num_parts;

  int    num_swapped_1 = 0;
  int    num_swapped_2 = 0;
  double max_v         = 0.0;
  double running_max_v = (double)num_verts;

  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int   queue_size    = num_verts;
  int   next_size     = 0;
  int   t             = 0;
  int   num_tries     = 0;

  #pragma omp parallel
  {
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

    #pragma omp barrier

    double* part_counts  = new double[num_parts];
    double* part_weights = new double[num_parts];

    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;


    for (int p = 0; p < num_parts; ++p)
    {
      if (part_sizes[p] == 0)
        part_weights[p] = 1; // numeric_limits<double>::max();
      else if (g.max_partition_size > 0 &&
               part_sizes[p] > g.max_partition_size)
      {
        part_weights[p] = 0.0;
        printf("Partition %d is overloaded (max_partition_size) at %ld\n", p, part_sizes[p]);
      }
      else
      {
        part_weights[p] = (vert_balance * avg_size / (double)part_sizes[p]) - 1.0;
                if (part_weights[p] < 0.0)
                  part_weights[p] = 0.0;
      }
    }

    // ======================================================
    // Outer loop
    // ======================================================
    while(t < vert_outer_iter)
    {

      #pragma omp for schedule(static) nowait
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;
      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        in_queue_next[i] = false;

      #pragma omp single
      {
        num_swapped_1 = 0;
        queue_size    = num_verts;
        next_size     = 0;
      }

      // ======================================================
      // ######
      // #     #   ##   #        ##   #    #  ####  ######
      // #     #  #  #  #       #  #  ##   # #    # #
      // ######  #    # #      #    # # #  # #      #####
      // #     # ###### #      ###### #  # # #      #
      // #     # #    # #      #    # #   ## #    # #
      // ######  #    # ###### #    # #    #  ####  ######
      // ======================================================
      int num_iter = 0;
      while (/*swapped &&*/ num_iter < vert_balance_iter)
      {
        #pragma omp for schedule(guided) reduction(+:num_swapped_1) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v        = queue[i];
          in_queue[v]  = false;
          int part     = parts[v];
          int v_weight = 1;
          if (has_vwgts) v_weight = g.vertex_weights[v];

          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0.0;

          unsigned out_degree = out_degree(g, v);
          int*     outs       = out_vertices(g, v);
          int*     weights    = out_weights(g, v);
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int    out        = outs[j];
            int    part_out   = parts[out];  // partition of neighbor vertex
            double weight_out = 1.0;
            if (has_ewgts)
              weight_out = (double)weights[j];
            part_counts[part_out] += (double)out_degree(g, out) * weight_out;
          }

          // -----------------------------------------------------
          // Inter-partition weight adjustment
          // -----------------------------------------------------
          int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
          for (int p = 0; p < num_parts; ++p)
              part_counts[p] *= partition_comm_weights[p]; // adjust part_counts based on partition_comm_weights for each partition pairs

          // -----------------------------------------------------
          // Find the partition with the maximum count
          // -----------------------------------------------------
          int max_part = part;
          double max_val = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            // Adjusts the count to reflect the desirability of adding more vertices to this
            // partition.
            //
            // - If a partition is underloaded (i.e., its size is less than the average
            // size), its weight will be greater than 0, making it more attractive.
            //
            // - Conversely, if a partition is overloaded, its weight will be 0, making it
            // less attractive.
            //
            // Biases the decision towards balancing the partitions according to the
            // specified balance constraint.

            part_counts[p] *= part_weights[p]; // adjust part_counts based on part_weights

            if (part_counts[p] > max_val)
            {
              max_val  = part_counts[p];
              max_part = p;
            }
          }

          // -----------------------------------------------------
          // Swap vertex v to the partition with the maximum count
          // -----------------------------------------------------
          if (max_part != part)
          {
            // Unless in bin packing mode, do not move a vertex if it would empty its partition
            if (!g.do_bin_packing && (part_sizes[part] - v_weight <= 0))
              continue;

            if (g.max_partition_size > 0 &&
                part_sizes[max_part] + v_weight > g.max_partition_size)
              continue;

            parts[v] = max_part; // reassign vertex v to the largest partition
            ++num_swapped_1;     // increment the number of vertices swapped

            #pragma omp atomic
            part_sizes[max_part] += v_weight;
            #pragma omp atomic
            part_sizes[part]     -= v_weight;

            part_weights[part]     = vert_balance * avg_size / (double)part_sizes[part] - 1.0;
            part_weights[max_part] = vert_balance * avg_size / (double)part_sizes[max_part] - 1.0;

            if (part_weights[part] < 0.0)
              part_weights[part] = 0.0;
            if (part_weights[max_part] < 0.0)
              part_weights[max_part] = 0.0;

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
				    printf("num_swapped_1 (vert, balance): %d\n", num_swapped_1);
          #endif

          int*  temp    = queue;
          queue         = queue_next;
          queue_next    = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          num_swapped_1 = 0;

          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        } // end single
      } // end while


      // ================================================================
      // ######
      // #     # ###### ###### # #    # ###### #    # ###### #    # #####
      // #     # #      #      # ##   # #      ##  ## #      ##   #   #
      // ######  #####  #####  # # #  # #####  # ## # #####  # #  #   #
      // #   #   #      #      # #  # # #      #    # #      #  # #   #
      // #    #  #      #      # #   ## #      #    # #      #   ##   #
      // #     # ###### #      # #    # ###### #    # ###### #    #   #
      // ================================================================
      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;

      #pragma omp single
      {
        num_swapped_2 = 0;
        queue_size    = num_verts;
        next_size     = 0;
      }

      num_iter = 0;
      while (/*swapped &&*/ num_iter < vert_refine_iter)
      {
        #pragma omp for schedule(guided) reduction(+:num_swapped_2) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v        = queue[i];
          int part     = parts[v];
          in_queue[v]  = false;

          int v_weight = 1;
          if (has_vwgts) v_weight = g.vertex_weights[v];

          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0;

          unsigned out_degree = out_degree(g, v);
          int*     outs       = out_vertices(g, v);
          int*     weights    = out_weights(g, v);
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out        = outs[j];
            int part_out   = parts[out];
            int out_weight = 1;
            if (has_ewgts) out_weight = weights[j];
            part_counts[part_out] += out_weight;
          }

          // -----------------------------------------------------
          // Inter-partition weight adjustment
          // -----------------------------------------------------
          int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
          for (int p = 0; p < num_parts; ++p)
            part_counts[p] *= partition_comm_weights[p];

          // -----------------------------------------------------
          // Find the partition with the maximum count
          // -----------------------------------------------------
          int max_part  = -1;
          int max_count = -1;
          for (int p = 0; p < num_parts; ++p)
            if (part_counts[p] > max_count)
            {
              max_count = part_counts[p];
              max_part  = p;
            }

          // -----------------------------------------------------
          // Swap vertex v to the partition with the maximum count
          // -----------------------------------------------------
          if (max_part != part)
          {

            // Unless in bin packing mode, do not move a vertex if it would empty its partition
            if (!g.do_bin_packing && (part_sizes[part] - v_weight <= 0))
              continue;

            if (g.max_partition_size > 0 &&
                (part_sizes[max_part]+v_weight) > g.max_partition_size)
              continue;

            double new_max_imb = (double)(part_sizes[max_part] + v_weight) / avg_size;

            if (new_max_imb < vert_balance)
            {
              ++num_swapped_2;
              parts[v] = max_part;

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
            } // end if new_max_imb
          } // end if max_part
        } // end for

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
					  printf("num_swapped_2 (vert, refine): %d\n", num_swapped_2);
          #endif
          int*  temp    = queue;
          queue         = queue_next;
          queue_next    = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          num_swapped_2 = 0;

          max_v = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_sizes[p] / avg_size > max_v)
              max_v = (double)part_sizes[p] / avg_size;
          }
          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        } // end single
      } // end while

      #pragma omp single
      {
        if (max_v > vert_balance*1.01 &&
            t == vert_outer_iter-1 &&
            num_tries < MAX_TRIES)
        {
          --t;
          if (max_v < running_max_v)
          {
            running_max_v = max_v;
            printf("Vertex balance missed, attempting further iterations: (%2.3lf)\n", max_v);
          }
          else
            ++num_tries;
        }
        else
          ++t;
      } // end single
    } // end for

    delete [] part_counts;
    delete [] part_weights;
  } // end omp parallel

  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}





/*
##      ##         ######     ###    ########         ####        ##     ## ######## ########  ########        ########     ###    ##
##  ##  ##        ##    ##   ## ##   ##     ##         ##         ##     ## ##       ##     ##    ##           ##     ##   ## ##   ##
##  ##  ##        ##        ##   ##  ##     ##         ##         ##     ## ##       ##     ##    ##           ##     ##  ##   ##  ##
##  ##  ##        ##       ##     ## ########          ##         ##     ## ######   ########     ##           ########  ##     ## ##
##  ##  ##        ##       ######### ##                ##          ##   ##  ##       ##   ##      ##           ##     ## ######### ##
##  ##  ## ###    ##    ## ##     ## ##        ###     ##  ###      ## ##   ##       ##    ##     ##    ###    ##     ## ##     ## ##       ###
 ###  ###  ###     ######  ##     ## ##        ###    #### ###       ###    ######## ##     ##    ##    ###    ########  ##     ## ######## ###
*/

void
label_balance_verts_weighted_interpart_capacity(pulp_graph_t& g, int num_parts, int* parts,
  int vert_outer_iter, int vert_balance_iter, int vert_refine_iter, double vert_balance)
{
  printf("running label_balance_verts_weighted_interpart_capacity\n");

  bool has_ipwgts        = (g.interpartition_weights != NULL);
  bool has_p_capacities  = (g.partition_capacities != NULL);
  if (!has_p_capacities || !has_ipwgts)
  {
    printf("Error: partition capacities required for weighted interpart label prop\n");
    exit(0);
  }

  int  num_verts = g.n;
  bool has_ewgts = (g.edge_weights != NULL);
  bool has_vwgts = (g.vertex_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;

  long*  part_sizes = new long[num_parts];
  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;

  double* avg_sizes     = new double[num_parts];
  for (int i = 0; i < num_parts; ++i)
  {
    if (g.max_partition_size > 0) // if max_partition_size is set
      avg_sizes[i] = (double)g.max_partition_size * (double)g.partition_capacities[i];
    else
    {
      double unit_avg_size = (double)g.vertex_weights_sum / (double)g.partition_capacities_sum;
      avg_sizes[i] = unit_avg_size * g.partition_capacities[i];
    }
    printf("avg_sizes[%d] = %lf\n", i, avg_sizes[i]);
  }

  int    num_swapped_1 = 0;
  int    num_swapped_2 = 0;
  double max_v         = 0.0;
  double running_max_v = (double)num_verts;

  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int   queue_size    = num_verts;
  int   next_size     = 0;
  int   t             = 0;
  int   num_tries     = 0;

  #pragma omp parallel
  {
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

    #pragma omp barrier

    double* part_counts  = new double[num_parts];
    double* part_weights = new double[num_parts];

    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;

    printf("%10s %10s %10s\n", "part_sizes", "part_weights", "avg_sizes");
    for (int p = 0; p < num_parts; ++p)
    {
      if (part_sizes[p] == 0)
        part_weights[p] = 1; //numeric_limits<double>::max();
      else
      {
        part_weights[p] = (vert_balance * avg_sizes[p] / (double)part_sizes[p]) - 1.0;
        if (part_weights[p] < 0.0) part_weights[p] = 0.0;
      }
      printf("%10ld %10lf %10lf\n", part_sizes[p], part_weights[p], avg_sizes[p]);
    }

    // ======================================================
    // Outer loop
    // ======================================================
    while(t < vert_outer_iter)
    {

      #pragma omp for schedule(static) nowait
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;
      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        in_queue_next[i] = false;

      #pragma omp single
      {
        num_swapped_1 = 0;
        queue_size    = num_verts;
        next_size     = 0;
      }

      // ======================================================
      // ######
      // #     #   ##   #        ##   #    #  ####  ######
      // #     #  #  #  #       #  #  ##   # #    # #
      // ######  #    # #      #    # # #  # #      #####
      // #     # ###### #      ###### #  # # #      #
      // #     # #    # #      #    # #   ## #    # #
      // ######  #    # ###### #    # #    #  ####  ######
      // ======================================================
      int num_iter = 0;
      while (/*swapped &&*/ num_iter < vert_balance_iter)
      {
        #pragma omp for schedule(guided) reduction(+:num_swapped_1) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v        = queue[i];
          in_queue[v]  = false;
          int part     = parts[v];
          int v_weight = 1;
          if (has_vwgts) v_weight = g.vertex_weights[v];

          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0.0;

          unsigned out_degree = out_degree(g, v);
          int*     outs       = out_vertices(g, v);
          int*     weights    = out_weights(g, v);
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int    out        = outs[j];
            int    part_out   = parts[out];  // partition of neighbor vertex
            double weight_out = 1.0;
            if (has_ewgts)
              weight_out = (double)weights[j];
            part_counts[part_out] += (double)out_degree(g, out) * weight_out;
          }

          // -----------------------------------------------------
          // Inter-partition weight adjustment
          // -----------------------------------------------------
          int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
          for (int p = 0; p < num_parts; ++p)
              part_counts[p] *= partition_comm_weights[p]; // adjust part_counts based on partition_comm_weights for each partition pairs

          // -----------------------------------------------------
          // Find the partition with the maximum count
          // -----------------------------------------------------
          int max_part = part;
          double max_val = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            // Adjusts the count to reflect the desirability of adding more vertices to this
            // partition.
            //
            // - If a partition is underloaded (i.e., its size is less than the average
            //   size), its weight will be greater than 0, making it more attractive.
            //
            // - Conversely, if a partition is overloaded, its weight will be 0, making it
            //   impossible to add more vertices.
            //
            // Biases the decision towards balancing the partitions according to the
            // specified balance constraint.

            part_counts[p] *= part_weights[p]; // adjust part_counts based on part_weights

            if (part_counts[p] > max_val)
            {
              max_val  = part_counts[p];
              max_part = p;
            }
          }

          // -----------------------------------------------------
          // Swap vertex v to the partition with the maximum count
          // -----------------------------------------------------
          if (max_part != part)
          {
            // Unless in bin packing mode, do not move a vertex if it would empty its partition
            if (!g.do_bin_packing && (part_sizes[part] - v_weight <= 0))
            {
              continue;
            }

            if (g.max_partition_size > 0 &&
                (part_sizes[max_part]+v_weight) > (g.max_partition_size*g.partition_capacities[max_part]))
            {
              int thread_num = omp_get_thread_num();
              printf("Thread %d: WARNING: Max partition size exceeded on part %d: %ld + %d > %d * %d\n",
                    thread_num, max_part, part_sizes[max_part], v_weight, g.max_partition_size, g.partition_capacities[max_part]);
              continue;
            }

            parts[v] = max_part; // reassign vertex v to the largest partition
            ++num_swapped_1;     // increment the number of vertices swapped

            #pragma omp atomic
            part_sizes[max_part] += v_weight;
            #pragma omp atomic
            part_sizes[part]     -= v_weight;

            part_weights[part]     = (vert_balance * avg_sizes[part] / (double)part_sizes[part]) - 1.0;
            part_weights[max_part] = (vert_balance * avg_sizes[part] / (double)part_sizes[max_part]) - 1.0;

            if (part_weights[part] < 0.0)
              part_weights[part] = 0.0;
            if (part_weights[max_part] < 0.0)
              part_weights[max_part] = 0.0;

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
				    printf("num_swapped_1 (vert, balance): %d\n", num_swapped_1);
          #endif

          int*  temp    = queue;
          queue         = queue_next;
          queue_next    = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          num_swapped_1 = 0;

          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        } // end single
      } // end while


      // ================================================================
      // ######
      // #     # ###### ###### # #    # ###### #    # ###### #    # #####
      // #     # #      #      # ##   # #      ##  ## #      ##   #   #
      // ######  #####  #####  # # #  # #####  # ## # #####  # #  #   #
      // #   #   #      #      # #  # # #      #    # #      #  # #   #
      // #    #  #      #      # #   ## #      #    # #      #   ##   #
      // #     # ###### #      # #    # ###### #    # ###### #    #   #
      // ================================================================
      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;

      #pragma omp single
      {
        num_swapped_2 = 0;
        queue_size    = num_verts;
        next_size     = 0;
      }

      num_iter = 0;
      while (/*swapped &&*/ num_iter < vert_refine_iter)
      {
        #pragma omp for schedule(guided) reduction(+:num_swapped_2) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v        = queue[i];
          int part     = parts[v];
          in_queue[v]  = false;

          int v_weight = 1;
          if (has_vwgts) v_weight = g.vertex_weights[v];

          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0;

          unsigned out_degree = out_degree(g, v);
          int*     outs       = out_vertices(g, v);
          int*     weights    = out_weights(g, v);
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out        = outs[j];
            int part_out   = parts[out];
            int out_weight = 1;
            if (has_ewgts) out_weight = weights[j];
            part_counts[part_out] += out_weight;
          }

          // -----------------------------------------------------
          // Inter-partition weight adjustment
          // -----------------------------------------------------
          int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
          for (int p = 0; p < num_parts; ++p)
            part_counts[p] *= partition_comm_weights[p];

          // -----------------------------------------------------
          // Find the partition with the maximum count
          // -----------------------------------------------------
          int max_part  = -1;
          int max_count = -1;
          for (int p = 0; p < num_parts; ++p)
            if (part_counts[p] > max_count)
            {
              max_count = part_counts[p];
              max_part  = p;
            }

          // -----------------------------------------------------
          // Swap vertex v to the partition with the maximum count
          // -----------------------------------------------------
          if (max_part != part)
          {
            if (!g.do_bin_packing &&
                (part_sizes[part] - v_weight <= 0))
            {
              continue;
            }

            if (g.max_partition_size > 0 &&
                (part_sizes[max_part]+v_weight) > (g.max_partition_size*g.partition_capacities[max_part]))
            {
              int thread_num = omp_get_thread_num();
              printf("Thread %d: WARNING: Max partition size exceeded on part %d: %ld + %d > %d * %d\n",
                    thread_num, max_part, part_sizes[max_part], v_weight, g.max_partition_size, g.partition_capacities[max_part]);
              continue;
            }

            double new_max_imb = (double)(part_sizes[max_part] + v_weight)
                                         / avg_sizes[max_part];

            if (new_max_imb < vert_balance)
            {
              ++num_swapped_2;
              parts[v] = max_part;

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
            } // end if new_max_imb
          } // end if max_part
        } // end for

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
					  printf("num_swapped_2 (vert, refine): %d\n", num_swapped_2);
          #endif
          int*  temp    = queue;
          queue         = queue_next;
          queue_next    = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          num_swapped_2 = 0;

          max_v = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_sizes[p] / avg_sizes[p] > max_v)
              max_v = (double)part_sizes[p] / avg_sizes[p];
          }
          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        } // end single
      } // end while

      #pragma omp single
      {
        if (max_v > vert_balance*1.01 &&
            t == vert_outer_iter-1 &&
            num_tries < MAX_TRIES)
        {
          --t;
          if (max_v < running_max_v)
          {
            running_max_v = max_v;
            printf("Vertex balance missed, attempting further iterations: (%2.3lf)\n", max_v);
          }
          else
            ++num_tries;
        }
        else
          ++t;
      } // end omp single
    } // end while t < vert_outer_iter

    delete [] part_counts;
    delete [] part_weights;

  } // end omp parallel

  delete [] avg_sizes;

  delete [] part_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}



