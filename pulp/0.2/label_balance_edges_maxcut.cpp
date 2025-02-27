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
'##::::'##::::'###::::'##::::'##::'######::'##::::'##:'########::::'########:'########:::'######:::'########::::'########:::::'###::::'##::::::::::'###::::'##::: ##::'######::'########:
 ###::'###:::'## ##:::. ##::'##::'##... ##: ##:::: ##:... ##..::::: ##.....:: ##.... ##:'##... ##:: ##.....::::: ##.... ##:::'## ##::: ##:::::::::'## ##::: ###:: ##:'##... ##: ##.....::
 ####'####::'##:. ##:::. ##'##::: ##:::..:: ##:::: ##:::: ##::::::: ##::::::: ##:::: ##: ##:::..::: ##:::::::::: ##:::: ##::'##:. ##:: ##::::::::'##:. ##:: ####: ##: ##:::..:: ##:::::::
 ## ### ##:'##:::. ##:::. ###:::: ##::::::: ##:::: ##:::: ##::::::: ######::: ##:::: ##: ##::'####: ######:::::: ########::'##:::. ##: ##:::::::'##:::. ##: ## ## ##: ##::::::: ######:::
 ##. #: ##: #########::: ## ##::: ##::::::: ##:::: ##:::: ##::::::: ##...:::: ##:::: ##: ##::: ##:: ##...::::::: ##.... ##: #########: ##::::::: #########: ##. ####: ##::::::: ##...::::
 ##:.:: ##: ##.... ##:: ##:. ##:: ##::: ##: ##:::: ##:::: ##::::::: ##::::::: ##:::: ##: ##::: ##:: ##:::::::::: ##:::: ##: ##.... ##: ##::::::: ##.... ##: ##:. ###: ##::: ##: ##:::::::
 ##:::: ##: ##:::: ##: ##:::. ##:. ######::. #######::::: ##::::::: ########: ########::. ######::: ########:::: ########:: ##:::: ##: ########: ##:::: ##: ##::. ##:. ######:: ########:
..:::::..::..:::::..::..:::::..:::......::::.......::::::..::::::::........::........::::......::::........:::::........:::..:::::..::........::..:::::..::..::::..:::......:::........::*/

void label_balance_edges_maxcut(pulp_graph_t& g, int num_parts, int* parts,
  int edge_outer_iter, int edge_balance_iter, int edge_refine_iter,
  double vert_balance, double edge_balance)
{
  int       num_verts       = g.n;
  unsigned  num_edges       = g.m;
  unsigned  cut_size        = 0;
  int*      part_sizes      = new int[num_parts];
  unsigned* part_edge_sizes = new unsigned[num_parts];
  int*      part_cut_sizes  = new int[num_parts];

  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;
  for (int i = 0; i < num_parts; ++i)
    part_edge_sizes[i] = 0;
  for (int i = 0; i < num_parts; ++i)
    part_cut_sizes[i] = 0;

  double avg_size      = num_verts / num_parts;
  double avg_edge_size = num_edges / num_parts;
  int    num_swapped_1 = 0;
  int    num_swapped_2 = 0;
  double max_e         = 0.0;
  double max_c         = 0.0;
  double running_max_e = (double)num_edges;
  double weight_exponent_e = 1.0;
  double weight_exponent_c = 1.0;

  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int queue_size;
  int next_size;
  int t         = 0;
  int num_tries = 0;

  #pragma omp parallel
  {
    int*      part_sizes_thread      = new int[num_parts];
    unsigned* part_edge_sizes_thread = new unsigned[num_parts];
    int*      part_cut_sizes_thread  = new int[num_parts];
    unsigned  cut_size_thread = 0;
    for (int i = 0; i < num_parts; ++i)
      part_sizes_thread[i] = 0;
    for (int i = 0; i < num_parts; ++i)
      part_edge_sizes_thread[i] = 0;
    for (int i = 0; i < num_parts; ++i)
      part_cut_sizes_thread[i] = 0;

    #pragma omp for schedule(guided) nowait
    for (int i = 0; i < num_verts; ++i)
    {
      int      part       = parts[i];
      unsigned out_degree = out_degree(g, i);
      ++part_sizes_thread[part];
      part_edge_sizes_thread[part] += out_degree;

      int* outs = out_vertices(g, i);
      for (unsigned j = 0; j < out_degree; ++j)
      {
        int out      = outs[j];
        int out_part = parts[out];
        if (out_part != part)
        {
          ++part_cut_sizes_thread[part];
          ++cut_size_thread;
        }
      }
    }

    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_sizes[i] += part_sizes_thread[i];
    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_edge_sizes[i] += part_edge_sizes_thread[i];
    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_cut_sizes[i] += part_cut_sizes_thread[i];
    #pragma omp atomic
    cut_size += cut_size_thread;

    delete [] part_sizes_thread;
    delete [] part_edge_sizes_thread;
    delete [] part_cut_sizes_thread;


    double  avg_cut_size      = (double)cut_size / (double)num_parts;
    double* part_counts       = new double[num_parts];
    double* part_weights      = new double[num_parts];
    double* part_edge_weights = new double[num_parts];
    double* part_cut_weights  = new double[num_parts];

    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;


    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      in_queue_next[i] = false;

    while(t < edge_outer_iter)
    {

      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;

      #pragma omp single
      {
        num_swapped_1 = 0;
        queue_size    = num_verts;
        next_size     = 0;

        max_e = 0.0;
        max_c = 0.0;
        for (int p = 0; p < num_parts; ++p)
        {
          if ((double)part_edge_sizes[p] / avg_edge_size > max_e)
            max_e = (double)part_edge_sizes[p] / avg_edge_size;
          if ((double)part_cut_sizes[p] / avg_cut_size > max_c)
            max_c = (double)part_cut_sizes[p] / avg_cut_size;
        }
        if (max_e < edge_balance)
        {
          max_e              = edge_balance;
          weight_exponent_e  = 1.0;
          weight_exponent_c *= max_c;
        }
        else
        {
          weight_exponent_e *= max_e / edge_balance;
          weight_exponent_c  = 1.0;
        }
      }

      int num_iter = 0;
      while (/*swapped &&*/ num_iter < edge_balance_iter)
      {
        for (int p = 0; p < num_parts; ++p)
        {
          part_weights[p]      = vert_balance * avg_size / (double)part_sizes[p] - 1.0;
          part_edge_weights[p] = max_e * avg_edge_size / (double)part_edge_sizes[p] - 1.0;
          part_cut_weights[p]  = max_c * avg_cut_size / (double)part_cut_sizes[p] - 1.0;
          if (part_weights[p] < 0.0)
            part_weights[p] = 0.0;
          if (part_edge_weights[p] < 0.0)
            part_edge_weights[p] = 0.0;
          if (part_cut_weights[p] < 0.0)
            part_cut_weights[p] = 0.0;
        }

        #pragma omp for schedule(guided) reduction(+:num_swapped_1) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v = queue[i];
          in_queue[v] = false;
          int part = parts[v];
          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0.0;

          unsigned out_degree = out_degree(g, v);
          int*     outs       = out_vertices(g, v);
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out      = outs[j];
            int part_out = parts[out];
            part_counts[part_out] += 1.0;
          }

          int    max_part   = part;
          double max_val    = 0.0;
          int    part_count = (int)part_counts[part];
          int    max_count  = 0;
          for (int p = 0; p < num_parts; ++p)
          {
            int count_init = (int)part_counts[p];
            if (part_weights[p] > 0.0 && part_edge_weights[p] > 0.0 && part_cut_weights[p] > 0.0)
              part_counts[p] *= (part_edge_weights[p]*weight_exponent_e * part_cut_weights[p]*weight_exponent_c);
            else
              part_counts[p] = 0.0;

            if (part_counts[p] > max_val)
            {
              max_val   = part_counts[p];
              max_count = count_init;
              max_part  = p;
            }
          }

          if (max_part != part)
          {
            parts[v] = max_part;
            ++num_swapped_1;
            int diff_part     = 2*part_count - out_degree;
            int diff_max_part = out_degree - 2*max_count;
            int diff_cut      = diff_part + diff_max_part;
            #pragma omp atomic
            cut_size += diff_cut;
            #pragma omp atomic
            part_cut_sizes[part] += diff_part;
            #pragma omp atomic
            part_cut_sizes[max_part] += diff_max_part;
            #pragma omp atomic
            --part_sizes[part];
            #pragma omp atomic
            ++part_sizes[max_part];
            #pragma omp atomic
            part_edge_sizes[part] -= out_degree;
            #pragma omp atomic
            part_edge_sizes[max_part] += out_degree;

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

            avg_cut_size            = cut_size / num_parts;
            part_weights[part]      = vert_balance * avg_size / (double)part_sizes[part] - 1.0;
            part_edge_weights[part] = max_e * avg_edge_size / (double)part_edge_sizes[part] - 1.0;
            part_cut_weights[part]  = max_c * avg_cut_size / (double)part_cut_sizes[part] - 1.0;

            part_weights[max_part]      = vert_balance * avg_size / (double)part_sizes[max_part]  - 1.0;
            part_edge_weights[max_part] = max_e * avg_edge_size / (double)part_edge_sizes[max_part] - 1.0;
            part_cut_weights[max_part]  = max_c * avg_cut_size / (double)part_cut_sizes[max_part] - 1.0;

            if (part_weights[part] < 0.0)
              part_weights[part] = 0.0;
            if (part_edge_weights[part] < 0.0)
              part_edge_weights[part] = 0.0;
            if (part_cut_weights[part] < 0.0)
              part_cut_weights[part] = 0.0;

            if (part_weights[max_part] < 0.0)
              part_weights[max_part] = 0.0;
            if (part_edge_weights[max_part] < 0.0)
              part_edge_weights[max_part] = 0.0;
            if (part_cut_weights[max_part] < 0.0)
              part_cut_weights[max_part] = 0.0;
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
            printf("%d -- V: %2.2lf   E: %2.2lf, %lf   C: %2.2lf, %lf\n", num_swapped_1, vert_balance, max_e, weight_exponent_e, max_c, weight_exponent_c);
          #endif
          int* temp = queue;
          queue      = queue_next;
          queue_next = temp;
          bool* temp_b = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          max_e = 0.0;
          max_c = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_edge_sizes[p] / avg_edge_size > max_e)
              max_e = (double)part_edge_sizes[p] / avg_edge_size;
            if ((double)part_cut_sizes[p] / avg_cut_size > max_c)
              max_c = (double)part_cut_sizes[p] / avg_cut_size;
          }
          if (max_e < edge_balance)
          {
            max_e              = edge_balance;
            weight_exponent_e  = 1.0;
            weight_exponent_c *= max_c;
          }
          else
          {
            weight_exponent_e *= max_e / edge_balance;
            weight_exponent_c  = 1.0;
          }

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
      while (/*swapped &&*/ num_iter < edge_refine_iter)
      {
        for (int p = 0; p < num_parts; ++p)
        {
          part_weights[p]      = vert_balance * avg_size / (double)part_sizes[p] - 1.0;
          part_edge_weights[p] = max_e * avg_edge_size / (double)part_edge_sizes[p] - 1.0;
          part_cut_weights[p]  = max_c * avg_cut_size / (double)part_cut_sizes[p] - 1.0;
          if (part_weights[p] < 0.0)
            part_weights[p] = 0.0;
          if (part_edge_weights[p] < 0.0)
            part_edge_weights[p] = 0.0;
          if (part_cut_weights[p] < 0.0)
            part_cut_weights[p] = 0.0;
        }

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
          int part_count = part_counts[part];
          for (int p = 0; p < num_parts; ++p)
            if (part_counts[p] > max_count)
            {
              max_count = part_counts[p];
              max_part = p;
            }

          if (max_part != part)
          {
            double new_max_imb      = (double)(part_sizes[max_part] + 1) / avg_size;
            double new_max_edge_imb = (double)(part_edge_sizes[max_part] + out_degree) / avg_edge_size;
            double new_max_cut_imb  = (double)(part_cut_sizes[max_part] + out_degree - 2*max_count) / avg_cut_size;
            double new_cut_imb      = (double)(part_cut_sizes[part] + 2*part_count - out_degree) / avg_cut_size;
            if ( new_max_imb < vert_balance &&
              new_max_edge_imb < max_e &&
              new_max_cut_imb < max_c && new_cut_imb < max_c)
            {
              ++num_swapped_2;
              parts[v] = max_part;
              int diff_part     = 2*part_count - out_degree;
              int diff_max_part = out_degree - 2*max_count;
              int diff_cut      = diff_part + diff_max_part;
              #pragma omp atomic
              cut_size += diff_cut;
              #pragma omp atomic
              part_cut_sizes[part] += diff_part;
              #pragma omp atomic
              part_cut_sizes[max_part] += diff_max_part;
              #pragma omp atomic
              ++part_sizes[max_part];
              #pragma omp atomic
              --part_sizes[part];
              #pragma omp atomic
              part_edge_sizes[max_part] += out_degree;
              #pragma omp atomic
              part_edge_sizes[part] -= out_degree;

              avg_cut_size = cut_size / num_parts;

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
            printf("%d -- V: %2.2lf  E: %2.2lf  C: %2.2lf\n", num_swapped_2, vert_balance, max_e, max_c);
          #endif
          int* temp = queue;
          queue      = queue_next;
          queue_next = temp;
          bool* temp_b = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          num_swapped_2 = 0;

          max_e = 0.0;
          max_c = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_edge_sizes[p] / avg_edge_size > max_e)
              max_e = (double)part_edge_sizes[p] / avg_edge_size;
            if ((double)part_cut_sizes[p] / avg_cut_size > max_c)
              max_c = (double)part_cut_sizes[p] / avg_cut_size;
          }
          if (max_e < edge_balance)
          {
            max_e              = edge_balance;
            weight_exponent_e  = 1.0;
            weight_exponent_c *= max_c;
          }
          else
          {
            weight_exponent_e *= max_e / edge_balance;
            weight_exponent_c = 1.0;
          }

        #if OUTPUT_STEP
          evaluate_quality(g, num_parts, parts);
        #endif
        }
      }

      #pragma omp single
      {
        if (max_e > edge_balance*1.01 &&
            t == edge_outer_iter-1 &&
            num_tries < MAX_TRIES)
        {
          --t;
          if (max_e < running_max_e*0.99)
          {
            running_max_e = max_e;
            printf("Edge balance missed, attempting further iterations: (%2.3lf)\n", max_e);
          }
          else
            ++num_tries;
        }
        else
          ++t;
      }

    } // end while

    delete [] part_counts;
    delete [] part_weights;
    delete [] part_edge_weights;
    delete [] part_cut_weights;

  } // end par

  delete [] part_sizes;
  delete [] part_edge_sizes;
  delete [] part_cut_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}






/*
'##:::::'##::::'########:'########:::'######:::'########::::'########:::::'###::::'##:::::::
 ##:'##: ##:::: ##.....:: ##.... ##:'##... ##:: ##.....::::: ##.... ##:::'## ##::: ##:::::::
 ##: ##: ##:::: ##::::::: ##:::: ##: ##:::..::: ##:::::::::: ##:::: ##::'##:. ##:: ##:::::::
 ##: ##: ##:::: ######::: ##:::: ##: ##::'####: ######:::::: ########::'##:::. ##: ##:::::::
 ##: ##: ##:::: ##...:::: ##:::: ##: ##::: ##:: ##...::::::: ##.... ##: #########: ##:::::::
 ##: ##: ##:::: ##::::::: ##:::: ##: ##::: ##:: ##:::::::::: ##:::: ##: ##.... ##: ##:::::::
. ###. ###::::: ########: ########::. ######::: ########:::: ########:: ##:::: ##: ########:
:...::...::::::........::........::::......::::........:::::........:::..:::::..::........::
*/

void label_balance_edges_maxcut_weighted(
  pulp_graph_t& g, int num_parts, int* parts,
  int edge_outer_iter, int edge_balance_iter, int edge_refine_iter,
  double vert_balance, double edge_balance)
{
  int      num_verts = g.n;
  unsigned num_edges = g.m;
  unsigned cut_size  = 0;

  bool has_vwgts = (g.vertex_weights != NULL);
  bool has_ewgts = (g.edge_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;

  int*      part_sizes      = new int[num_parts];
  unsigned* part_edge_sizes = new unsigned[num_parts];
  int*      part_cut_sizes  = new int[num_parts];

  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;
  for (int i = 0; i < num_parts; ++i)
    part_edge_sizes[i] = 0;
  for (int i = 0; i < num_parts; ++i)
    part_cut_sizes[i] = 0;

  double avg_size      = (double)g.vertex_weights_sum / (double)num_parts;
  double avg_edge_size = num_edges / num_parts;
  int    num_swapped_1 = 0;
  int    num_swapped_2 = 0;
  double max_e         = 0.0;
  double max_c         = 0.0;
  double running_max_e     = (double)num_edges;
  double weight_exponent_e = 1.0;
  double weight_exponent_c = 1.0;

  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int queue_size;
  int next_size;
  int t = 0;
  int num_tries = 0;

  #pragma omp parallel
  {

    int*      part_sizes_thread      = new int[num_parts];
    unsigned* part_edge_sizes_thread = new unsigned[num_parts];
    int*      part_cut_sizes_thread  = new int[num_parts];
    unsigned  cut_size_thread = 0;

    for (int i = 0; i < num_parts; ++i)
      part_sizes_thread[i] = 0;
    for (int i = 0; i < num_parts; ++i)
      part_edge_sizes_thread[i] = 0;
    for (int i = 0; i < num_parts; ++i)
      part_cut_sizes_thread[i] = 0;

    #pragma omp for schedule(guided) nowait
    for (int i = 0; i < num_verts; ++i)
    {
      int part = parts[i];
      unsigned out_degree = out_degree(g, i);
      if (has_vwgts) part_sizes_thread[part] += g.vertex_weights[i];
      else ++part_sizes_thread[part];
      part_edge_sizes_thread[part] += out_degree;

      int* outs = out_vertices(g, i);
      int* weights = out_weights(g, i);
      for (unsigned j = 0; j < out_degree; ++j)
      {
        int out = outs[j];
        int out_part = parts[out];
        if (out_part != part)
        {
          if (has_ewgts)
          {
            part_cut_sizes_thread[part] += weights[j];
            cut_size_thread += weights[j];
          }
          else
          {
            ++part_cut_sizes_thread[part];
            ++cut_size_thread;
          }
        }
      }
    }

    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_sizes[i] += part_sizes_thread[i];
    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_edge_sizes[i] += part_edge_sizes_thread[i];
    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_cut_sizes[i] += part_cut_sizes_thread[i];
    #pragma omp atomic
    cut_size += cut_size_thread;

    delete [] part_sizes_thread;
    delete [] part_edge_sizes_thread;
    delete [] part_cut_sizes_thread;


    double avg_cut_size = (double)cut_size / (double)num_parts;
    double* part_counts       = new double[num_parts];
    double* part_weights      = new double[num_parts];
    double* part_edge_weights = new double[num_parts];
    double* part_cut_weights  = new double[num_parts];

    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;


    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      in_queue_next[i] = false;

    while(t < edge_outer_iter)
    {

      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;

      #pragma omp single
      {
        num_swapped_1 = 0;
        queue_size    = num_verts;
        next_size     = 0;

        max_e = 0.0;
        max_c = 0.0;
        for (int p = 0; p < num_parts; ++p)
        {
          if ((double)part_edge_sizes[p] / avg_edge_size > max_e)
            max_e = (double)part_edge_sizes[p] / avg_edge_size;
          if ((double)part_cut_sizes[p] / avg_cut_size > max_c)
            max_c = (double)part_cut_sizes[p] / avg_cut_size;
        }
        if (max_e < edge_balance)
        {
          max_e              = edge_balance;
          weight_exponent_e  = 1.0;
          weight_exponent_c *= max_c;
        }
        else
        {
          weight_exponent_e *= max_e / edge_balance;
          weight_exponent_c  = 1.0;
        }
      }

      int num_iter = 0;
      while (/*swapped &&*/ num_iter < edge_balance_iter)
      {
        for (int p = 0; p < num_parts; ++p)
        {
          part_weights[p]      = vert_balance * avg_size / (double)part_sizes[p] - 1.0;
          part_edge_weights[p] = max_e * avg_edge_size / (double)part_edge_sizes[p] - 1.0;
          part_cut_weights[p]  = max_c * avg_cut_size / (double)part_cut_sizes[p] - 1.0;
          if (part_weights[p] < 0.0)
            part_weights[p] = 0.0;
          if (part_edge_weights[p] < 0.0)
            part_edge_weights[p] = 0.0;
          if (part_cut_weights[p] < 0.0)
            part_cut_weights[p] = 0.0;
        }

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

          unsigned out_degree  = out_degree(g, v);
          int*     outs        = out_vertices(g, v);
          int*     weights     = out_weights(g, v);
          int      sum_weights = 0;
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int    out        = outs[j];
            int    part_out   = parts[out];
            double weight_out = 1.0;
            if (has_ewgts) weight_out = (double)weights[j];
            part_counts[part_out] += weight_out;
            sum_weights           += (int)weight_out;
          }

          int    max_part   = part;
          double max_val    = 0.0;
          int    part_count = (int)part_counts[part];
          int    max_count  = 0;
          for (int p = 0; p < num_parts; ++p)
          {
            int count_init = (int)part_counts[p];
            if (part_weights[p] > 0.0 && part_edge_weights[p] > 0.0 && part_cut_weights[p] > 0.0)
              part_counts[p] *= (part_edge_weights[p]*weight_exponent_e * part_cut_weights[p]*weight_exponent_c);
            else
              part_counts[p] = 0.0;

            if (part_counts[p] > max_val)
            {
              max_val   = part_counts[p];
              max_count = count_init;
              max_part  = p;
            }
          }

          if (max_part != part)
          {
            parts[v] = max_part;
            ++num_swapped_1;
            int diff_part     = 2*part_count - sum_weights;
            int diff_max_part = sum_weights - 2*max_count;
            int diff_cut      = diff_part + diff_max_part;
            #pragma omp atomic
            cut_size += diff_cut;
            #pragma omp atomic
            part_cut_sizes[part] += diff_part;
            #pragma omp atomic
            part_cut_sizes[max_part] += diff_max_part;
            #pragma omp atomic
            part_sizes[part] -= v_weight;
            #pragma omp atomic
            part_sizes[max_part] += v_weight;
            #pragma omp atomic
            part_edge_sizes[part] -= out_degree;
            #pragma omp atomic
            part_edge_sizes[max_part] += out_degree;

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

            avg_cut_size            = cut_size / num_parts;
            part_weights[part]      = vert_balance * avg_size / (double)part_sizes[part] - 1.0;
            part_edge_weights[part] = max_e * avg_edge_size / (double)part_edge_sizes[part] - 1.0;
            part_cut_weights[part]  = max_c * avg_cut_size / (double)part_cut_sizes[part] - 1.0;

            part_weights[max_part]      = vert_balance * avg_size / (double)part_sizes[max_part]  - 1.0;
            part_edge_weights[max_part] = max_e * avg_edge_size / (double)part_edge_sizes[max_part] - 1.0;
            part_cut_weights[max_part]  = max_c * avg_cut_size / (double)part_cut_sizes[max_part] - 1.0;

            if (part_weights[part] < 0.0)
              part_weights[part] = 0.0;
            if (part_edge_weights[part] < 0.0)
              part_edge_weights[part] = 0.0;
            if (part_cut_weights[part] < 0.0)
              part_cut_weights[part] = 0.0;

            if (part_weights[max_part] < 0.0)
              part_weights[max_part] = 0.0;
            if (part_edge_weights[max_part] < 0.0)
              part_edge_weights[max_part] = 0.0;
            if (part_cut_weights[max_part] < 0.0)
              part_cut_weights[max_part] = 0.0;
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
            printf("%d -- V: %2.2lf   E: %2.2lf, %lf   C: %2.2lf, %lf\n", num_swapped_1, vert_balance, max_e, weight_exponent_e, max_c, weight_exponent_c);
          #endif
          int* temp  = queue;
          queue      = queue_next;
          queue_next = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          max_e = 0.0;
          max_c = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_edge_sizes[p] / avg_edge_size > max_e)
              max_e = (double)part_edge_sizes[p] / avg_edge_size;
            if ((double)part_cut_sizes[p] / avg_cut_size > max_c)
              max_c = (double)part_cut_sizes[p] / avg_cut_size;
          }
          if (max_e < edge_balance)
          {
            max_e              = edge_balance;
            weight_exponent_e  = 1.0;
            weight_exponent_c *= max_c;
          }
          else
          {
            weight_exponent_e *= max_e / edge_balance;
            weight_exponent_c  = 1.0;
          }

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
      while (/*swapped &&*/ num_iter < edge_refine_iter)
      {
        for (int p = 0; p < num_parts; ++p)
        {
          part_weights[p]      = vert_balance * avg_size / (double)part_sizes[p] - 1.0;
          part_edge_weights[p] = max_e * avg_edge_size / (double)part_edge_sizes[p] - 1.0;
          part_cut_weights[p]  = max_c * avg_cut_size / (double)part_cut_sizes[p] - 1.0;
          if (part_weights[p] < 0.0)
            part_weights[p] = 0.0;
          if (part_edge_weights[p] < 0.0)
            part_edge_weights[p] = 0.0;
          if (part_cut_weights[p] < 0.0)
            part_cut_weights[p] = 0.0;
        }

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

          unsigned out_degree  = out_degree(g, v);
          int*     outs        = out_vertices(g, v);
          int*     weights     = out_weights(g, v);
          int      sum_weights = 0;
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int    out        = outs[j];
            int    part_out   = parts[out];
            double weight_out = 1.0;
            if (has_ewgts) weight_out = (double)weights[j];
            part_counts[part_out] += weight_out;
            sum_weights += (int)weight_out;
          }

          int max_part   = -1;
          int max_count  = -1;
          int part_count = part_counts[part];
          for (int p = 0; p < num_parts; ++p)
            if (part_counts[p] > max_count)
            {
              max_count = part_counts[p];
              max_part = p;
            }

          if (max_part != part)
          {
            int    diff_part        = 2*part_count - sum_weights;
            int    diff_max_part    = sum_weights - 2*max_count;
            int    diff_cut         = diff_part + diff_max_part;
            double new_max_imb      = (double)(part_sizes[max_part] + v_weight) / avg_size;
            double new_max_edge_imb = (double)(part_edge_sizes[max_part] + out_degree) / avg_edge_size;
            double new_max_cut_imb  = (double)(part_cut_sizes[max_part] + diff_max_part) / avg_cut_size;
            double new_cut_imb      = (double)(part_cut_sizes[part] + diff_part) / avg_cut_size;

            if (new_max_imb < vert_balance &&
                new_max_edge_imb < max_e &&
                new_max_cut_imb < max_c &&
                new_cut_imb < max_c)
            {
              ++num_swapped_2;
              parts[v] = max_part;

              #pragma omp atomic
              cut_size += diff_cut;

              #pragma omp atomic
              part_cut_sizes[part]     += diff_part;
              #pragma omp atomic
              part_cut_sizes[max_part] += diff_max_part;

              #pragma omp atomic
              ++part_sizes[max_part];
              #pragma omp atomic
              --part_sizes[part];

              #pragma omp atomic
              part_edge_sizes[max_part] += out_degree;
              #pragma omp atomic
              part_edge_sizes[part]     -= out_degree;

              avg_cut_size = cut_size / num_parts;

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
            printf("%d -- V: %2.2lf  E: %2.2lf  C: %2.2lf\n", num_swapped_2, vert_balance, max_e, max_c);
          #endif
          int* temp  = queue;
          queue      = queue_next;
          queue_next = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          num_swapped_2 = 0;

          max_e = 0.0;
          max_c = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_edge_sizes[p] / avg_edge_size > max_e)
              max_e = (double)part_edge_sizes[p] / avg_edge_size;
            if ((double)part_cut_sizes[p] / avg_cut_size > max_c)
              max_c = (double)part_cut_sizes[p] / avg_cut_size;
          }
          if (max_e < edge_balance)
          {
            max_e              = edge_balance;
            weight_exponent_e  = 1.0;
            weight_exponent_c *= max_c;
          }
          else
          {
            weight_exponent_e *= max_e / edge_balance;
            weight_exponent_c  = 1.0;
          }

          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        }
      }

      #pragma omp single
      {
        if (max_e > edge_balance*1.01 && t == edge_outer_iter-1 && num_tries < MAX_TRIES)
        {
          --t;
          if (max_e < running_max_e*0.99)
          {
            running_max_e = max_e;
            printf("Edge balance missed, attempting further iterations: (%2.3lf)\n", max_e);
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
    delete [] part_edge_weights;
    delete [] part_cut_weights;

  } // end par

  delete [] part_sizes;
  delete [] part_edge_sizes;
  delete [] part_cut_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}








/*
'##:::::'##::::'####::::'########:'########:::'######:::'########::::'########:::::'###::::'##:::::::
 ##:'##: ##::::. ##::::: ##.....:: ##.... ##:'##... ##:: ##.....::::: ##.... ##:::'## ##::: ##:::::::
 ##: ##: ##::::: ##::::: ##::::::: ##:::: ##: ##:::..::: ##:::::::::: ##:::: ##::'##:. ##:: ##:::::::
 ##: ##: ##::::: ##::::: ######::: ##:::: ##: ##::'####: ######:::::: ########::'##:::. ##: ##:::::::
 ##: ##: ##::::: ##::::: ##...:::: ##:::: ##: ##::: ##:: ##...::::::: ##.... ##: #########: ##:::::::
 ##: ##: ##::::: ##::::: ##::::::: ##:::: ##: ##::: ##:: ##:::::::::: ##:::: ##: ##.... ##: ##:::::::
. ###. ###:::::'####:::: ########: ########::. ######::: ########:::: ########:: ##:::: ##: ########:
:...::...::::::....:::::........::........::::......::::........:::::........:::..:::::..::........::
*/

void label_balance_edges_maxcut_weighted_interpart(
  pulp_graph_t& g, int num_parts, int* parts,
  int edge_outer_iter, int edge_balance_iter, int edge_refine_iter,
  double vert_balance, double edge_balance)
{
  int      num_verts = g.n;
  unsigned num_edges = g.m;
  unsigned cut_size  = 0;

  bool has_ewgts = (g.edge_weights != NULL);
  bool has_vwgts = (g.vertex_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;

  int*      part_sizes      = new int[num_parts];
  unsigned* part_edge_sizes = new unsigned[num_parts];
  int*      part_cut_sizes  = new int[num_parts];

  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;
  for (int i = 0; i < num_parts; ++i)
    part_edge_sizes[i] = 0;
  for (int i = 0; i < num_parts; ++i)
    part_cut_sizes[i] = 0;

  double avg_size          = (double)g.vertex_weights_sum / (double)num_parts;
  double avg_edge_size     = num_edges / num_parts;
  int    num_swapped_1     = 0;
  int    num_swapped_2     = 0;
  double max_e             = 0.0;
  double max_c             = 0.0;
  double running_max_e     = (double)num_edges; // used to track the max edge balance
  double weight_exponent_e = 1.0;
  double weight_exponent_c = 1.0;

  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int queue_size;
  int next_size;
  int t         = 0;
  int num_tries = 0;

  #pragma omp parallel
  {

    int*      part_sizes_thread      = new int[num_parts];
    unsigned* part_edge_sizes_thread = new unsigned[num_parts];
    int*      part_cut_sizes_thread  = new int[num_parts];
    unsigned  cut_size_thread        = 0;

    for (int i = 0; i < num_parts; ++i)
      part_sizes_thread[i] = 0;
    for (int i = 0; i < num_parts; ++i)
      part_edge_sizes_thread[i] = 0;
    for (int i = 0; i < num_parts; ++i)
      part_cut_sizes_thread[i] = 0;

    #pragma omp for schedule(guided) nowait
    for (int i = 0; i < num_verts; ++i)
    {
      int      part       = parts[i];
      unsigned out_degree = out_degree(g, i); // get the out degree of vertex i
      if (has_vwgts)
        part_sizes_thread[part] += g.vertex_weights[i];
      else
        ++part_sizes_thread[part];

      // Update the edge sizes
      part_edge_sizes_thread[part] += out_degree;

      int* outs                   = out_vertices(g, i);
      int* weights                = out_weights(g, i);
      int* partition_comm_weights = out_interpart_weights(g, part, num_parts);  // Notes(julie): added

      for (unsigned j = 0; j < out_degree; ++j)
      {
        int out      = outs[j];
        int out_part = parts[out];
        if (out_part != part)
        {
          if (has_ewgts /* && has_label_weight*/) // always have label weight
          {
            part_cut_sizes_thread[part] += weights[j] * partition_comm_weights[out_part];
            cut_size_thread             += weights[j] * partition_comm_weights[out_part];
          }
          else
          {
            ++part_cut_sizes_thread[part];
            ++cut_size_thread;
          }
        }
      }
    }

    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_sizes[i] += part_sizes_thread[i];
    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_edge_sizes[i] += part_edge_sizes_thread[i];
    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_cut_sizes[i] += part_cut_sizes_thread[i];
    #pragma omp atomic
    cut_size += cut_size_thread;

    delete [] part_sizes_thread;
    delete [] part_edge_sizes_thread;
    delete [] part_cut_sizes_thread;


    double  avg_cut_size      = (double)cut_size / (double)num_parts;

	// ========================================================================
	// OUTER LOOP ITERATION
	// ========================================================================
    double* part_counts       = new double[num_parts];
    double* part_weights      = new double[num_parts];
    double* part_edge_weights = new double[num_parts];
    double* part_cut_weights  = new double[num_parts];

    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;


    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      in_queue_next[i] = false;

    while(t < edge_outer_iter)
    {

      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;

      // ==================================================
      //
      //  #####    ##   #        ##   #    #  ####  ######
      //  #    #  #  #  #       #  #  ##   # #    # #
      //  #####  #    # #      #    # # #  # #      #####
      //  #    # ###### #      ###### #  # # #      #
      //  #    # #    # #      #    # #   ## #    # #
      //  #####  #    # ###### #    # #    #  ####  ######
      //
      // ==================================================
      #pragma omp single
      {
        num_swapped_1 = 0;
        queue_size    = num_verts;
        next_size     = 0;

        max_e = 0.0;
        max_c = 0.0;
        for (int p = 0; p < num_parts; ++p)
        {
          if ((double)part_edge_sizes[p] / avg_edge_size > max_e)
            max_e = (double)part_edge_sizes[p] / avg_edge_size;

          if ((double)part_cut_sizes[p] / avg_cut_size > max_c)
            max_c = (double)part_cut_sizes[p] / avg_cut_size;
        }

        if (max_e < edge_balance)
        {
          max_e              = edge_balance;
          weight_exponent_e  = 1.0;
          weight_exponent_c *= max_c;
        }
        else
        {
          weight_exponent_e *= max_e / edge_balance;
          weight_exponent_c  = 1.0;
        }
      } // end single

      int num_iter = 0;
      while (/*swapped &&*/ num_iter < edge_balance_iter)
      {
        for (int p = 0; p < num_parts; ++p)
        {
          part_weights[p]      = vert_balance * avg_size / (double)part_sizes[p] - 1.0;
          part_edge_weights[p] = max_e * avg_edge_size / (double)part_edge_sizes[p] - 1.0;
          part_cut_weights[p]  = max_c * avg_cut_size / (double)part_cut_sizes[p] - 1.0;
          if (part_weights[p] < 0.0)
            part_weights[p] = 0.0;
          if (part_edge_weights[p] < 0.0)
            part_edge_weights[p] = 0.0;
          if (part_cut_weights[p] < 0.0)
            part_cut_weights[p] = 0.0;
        }

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

          unsigned out_degree  = out_degree(g, v);
          int*     outs        = out_vertices(g, v);
          int*     weights     = out_weights(g, v);
          int      sum_weights = 0;
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out = outs[j];
            int part_out = parts[out];
            double weight_out = 1.0;
            if (has_ewgts) weight_out = (double)weights[j];
            part_counts[part_out] += weight_out;
            sum_weights += (int)weight_out;
          }

          // -----------------------------------------------------
          // Inter-partition weight adjustment
          // -----------------------------------------------------
          int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
          for (int p = 0; p < num_parts; ++p)
            part_counts[p] *= partition_comm_weights[p];


          int max_part = part;
          double max_val = 0.0;
          int part_count = (int)part_counts[part];
          int max_count = 0;
          for (int p = 0; p < num_parts; ++p)
          {
            int count_init = (int)part_counts[p];
            if (part_weights[p] > 0.0 &&
                part_edge_weights[p] > 0.0 &&
                part_cut_weights[p] > 0.0)
              part_counts[p] *= (part_edge_weights[p]*weight_exponent_e * part_cut_weights[p]*weight_exponent_c);
            else
              part_counts[p] = 0.0;

            if (part_counts[p] > max_val)
            {
              max_val = part_counts[p];
              max_count = count_init;
              max_part = p;
            }
          }

          if (max_part != part)
          {

            // Unless in bin packing mode, do not move a vertex if it would empty its partition
            if (!g.do_bin_packing && (part_sizes[part] - v_weight <= 0))
              continue;

            if (g.max_partition_size > 0 &&
                (part_sizes[max_part]+v_weight) > (g.max_partition_size*g.partition_capacities[max_part]))
              continue;

            parts[v] = max_part;
            ++num_swapped_1;
            int diff_part     = 2 * part_count - sum_weights;
            int diff_max_part = sum_weights - 2 * max_count;
            int diff_cut      = diff_part + diff_max_part;
            #pragma omp atomic
            cut_size += diff_cut;
            #pragma omp atomic
            part_cut_sizes[part]     += diff_part;
            #pragma omp atomic
            part_cut_sizes[max_part] += diff_max_part;
            #pragma omp atomic
            part_sizes[part]     -= v_weight;
            #pragma omp atomic
            part_sizes[max_part] += v_weight;
            #pragma omp atomic
            part_edge_sizes[part] -= out_degree;
            #pragma omp atomic
            part_edge_sizes[max_part] += out_degree;

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

            avg_cut_size            = cut_size / num_parts;
            part_weights[part]      = vert_balance * avg_size / (double)part_sizes[part] - 1.0;
            part_edge_weights[part] = max_e * avg_edge_size / (double)part_edge_sizes[part] - 1.0;
            part_cut_weights[part]  = max_c * avg_cut_size / (double)part_cut_sizes[part] - 1.0;

            part_weights[max_part]      = vert_balance * avg_size / (double)part_sizes[max_part]  - 1.0;
            part_edge_weights[max_part] = max_e * avg_edge_size / (double)part_edge_sizes[max_part] - 1.0;
            part_cut_weights[max_part]  = max_c * avg_cut_size / (double)part_cut_sizes[max_part] - 1.0;

            if (part_weights[part] < 0.0)
              part_weights[part] = 0.0;
            if (part_edge_weights[part] < 0.0)
              part_edge_weights[part] = 0.0;
            if (part_cut_weights[part] < 0.0)
              part_cut_weights[part] = 0.0;

            if (part_weights[max_part] < 0.0)
              part_weights[max_part] = 0.0;
            if (part_edge_weights[max_part] < 0.0)
              part_edge_weights[max_part] = 0.0;
            if (part_cut_weights[max_part] < 0.0)
              part_cut_weights[max_part] = 0.0;
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
			      printf("\nnum_swapped_1 (edge, balance): %d -- V: %2.2lf   E: %2.2lf, %lf   C: %2.2lf, %lf\n", num_swapped_1, vert_balance, max_e, weight_exponent_e, max_c, weight_exponent_c);
          #endif

          int*  temp    = queue;
          queue         = queue_next;
          queue_next    = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          max_e = 0.0;
          max_c = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_edge_sizes[p] / avg_edge_size > max_e)
              max_e = (double)part_edge_sizes[p] / avg_edge_size;
            if ((double)part_cut_sizes[p] / avg_cut_size > max_c)
              max_c = (double)part_cut_sizes[p] / avg_cut_size;
          }
          if (max_e < edge_balance)
          {
            max_e = edge_balance;
            weight_exponent_e = 1.0;
            weight_exponent_c *= max_c;
          }
          else
          {
            weight_exponent_e *= max_e / edge_balance;
            weight_exponent_c = 1.0;
          }

          num_swapped_1 = 0;

          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        }
      } // end while



		// ========================================================================
		//  ######
		//  #     # ###### ###### # #    # ###### #    # ###### #    # #####
		//  #     # #      #      # ##   # #      ##  ## #      ##   #   #
		//  ######  #####  #####  # # #  # #####  # ## # #####  # #  #   #
		//  #   #   #      #      # #  # # #      #    # #      #  # #   #
		//  #    #  #      #      # #   ## #      #    # #      #   ##   #
		//  #     # ###### #      # #    # ###### #    # ###### #    #   #
		// ========================================================================

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
      while (/*swapped &&*/ num_iter < edge_refine_iter)
      {
        for (int p = 0; p < num_parts; ++p)
        {
          part_weights[p]      = vert_balance * avg_size / (double)part_sizes[p] - 1.0;
          part_edge_weights[p] = max_e * avg_edge_size / (double)part_edge_sizes[p] - 1.0;
          part_cut_weights[p]  = max_c * avg_cut_size / (double)part_cut_sizes[p] - 1.0;

          if (part_weights[p] < 0.0)
            part_weights[p] = 0.0;
          if (part_edge_weights[p] < 0.0)
            part_edge_weights[p] = 0.0;
          if (part_cut_weights[p] < 0.0)
            part_cut_weights[p] = 0.0;
        }

        #pragma omp for schedule(guided) reduction(+:num_swapped_2) nowait
        for (int i = 0; i < queue_size; ++i)
        {
          int v       = queue[i];
          in_queue[v] = false;
          int part    = parts[v];
          int v_weight = 1;
          if (has_vwgts) v_weight = g.vertex_weights[v];

          for (int p = 0; p < num_parts; ++p)
            part_counts[p] = 0;

          unsigned out_degree  = out_degree(g, v);
          int*     outs        = out_vertices(g, v);
          int*     weights     = out_weights(g, v);
          int      sum_weights = 0;
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out = outs[j];
            int part_out = parts[out];
            double weight_out = 1.0;
            if (has_ewgts) weight_out = (double)weights[j];
            part_counts[part_out] += weight_out;
            sum_weights += (int)weight_out;
          }

          // -----------------------------------------------------
          // Inter-partition weight adjustment
          // -----------------------------------------------------
          int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
          for (int p = 0; p < num_parts; ++p)
            part_counts[p] *= partition_comm_weights[p];

          int max_part   = -1;
          int max_count  = -1;
          int part_count = part_counts[part];
          for (int p = 0; p < num_parts; ++p)
          {
            if (part_counts[p] > max_count)
            {
              max_count = part_counts[p];
              max_part = p;
            }
          }

  				// SWAP IF IMPROVEMENT (only if it improves the balance ratios)
          if (max_part != part)
          {
            // Unless in bin packing mode, do not move a vertex if it would empty its partition
            if (!g.do_bin_packing && (part_sizes[part] - v_weight <= 0))
              continue;

            if (g.max_partition_size > 0 &&
                (part_sizes[max_part]+v_weight) > (g.max_partition_size*g.partition_capacities[max_part]))
              continue;

            int    diff_part        = 2*part_count - sum_weights;
            int    diff_max_part    = sum_weights - 2*max_count;
            int    diff_cut         = diff_part + diff_max_part;
            double new_max_imb      = (double)(part_sizes[max_part] + v_weight) / avg_size;
            double new_max_edge_imb = (double)(part_edge_sizes[max_part] + out_degree) / avg_edge_size;
            double new_max_cut_imb  = (double)(part_cut_sizes[max_part] + diff_max_part) / avg_cut_size;
            double new_cut_imb      = (double)(part_cut_sizes[part] + diff_part) / avg_cut_size;
            if ( new_max_imb < vert_balance &&
                new_max_edge_imb < max_e &&
                 new_max_cut_imb < max_c &&
                 new_cut_imb < max_c)
            {
              ++num_swapped_2;
              parts[v] = max_part;
              #pragma omp atomic
              cut_size += diff_cut;
              #pragma omp atomic
              part_cut_sizes[part] += diff_part;
              #pragma omp atomic
              part_cut_sizes[max_part] += diff_max_part;
              #pragma omp atomic
              ++part_sizes[max_part];
              #pragma omp atomic
              --part_sizes[part];
              #pragma omp atomic
              part_edge_sizes[max_part] += out_degree;
              #pragma omp atomic
              part_edge_sizes[part] -= out_degree;

              avg_cut_size = cut_size / num_parts;

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
            printf("\nnum_swapped_2 (edge, refine): %d -- V: %2.2lf  E: %2.2lf  C: %2.2lf\n", num_swapped_2, vert_balance, max_e, max_c);
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

          max_e = 0.0;
          max_c = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_edge_sizes[p] / avg_edge_size > max_e)
              max_e = (double)part_edge_sizes[p] / avg_edge_size;
            if ((double)part_cut_sizes[p] / avg_cut_size > max_c)
              max_c = (double)part_cut_sizes[p] / avg_cut_size;
          }
          if (max_e < edge_balance)
          {
            max_e = edge_balance;
            weight_exponent_e = 1.0;
            weight_exponent_c *= max_c;
          }
          else
          {
            weight_exponent_e *= max_e / edge_balance;
            weight_exponent_c = 1.0;
          }

          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        }
      }

      #pragma omp single
      {
        if (max_e > edge_balance*1.01 && t == edge_outer_iter-1 && num_tries < MAX_TRIES)
        {
          --t;
          if (max_e < running_max_e*0.99)
          {
            running_max_e = max_e;
            printf("Edge balance missed, attempting further iterations: (%2.3lf)\n", max_e);
          }
          else
            ++num_tries;
        }
        else
          ++t;
      }

    } // end while

    delete [] part_counts;
    delete [] part_weights;
    delete [] part_edge_weights;
    delete [] part_cut_weights;

  } // end omp parallel


  delete [] part_sizes;
  delete [] part_edge_sizes;
  delete [] part_cut_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}








/*
'##:::::'##:::::::::'####::::::::::'######:::::'###::::'########::::::::::'########:'########:::'######:::'########::::'########:::::'###::::'##:::::::
 ##:'##: ##:::::::::. ##::::::::::'##... ##:::'## ##::: ##.... ##::::::::: ##.....:: ##.... ##:'##... ##:: ##.....::::: ##.... ##:::'## ##::: ##:::::::
 ##: ##: ##:::::::::: ##:::::::::: ##:::..:::'##:. ##:: ##:::: ##::::::::: ##::::::: ##:::: ##: ##:::..::: ##:::::::::: ##:::: ##::'##:. ##:: ##:::::::
 ##: ##: ##:::::::::: ##:::::::::: ##:::::::'##:::. ##: ########:::::::::: ######::: ##:::: ##: ##::'####: ######:::::: ########::'##:::. ##: ##:::::::
 ##: ##: ##:::::::::: ##:::::::::: ##::::::: #########: ##.....::::::::::: ##...:::: ##:::: ##: ##::: ##:: ##...::::::: ##.... ##: #########: ##:::::::
 ##: ##: ##:'###::::: ##::'###:::: ##::: ##: ##.... ##: ##::::::::'###:::: ##::::::: ##:::: ##: ##::: ##:: ##:::::::::: ##:::: ##: ##.... ##: ##:::::::
. ###. ###:: ###::::'####: ###::::. ######:: ##:::: ##: ##:::::::: ###:::: ########: ########::. ######::: ########:::: ########:: ##:::: ##: ########:
:...::...:::...:::::....::...::::::......:::..:::::..::..:::::::::...:::::........::........::::......::::........:::::........:::..:::::..::........::
*/

void
label_balance_edges_maxcut_weighted_interpart_capacitated(pulp_graph_t& g, int num_parts, int* parts,
                                                          int edge_outer_iter, int edge_balance_iter, int edge_refine_iter,
                                                          double vert_balance, double edge_balance)
{
  bool has_ipwgts        = (g.interpartition_weights != NULL);
  bool has_p_capacities  = (g.partition_capacities != NULL);
  if (!has_p_capacities || !has_ipwgts)
  {
    printf("Error: partition (label) capacities and inter-part weights required.\n");
    exit(0);
  }

  int      num_verts = g.n;
  unsigned num_edges = g.m;
  unsigned cut_size  = 0;

  bool has_ewgts = (g.edge_weights != NULL);
  bool has_vwgts = (g.vertex_weights != NULL);
  if (!has_vwgts) g.vertex_weights_sum = g.n;

  int*      part_sizes      = new int[num_parts];
  unsigned* part_edge_sizes = new unsigned[num_parts];
  int*      part_cut_sizes  = new int[num_parts];

  for (int i = 0; i < num_parts; ++i)
    part_sizes[i] = 0;
  for (int i = 0; i < num_parts; ++i)
    part_edge_sizes[i] = 0;
  for (int i = 0; i < num_parts; ++i)
    part_cut_sizes[i] = 0;

  // Note(julie9): Part size and Edge cut count set according to the partition capacities
  double  unit_avg_size      = (double)g.vertex_weights_sum / (double)g.partition_capacities_sum;
  double  unit_avg_edge_size = (double)num_edges / (double)g.partition_capacities_sum;
  double  unit_avg_cut_size  = (double)cut_size / (double)g.partition_capacities_sum;
  double* avg_sizes          = new double[num_parts];
  double* avg_edge_sizes     = new double[num_parts];
  double* avg_cut_sizes      = new double[num_parts];
  for (int i = 0; i < num_parts; ++i)
    avg_sizes[i] = unit_avg_size * g.partition_capacities[i];
  for (int i = 0; i < num_parts; ++i)
    avg_edge_sizes[i] = unit_avg_edge_size * g.partition_capacities[i];
  for (int i = 0; i < num_parts; ++i)
    avg_cut_sizes[i] = unit_avg_cut_size * g.partition_capacities[i];

  int    num_swapped_1     = 0;
  int    num_swapped_2     = 0;
  double max_e             = 0.0;
  double max_c             = 0.0;
  double running_max_e     = (double)num_edges; // used to track the max edge balance
  double weight_exponent_e = 1.0;
  double weight_exponent_c = 1.0;

  int*  queue         = new int[num_verts*QUEUE_MULTIPLIER];
  int*  queue_next    = new int[num_verts*QUEUE_MULTIPLIER];
  bool* in_queue      = new bool[num_verts];
  bool* in_queue_next = new bool[num_verts];
  int queue_size;
  int next_size;
  int t         = 0;
  int num_tries = 0;

  #pragma omp parallel
  {

    int*      part_sizes_thread      = new int[num_parts];
    unsigned* part_edge_sizes_thread = new unsigned[num_parts];
    int*      part_cut_sizes_thread  = new int[num_parts];
    unsigned  cut_size_thread        = 0;

    for (int i = 0; i < num_parts; ++i)
      part_sizes_thread[i] = 0;
    for (int i = 0; i < num_parts; ++i)
      part_edge_sizes_thread[i] = 0;
    for (int i = 0; i < num_parts; ++i)
      part_cut_sizes_thread[i] = 0;

    #pragma omp for schedule(guided) nowait
    for (int i = 0; i < num_verts; ++i)
    {
      int      part       = parts[i];
      unsigned out_degree = out_degree(g, i); // get the out degree of vertex i
      if (has_vwgts)
      part_sizes_thread[part] += g.vertex_weights[i];
      else
        ++part_sizes_thread[part];

      // Update the edge sizes
      part_edge_sizes_thread[part] += out_degree;

      int* outs          = out_vertices(g, i);
      int* weights       = out_weights(g, i);
      int* partition_comm_weights = out_interpart_weights(g, part, num_parts);  // Notes(julie): added

      for (unsigned j = 0; j < out_degree; ++j)
      {
        int out      = outs[j];
        int out_part = parts[out];
        if (out_part != part)
        {
          if (has_ewgts)
        {
            part_cut_sizes_thread[part] += weights[j] * partition_comm_weights[out_part];
            cut_size_thread             += weights[j] * partition_comm_weights[out_part];
          }
          else
          {
            ++part_cut_sizes_thread[part];
            ++cut_size_thread;
          }
        }
      }
    }

    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_sizes[i] += part_sizes_thread[i];
    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_edge_sizes[i] += part_edge_sizes_thread[i];
    for (int i = 0; i < num_parts; ++i)
      #pragma omp atomic
      part_cut_sizes[i] += part_cut_sizes_thread[i];
    #pragma omp atomic
    cut_size += cut_size_thread;

    delete [] part_sizes_thread;
    delete [] part_edge_sizes_thread;
    delete [] part_cut_sizes_thread;

    unit_avg_cut_size = (double)cut_size / (double)g.partition_capacities_sum;
    for (int i = 0; i < num_parts; ++i)
      avg_cut_sizes[i] = unit_avg_cut_size * g.partition_capacities[i];

    // ========================================================================
    // OUTER LOOP ITERATION
    // ========================================================================
    double* part_counts       = new double[num_parts];
    double* part_weights      = new double[num_parts];
    double* part_edge_weights = new double[num_parts];
    double* part_cut_weights  = new double[num_parts];

    int thread_queue[ THREAD_QUEUE_SIZE ];
    int thread_queue_size = 0;
    int thread_start;


    #pragma omp for schedule(static) nowait
    for (int i = 0; i < num_verts; ++i)
      in_queue_next[i] = false;

    while(t < edge_outer_iter)
    {

      #pragma omp for schedule(static)
      for (int i = 0; i < num_verts; ++i)
        queue[i] = i;

      // ==================================================
      //
      //  #####    ##   #        ##   #    #  ####  ######
      //  #    #  #  #  #       #  #  ##   # #    # #
      //  #####  #    # #      #    # # #  # #      #####
      //  #    # ###### #      ###### #  # # #      #
      //  #    # #    # #      #    # #   ## #    # #
      //  #####  #    # ###### #    # #    #  ####  ######
      //
      // ==================================================
      #pragma omp single
      {
        num_swapped_1 = 0;
        queue_size    = num_verts;
        next_size     = 0;

        max_e = 0.0;
        max_c = 0.0;
        for (int p = 0; p < num_parts; ++p)
        {
          if ((double)part_edge_sizes[p] / avg_edge_sizes[p] > max_e)
            max_e = (double)part_edge_sizes[p] / avg_edge_sizes[p];

          if ((double)part_cut_sizes[p] / avg_cut_sizes[p] > max_c)
            max_c = (double)part_cut_sizes[p] / avg_cut_sizes[p];
        }

        if (max_e < edge_balance)
        {
          max_e              = edge_balance;
          weight_exponent_e  = 1.0;
          weight_exponent_c *= max_c;
        }
        else
        {
          weight_exponent_e *= max_e / edge_balance;
          weight_exponent_c  = 1.0;
        }
      } // end single

      int num_iter = 0;
      while (/*swapped &&*/ num_iter < edge_balance_iter)
      {
        for (int p = 0; p < num_parts; ++p)
        {
          part_weights[p]      = vert_balance * avg_sizes[p] / (double)part_sizes[p] - 1.0;
          part_edge_weights[p] = max_e * avg_edge_sizes[p] / (double)part_edge_sizes[p] - 1.0;
          part_cut_weights[p]  = max_c * avg_cut_sizes[p] / (double)part_cut_sizes[p] - 1.0;
          if (part_weights[p] < 0.0)
            part_weights[p] = 0.0;
          if (part_edge_weights[p] < 0.0)
            part_edge_weights[p] = 0.0;
          if (part_cut_weights[p] < 0.0)
            part_cut_weights[p] = 0.0;
        }

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

          unsigned out_degree  = out_degree(g, v);
          int*     outs        = out_vertices(g, v);
          int*     weights     = out_weights(g, v);
          int      sum_weights = 0;
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out = outs[j];
            int part_out = parts[out];
            double weight_out = 1.0;
            if (has_ewgts) weight_out = (double)weights[j];
            part_counts[part_out] += weight_out;
            sum_weights += (int)weight_out;
          }

          // -----------------------------------------------------
          // Inter-partition weight adjustment
          // -----------------------------------------------------
          int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
          for (int p = 0; p < num_parts; ++p)
            part_counts[p] *= partition_comm_weights[p];


          int max_part = part;
          double max_val = 0.0;
          int part_count = (int)part_counts[part];
          int max_count = 0;
          for (int p = 0; p < num_parts; ++p)
          {
            int count_init = (int)part_counts[p];
            if (part_weights[p] > 0.0 &&
                part_edge_weights[p] > 0.0 &&
                part_cut_weights[p] > 0.0)
              part_counts[p] *= (part_edge_weights[p]*weight_exponent_e * part_cut_weights[p]*weight_exponent_c);
            else
              part_counts[p] = 0.0;

            if (part_counts[p] > max_val)
            {
              max_val = part_counts[p];
              max_count = count_init;
              max_part = p;
            }
          }

          if (max_part != part)
          {
            // Unless in bin packing mode, do not move a vertex if it would empty its partition
            if (!g.do_bin_packing && (part_sizes[part] - v_weight <= 0))
              continue;

            if (g.max_partition_size > 0 &&
                (part_sizes[max_part]+v_weight) > (g.max_partition_size*g.partition_capacities[max_part]))
            {
              printf("WARNING: Max partition size exceeded on part %d: %d + %d > %d * %d\n",
                     max_part, part_sizes[max_part], v_weight, g.max_partition_size, g.partition_capacities[max_part]);
              continue;
            }

            parts[v] = max_part;
            ++num_swapped_1;

            int diff_part     = 2 * part_count - sum_weights;
            int diff_max_part = sum_weights - 2 * max_count;
            int diff_cut      = diff_part + diff_max_part;

            #pragma omp atomic
            cut_size += diff_cut;

            #pragma omp atomic
            part_cut_sizes[part]     += diff_part;
            #pragma omp atomic
            part_cut_sizes[max_part] += diff_max_part;

            #pragma omp atomic
            part_sizes[part]     -= v_weight;
            #pragma omp atomic
            part_sizes[max_part] += v_weight;

            #pragma omp atomic
            part_edge_sizes[part]     -= out_degree;
            #pragma omp atomic
            part_edge_sizes[max_part] += out_degree;

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

            unit_avg_cut_size = (double)cut_size / (double)g.partition_capacities_sum;
            for (int i = 0; i < num_parts; ++i)
              avg_cut_sizes[i] = unit_avg_cut_size * g.partition_capacities[i];

            part_weights[part]      = vert_balance * avg_sizes[part] / (double)part_sizes[part] - 1.0;
            part_edge_weights[part] = max_e * avg_edge_sizes[part] / (double)part_edge_sizes[part] - 1.0;
            part_cut_weights[part]  = max_c * avg_cut_sizes[part] / (double)part_cut_sizes[part] - 1.0;

            part_weights[max_part]      = vert_balance * avg_sizes[part] / (double)part_sizes[max_part]  - 1.0;
            part_edge_weights[max_part] = max_e * avg_edge_sizes[part] / (double)part_edge_sizes[max_part] - 1.0;
            part_cut_weights[max_part]  = max_c * avg_cut_sizes[part] / (double)part_cut_sizes[max_part] - 1.0;

            if (part_weights[part] < 0.0)
              part_weights[part] = 0.0;
            if (part_edge_weights[part] < 0.0)
              part_edge_weights[part] = 0.0;
            if (part_cut_weights[part] < 0.0)
              part_cut_weights[part] = 0.0;

            if (part_weights[max_part] < 0.0)
              part_weights[max_part] = 0.0;
            if (part_edge_weights[max_part] < 0.0)
              part_edge_weights[max_part] = 0.0;
            if (part_cut_weights[max_part] < 0.0)
              part_cut_weights[max_part] = 0.0;
          }
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
			      printf("\nnum_swapped_1 (edge, balance): %d -- V: %2.2lf   E: %2.2lf, %lf   C: %2.2lf, %lf\n", num_swapped_1, vert_balance, max_e, weight_exponent_e, max_c, weight_exponent_c);
          #endif

          int*  temp    = queue;
          queue         = queue_next;
          queue_next    = temp;
          bool* temp_b  = in_queue;
          in_queue      = in_queue_next;
          in_queue_next = temp_b;
          queue_size    = next_size;
          next_size     = 0;

          max_e = 0.0;
          max_c = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_edge_sizes[p] / avg_edge_sizes[p] > max_e)
              max_e = (double)part_edge_sizes[p] / avg_edge_sizes[p];
            if ((double)part_cut_sizes[p] / avg_cut_sizes[p] > max_c)
              max_c = (double)part_cut_sizes[p] / avg_cut_sizes[p];
          }
          if (max_e < edge_balance)
          {
            max_e              = edge_balance;
            weight_exponent_e  = 1.0;
            weight_exponent_c *= max_c;
          }
          else
          {
            weight_exponent_e *= max_e / edge_balance;
            weight_exponent_c  = 1.0;
          }

          num_swapped_1 = 0;

          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        }
      } // end while



		// ========================================================================
		//  ######
		//  #     # ###### ###### # #    # ###### #    # ###### #    # #####
		//  #     # #      #      # ##   # #      ##  ## #      ##   #   #
		//  ######  #####  #####  # # #  # #####  # ## # #####  # #  #   #
		//  #   #   #      #      # #  # # #      #    # #      #  # #   #
		//  #    #  #      #      # #   ## #      #    # #      #   ##   #
		//  #     # ###### #      # #    # ###### #    # ###### #    #   #
		// ========================================================================

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
      while (/*swapped &&*/ num_iter < edge_refine_iter)
      {
        for (int p = 0; p < num_parts; ++p)
        {
          part_weights[p]      = vert_balance * avg_sizes[p] / (double)part_sizes[p] - 1.0;
          part_edge_weights[p] = max_e * avg_edge_sizes[p] / (double)part_edge_sizes[p] - 1.0;
          part_cut_weights[p]  = max_c * avg_cut_sizes[p] / (double)part_cut_sizes[p] - 1.0;

          if (part_weights[p] < 0.0)
            part_weights[p] = 0.0;
          if (part_edge_weights[p] < 0.0)
            part_edge_weights[p] = 0.0;
          if (part_cut_weights[p] < 0.0)
            part_cut_weights[p] = 0.0;
        }

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

          unsigned out_degree  = out_degree(g, v);
          int*     outs        = out_vertices(g, v);
          int*     weights     = out_weights(g, v);
          int      sum_weights = 0;
          for (unsigned j = 0; j < out_degree; ++j)
          {
            int out = outs[j];
            int part_out = parts[out];
            double weight_out = 1.0;
            if (has_ewgts) weight_out = (double)weights[j];
            part_counts[part_out] += weight_out;
            sum_weights += (int)weight_out;
          }

          // -----------------------------------------------------
          // Inter-partition weight adjustment
          // -----------------------------------------------------
          int* partition_comm_weights = out_interpart_weights(g, part, num_parts);
          for (int p = 0; p < num_parts; ++p)
            part_counts[p] *= partition_comm_weights[p];

          int max_part   = -1;
          int max_count  = -1;
          int part_count = part_counts[part];
          for (int p = 0; p < num_parts; ++p)
          {
            if (part_counts[p] > max_count)
            {
              max_count = part_counts[p];
              max_part = p;
            }
          }

  				// SWAP IF IMPROVEMENT (only if it improves the balance ratios)
          if (max_part != part)
          {
            // Unless in bin packing mode, do not move a vertex if it would empty its partition
            if (!g.do_bin_packing && (part_sizes[part] - v_weight <= 0))
              continue;

            if (g.max_partition_size > 0 &&
                (part_sizes[max_part]+v_weight) > (g.max_partition_size*g.partition_capacities[max_part]))
            {
              printf("WARNING: Max partition size exceeded on part %d: %d + %d > %d * %d\n",
                     max_part, part_sizes[max_part], v_weight, g.max_partition_size, g.partition_capacities[max_part]);
              continue;
            }

            int diff_part     = 2*part_count - sum_weights;
            int diff_max_part = sum_weights - 2*max_count;
            int diff_cut      = diff_part + diff_max_part;
            double new_max_imb      = (double)(part_sizes[max_part] + v_weight) / avg_sizes[max_part];
            double new_max_edge_imb = (double)(part_edge_sizes[max_part] + out_degree) / avg_edge_sizes[max_part];
            double new_max_cut_imb  = (double)(part_cut_sizes[max_part] + diff_max_part) / avg_cut_sizes[max_part];
            double new_cut_imb      = (double)(part_cut_sizes[part] + diff_part) / avg_cut_sizes[part];

            if (new_max_imb < vert_balance &&
                new_max_edge_imb < max_e &&
                new_max_cut_imb < max_c &&
                new_cut_imb < max_c)
            {
              ++num_swapped_2;
              parts[v] = max_part;

              #pragma omp atomic
              cut_size += diff_cut;

              #pragma omp atomic
              part_cut_sizes[part]     += diff_part;
              #pragma omp atomic
              part_cut_sizes[max_part] += diff_max_part;

              #pragma omp atomic
              part_sizes[max_part] += v_weight;
              #pragma omp atomic
              part_sizes[part]    -= v_weight;

              #pragma omp atomic
              part_edge_sizes[max_part] += out_degree;
              #pragma omp atomic
              part_edge_sizes[part]     -= out_degree;

              unit_avg_cut_size = (double)cut_size / (double)g.partition_capacities_sum;
              for (int i = 0; i < num_parts; ++i)
                avg_cut_sizes[i] = unit_avg_cut_size * g.partition_capacities[i];


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
            } // end if refinement improves balance
          } // end if (max_part != part)
        } // end for loop over vertices

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
            printf("\nnum_swapped_2 (edge, refine): %d -- V: %2.2lf  E: %2.2lf  C: %2.2lf\n", num_swapped_2, vert_balance, max_e, max_c);
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

          max_e = 0.0;
          max_c = 0.0;
          for (int p = 0; p < num_parts; ++p)
          {
            if ((double)part_edge_sizes[p] / avg_edge_sizes[p] > max_e)
              max_e = (double)part_edge_sizes[p] / avg_edge_sizes[p];
            if ((double)part_cut_sizes[p] / avg_cut_sizes[p] > max_c)
              max_c = (double)part_cut_sizes[p] / avg_cut_sizes[p];
          }
          if (max_e < edge_balance)
          {
            max_e = edge_balance;
            weight_exponent_e = 1.0;
            weight_exponent_c *= max_c;
          }
          else
          {
            weight_exponent_e *= max_e / edge_balance;
            weight_exponent_c = 1.0;
          }

          #if OUTPUT_STEP
            evaluate_quality(g, num_parts, parts);
          #endif
        }
      }

      #pragma omp single
      {
        if (max_e > edge_balance*1.01 && t == edge_outer_iter-1 && num_tries < MAX_TRIES)
        {
          --t;
          if (max_e < running_max_e*0.99)
          {
            running_max_e = max_e;
            printf("Edge balance missed, attempting further iterations: (%2.3lf)\n", max_e);
          }
          else
            ++num_tries;
        }
        else
          ++t;
      }

    } // end while

    delete [] part_counts;
    delete [] part_weights;
    delete [] part_edge_weights;
    delete [] part_cut_weights;

  } // end omp parallel

  delete [] avg_sizes;
  delete [] avg_edge_sizes;
  delete [] avg_cut_sizes;

  delete [] part_sizes;
  delete [] part_edge_sizes;
  delete [] part_cut_sizes;
  delete [] queue;
  delete [] queue_next;
  delete [] in_queue;
  delete [] in_queue_next;
}




