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

#include <cstdlib>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <getopt.h>
#include <string.h>
#include <omp.h>
#include <math.h>

#include "pulp.h"
#include "io.cpp"

#include "plot_gnu.c"

void print_usage_full(char** argv)
{
  printf("To run: %s [graphfile] [num parts] [options]\n\n", argv[0]);
  printf("Options:\n");
  printf("\t-v [#.#]:\n");
  printf("\t\tVertex balance constraint [default: 1.10 (10%%)]\n");
  printf("\t-e [#.#]:\n");
  printf("\t\tEdge balance constraint [default: off]\n");
  printf("\t-c:\n");
  printf("\t\tAttempt to minimize per-part cut\n");
  printf("\t-l:\n");
  printf("\t\tDo label propagation-based initialization\n");
  printf("\t-m [#]:\n");
  printf("\t\tGenerate multiple partitions [default: 1]\n");
  printf("\t-o [file]:\n");
  printf("\t\tOutput parts file [default: graphname.part.numparts]\n");
  printf("\t-i [file]:\n");
  printf("\t\tInput parts file [default: none]\n");
  printf("\t-s [seed]:\n");
  printf("\t\tSet seed integer [default: random int]\n");
  printf("\t-w [file]:\n");
  printf("\t\tInter-partition weights file\n");
  printf("\t-p [file]:\n");
  printf("\t\tPartition weights file (computational power)\n");
  exit(0);
}

/*
'##::::'##::::'###::::'####:'##::: ##:
 ###::'###:::'## ##:::. ##:: ###:: ##:
 ####'####::'##:. ##::: ##:: ####: ##:
 ## ### ##:'##:::. ##:: ##:: ## ## ##:
 ##. #: ##: #########:: ##:: ##. ####:
 ##:.:: ##: ##.... ##:: ##:: ##:. ###:
 ##:::: ##: ##:::: ##:'####: ##::. ##:
..:::::..::..:::::..::....::..::::..::
*/
int main(int argc, char** argv)
{
  setbuf(stdout, NULL);
  srand(time(0));
  if (argc < 3)
  {
    print_usage_full(argv);
    exit(0);
  }

  int   n                  = 0;
  long  m                  = 0;
  int*  out_array          = NULL;
  long* out_degree_list    = NULL;
  int*  vertex_weights     = NULL;
  long  vertex_weights_sum = 0;
  int*  edge_weights       = NULL;

  int* interpartition_weights   = NULL;
  int* partition_capacities     = NULL;
  long partition_capacities_sum = 0;

  char* graph_name      = strdup(argv[1]);
  char* num_parts_str   = strdup(argv[2]);
  int   num_parts       = atoi(num_parts_str);
  int*  parts;
  int   num_partitions  = 1;

  char parts_out[1024]; parts_out[0] = '\0';
  char parts_in[1024]; parts_in[0]   = '\0';

  char interpart_weights_file[1024] = "\0";
  char partition_weights_file[1024] = "\0";

  double vert_balance       = 1.10;
  double edge_balance       = 1.50;
  bool   do_bfs_init        = true;
  bool   do_lp_init         = false;
  bool   do_edge_balance    = false;
  bool   do_vert_balance    = false;
  bool   do_maxcut_balance  = false;
  bool   do_bin_packing     = false;
  int    max_partition_size = 0;
  bool   eval_quality       = false;
  int    pulp_seed          = rand();

  bool using_interpartition_weights = false;
  bool using_partition_capacities   = false;

  char c;
  while ((c = getopt (argc, argv, "v:e:i:o:cs:lm:qw:p:ba:")) != -1)
  {
    switch (c)
    {
      case 'v':                               // -v flag : vertex balance
        vert_balance = strtod(optarg, NULL);
        do_vert_balance = true;
        break;
      case 'e':                               // -e flag : edge balance
        edge_balance = strtod(optarg, NULL);
        do_edge_balance = true;
        break;
      case 'i':                               // -i flag : input parts file
        strcat(parts_in, optarg);
        break;
      case 'o':                               // -o flag : output parts file
        strcat(parts_out, optarg);
        break;
      case 'c':                               // -c flag : minimize per-part cut
        do_maxcut_balance = true;
        break;
      case 'm':                               // -m flag : generate multiple assignments
        num_partitions = atoi(optarg);
        break;
      case 's':                               // -s flag : set seed
        pulp_seed = atoi(optarg);
        break;
      case 'l':                               // -l flag : do label prop init
        do_lp_init = true;
        do_bfs_init = false;
        break;
      case 'q':                               // -q flag : evaluate quality
        eval_quality = true;
        break;
      case 'w':                               // -w flag : inter-partition weights file
        strcpy(interpart_weights_file, optarg);
        using_interpartition_weights = true;
        break;
      case 'p':                               // -p flag : partition capacity (ex: computational power)
        strcpy(partition_weights_file, optarg);
        using_partition_capacities = true;
        break;
      case 'b':                               // -b flag : do bin packing (minimize number of parts if possible)
        do_bin_packing = true;
        break;
      case 'a':                               // -a flag : cap the largest vertex weight to the partition capacity
        max_partition_size = atoi(optarg);
        break;
      case '?':
      {
        if (optopt == 'v' || optopt == 'e' || optopt == 'i' || optopt == 'o' || optopt == 'm' ||
            optopt == 's' || optopt == 'w' || optopt == 'p')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        print_usage_full(argv);
        break;
      }
      default:
        abort();
    }
  }

  double elt = 0.0;
  elt = timer();

  // ==========================================================
  // Read in graph data
  // ==========================================================
  printf("Reading in data for %s ... \n", graph_name);

  // TODO(julie9): Implement usage of partition weights without interpartition weights,
  //              but for now, partition weights require interpartition weights.
  //              (but the reverse is okay).
  if (using_partition_capacities && !using_interpartition_weights)
  {
      fprintf(stderr, "Error: The -p option requires the -w option to be specified as well.\n");
      print_usage_full(argv);
      exit(EXIT_FAILURE);
  }

  if (using_interpartition_weights)
    read_interpartition_weights(interpart_weights_file, num_parts,
                                interpartition_weights);

  if (using_partition_capacities)
    read_partition_capacities(partition_weights_file, num_parts,
                              partition_capacities, partition_capacities_sum);

  read_graph(graph_name, n, m, out_array, out_degree_list,
             vertex_weights, edge_weights, vertex_weights_sum);

  pulp_graph_t g = {
    .n                        = n,
    .m                        = m,
    .out_array                = out_array,
    .out_degree_list          = out_degree_list,
    .vertex_weights           = vertex_weights,
    .edge_weights             = edge_weights,
    .vertex_weights_sum       = vertex_weights_sum,
    .interpartition_weights   = interpartition_weights,
    .partition_capacities     = partition_capacities,
    .partition_capacities_sum = partition_capacities_sum,
    .max_partition_size       = max_partition_size,
    .do_bin_packing           = do_bin_packing,
  };

  elt = timer() - elt;
  printf("... Done reading input file(s): %9.6lf\n", elt);

  parts = new int[g.n];

  // .........................................................................
  for (int i = 0; i < num_partitions; ++i)
  {
    if (strlen(parts_in) != 0)
    {
      printf("Reading in parts file %s ... \n", parts_in);
      elt = timer();

      do_lp_init  = false;
      do_bfs_init = false;

      read_parts(parts_in, g.n, parts);

      elt = timer() - elt;
      printf("... Done reading parts file: %9.6lf\n", elt);
    }
    else if (do_bfs_init)
      for (int j = 0; j < g.n; ++j) parts[j] = -1;
    else
      for (int i = 0; i < g.n; ++i) parts[i] = rand() % num_parts;

    pulp_part_control_t ppc = {
      .vert_balance                 = vert_balance,
      .edge_balance                 = edge_balance,
      .do_lp_init                   = do_lp_init,
      .do_bfs_init                  = do_bfs_init,
      .do_vert_balance              = do_vert_balance,
      .do_edge_balance              = do_edge_balance,
      .do_maxcut_balance            = do_maxcut_balance,
      .using_interpartition_weights = using_interpartition_weights,
      .using_partition_capacities   = using_partition_capacities,
      .verbose_output               = true,
      .pulp_seed                    = pulp_seed
    };

    // ==========================================================
    // Run PULP
    // ==========================================================
    printf("Beginning partitioning ... \n");
    elt = timer();

    pulp_run(&g, &ppc, parts, num_parts);

    elt = timer() - elt;
    printf("Partitioning Time: %9.6lf seconds\n\n", elt);

    // ==========================================================
    // Write output
    // ==========================================================
    char temp_out[1024]; temp_out[0] = '\0';
    strcat(temp_out, parts_out);
    if (strlen(temp_out) == 0)
    {
      strcat(temp_out, graph_name);
      strcat(temp_out, ".parts.");
      strcat(temp_out, num_parts_str);
    }
    if (num_partitions > 1)
    {
      strcat(temp_out, ".");
      stringstream ss;
      ss << i;
      strcat(temp_out, ss.str().c_str());
    }
    printf("Writing parts file %s ... \n", temp_out);
    elt = timer();

    write_parts(temp_out, g.n, parts);

    elt = timer() - elt;
    printf(" Done: %9.6lf seconds\n", elt);

    if (eval_quality)
      evaluate_quality(g, num_parts, parts);


    // ==========================================================
    // Plotting
    // ==========================================================
    if (num_partitions == 1)
    {
      printf("Plotting ... \n");
      elt = timer();

      generate_gnuplot_file(g.n, num_parts, parts, g.vertex_weights);

      elt = timer() - elt;
      printf("Plotting Time: %9.6lf seconds\n\n", elt);
    }


  } // end for num_partitions

  delete [] interpartition_weights;
  delete [] partition_capacities;

  delete [] parts;
  delete [] out_array;
  delete [] out_degree_list;
  delete [] vertex_weights;
  delete [] edge_weights;

  free(graph_name);
  free(num_parts_str);

  return 0;
}
