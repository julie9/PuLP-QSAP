using namespace std;

/*
 ######
 #     # ######   ##   #####       ##   #####       #
 #     # #       #  #  #    #     #  #  #    #      #
 ######  #####  #    # #    #    #    # #    #      #
 #   #   #      ###### #    #    ###### #    #      #
 #    #  #      #    # #    #    #    # #    # #    #
 #     # ###### #    # #####     #    # #####   ####

*/
void
read_adj(char* filename, int& n, long& m,
         int*& out_array, long*& out_degree_list,
         bool has_vert_weights, bool has_edge_weights,
         int*& vertex_weights, int*& edge_weights, long& vertex_weights_sum)
{
  ifstream infile;
  string line;
  string val;

  out_array = new int[m];
  out_degree_list = new long[n+1];
  if (has_vert_weights || has_edge_weights) vertex_weights = new int[n];
  else vertex_weights = NULL;
  if (has_edge_weights || has_vert_weights) edge_weights = new int[m];
  else edge_weights = NULL;
  vertex_weights_sum = 0;

  #pragma omp parallel for
  for (int i = 0; i < n+1; ++i)
    out_degree_list[i] = 0;

  long count              = 0;
  int  cur_vert           = 0;
  bool is_intro_line_read = false;

  // =======================================================
  // Open file
  // =======================================================
  infile.open(filename);
  if (!infile.is_open())
  {
    fprintf(stderr, "Error opening file: %s\n", filename);
    abort();
  }

  // =======================================================
  // Read file
  // =======================================================
  while (getline(infile, line))
  {

    if (line[0] == '%') continue; // Skip lines starting with '%'
    // Skip the first line(n,m already read in read_graph)
    if (is_intro_line_read == false)
    {
      is_intro_line_read = true;
      continue;
    }

    stringstream ss(line);
    out_degree_list[cur_vert] = count;
    if (has_vert_weights)
    {
      getline(ss, val, ' ');
      // printf("val: %s\n", val.c_str());
      vertex_weights[cur_vert] = atoi(val.c_str());
      vertex_weights_sum += vertex_weights[cur_vert];
    }
    else if (has_edge_weights)
    {
      vertex_weights[cur_vert] = 1;
      vertex_weights_sum += vertex_weights[cur_vert];
    }
    /*else
    {
      vertex_weights[cur_vert] = rand() % 10;
      vertex_weights_sum += vertex_weights[cur_vert];
    }*/
    ++cur_vert;

    while (getline(ss, val, ' '))
    {
      // Skip empty strings
      if(!val.empty())
      {
        out_array[count] = atoi(val.c_str())-1;
        if (has_edge_weights)
        {
          getline(ss, val, ' ');
          edge_weights[count] = atoi(val.c_str());
        }
        else if (has_vert_weights)
        {
          edge_weights[count] = 1;
        }
        /*else
        {
          edge_weights[count] = rand() % 10;
        }*/
        ++count;
      }
    }
    // printf("Processed vertex %d, count: %ld\n", cur_vert, count);
  }
  out_degree_list[cur_vert] = count;
  // printf("cur_vert: %d, n: %d, count: %ld, m: %ld\n", cur_vert, n, count, m);
  assert(cur_vert == n);
  assert(count == m);

  infile.close();
}

/*
 ######
 #     # ######   ##   #####      ####  #####    ##   #####  #    #
 #     # #       #  #  #    #    #    # #    #  #  #  #    # #    #
 ######  #####  #    # #    #    #      #    # #    # #    # ######
 #   #   #      ###### #    #    #  ### #####  ###### #####  #    #
 #    #  #      #    # #    #    #    # #   #  #    # #      #    #
 #     # ###### #    # #####      ####  #    # #    # #      #    #
*/
void
read_graph(char* filename, int& n, long& m,
           int*& out_array, long*& out_degree_list,
           int*& vertex_weights, int*& edge_weights, long& vertex_weights_sum)
{
  printf("Reading graph from file %s\n", filename);

  ifstream infile;
  string line;
  int format = 0;

  infile.open(filename);

  if (!infile.is_open())
  {
    fprintf(stderr, "Error: could not open input file\n");
    abort();
  }

  // Skip comments lines starting with '%'
  while (getline(infile, line))
    if (line[0] != '%')  break;

  // printf("%s\n", line.c_str());

  sscanf(line.c_str(), "%d %li %d", &n, &m, &format);
  m *= 2;
  infile.close();

  bool has_vert_weights = false;
  bool has_edge_weights = false;
  switch(format)
  {
    case  0: break;
    case  1: has_edge_weights = true; break;
    case 10: has_vert_weights = true; break;
    case 11: has_vert_weights = true; has_edge_weights = true; break;
    default:
      fprintf (stderr, "Unknown format specification: '%d'\n", format);
      abort();
  }
  // printf("Reading graph with %d vertices and %li edges\n", n, m);

  read_adj(filename, n, m, out_array, out_degree_list,
    has_vert_weights, has_edge_weights,
    vertex_weights, edge_weights, vertex_weights_sum);
}




/*
 ######
 #     # ######   ##   #####      ####    ##   #####    ##    ####  # ##### #   #
 #     # #       #  #  #    #    #    #  #  #  #    #  #  #  #    # #   #    # #
 ######  #####  #    # #    #    #      #    # #    # #    # #      #   #     #
 #   #   #      ###### #    #    #      ###### #####  ###### #      #   #     #
 #    #  #      #    # #    #    #    # #    # #      #    # #    # #   #     #
 #     # ###### #    # #####      ####  #    # #      #    #  ####  #   #     #
*/
void
read_partition_capacities(const char* filename, int num_elements,
                          int*& partition_capacities, long& partition_capacities_sum)
{
  printf("Reading capacities from file %s\n", filename);

  std::ifstream file(filename);
  if (!file.is_open())
  {
    std::cerr << "Error: Could not open file " << filename << std::endl;
    std::abort();
  }

  partition_capacities     = new int[num_elements];
  partition_capacities_sum = 0;
  int count = 0;
  std::string line;

  while (count < num_elements && std::getline(file, line))
  {
    std::istringstream iss(line);
    int capacity;
    while (count < num_elements && iss >> capacity)
    {
      partition_capacities[count++] = capacity;
      partition_capacities_sum     += capacity;
    }
  }

  if (count != num_elements)
  {
    std::cerr << "Error: The number of integers in the file is less than expected." << std::endl;
    delete[] partition_capacities;
    std::abort();
  }

  // printf("\npartition_weights = ");
  // for (int i = 0; i < num_elements; ++i)
  //     printf("%d ", partition_capacities[i]);
  // printf("\n");

  file.close();
}






/*
 ######
 #     # ######   ##   #####     # #    # ##### ###### #####        #####    ##   #####  #####    #    # ###### #  ####  #    # #####  ####
 #     # #       #  #  #    #    # ##   #   #   #      #    #       #    #  #  #  #    #   #      #    # #      # #    # #    #   #   #
 ######  #####  #    # #    #    # # #  #   #   #####  #    # ##### #    # #    # #    #   #      #    # #####  # #      ######   #    ####
 #   #   #      ###### #    #    # #  # #   #   #      #####        #####  ###### #####    #      # ## # #      # #  ### #    #   #        #
 #    #  #      #    # #    #    # #   ##   #   #      #   #        #      #    # #   #    #      ##  ## #      # #    # #    #   #   #    #
 #     # ###### #    # #####     # #    #   #   ###### #    #       #      #    # #    #   #      #    # ###### #  ####  #    #   #    ####
*/
void
read_interpartition_weights(const char* filename, int num_parts,
                            int*& interpartition_weights)
{
  printf("Reading interpartition weights from file %s\n", filename);

  ifstream file(filename);
  if (!file.is_open())
  {
    fprintf(stderr, "Error: Could not open file %s\n", filename);
    abort();
  }

  interpartition_weights = new int[num_parts * num_parts];
  string line;
  int index = 0;

  while (getline(file, line))
  {
    istringstream iss(line);
    int weight;
    while (iss >> weight)
    {
      if (index >= num_parts * num_parts)
      {
        fprintf(stderr, "Error: More weights in file than expected for num_parts * num_parts.\n");
        delete[] interpartition_weights;
        abort();
      }
      interpartition_weights[index++] = weight;
    }
  }

  if (index != num_parts * num_parts)
  {
    fprintf(stderr, "Error: The number of weights in the file does not match num_parts * num_parts.\n");
    delete[] interpartition_weights;
    abort();
  }
  file.close();

  // printf("\ninterpartition_weights = ");
  // for (int i = 0; i < num_parts * num_parts; ++i)
  //   printf("%d ", interpartition_weights[i]);
  // printf("\n");
}


/*
 ######
 #     # ######   ##   #####     #####    ##   #####  #####  ####
 #     # #       #  #  #    #    #    #  #  #  #    #   #   #
 ######  #####  #    # #    #    #    # #    # #    #   #    ####
 #   #   #      ###### #    #    #####  ###### #####    #        #
 #    #  #      #    # #    #    #      #    # #   #    #   #    #
 #     # ###### #    # #####     #      #    # #    #   #    ####
*/
void read_parts(char* filename, int num_verts, int* parts)
{
  ifstream infile;
  string line;
  infile.open(filename);

  for (int i = 0; i < num_verts; ++i)
  {
    getline(infile, line);
    parts[i] = atoi(line.c_str());
  }

  infile.close();
}

/*
 #     #
 #  #  # #####  # ##### ######    #####    ##   #####  #####  ####
 #  #  # #    # #   #   #         #    #  #  #  #    #   #   #
 #  #  # #    # #   #   #####     #    # #    # #    #   #    ####
 #  #  # #####  #   #   #         #####  ###### #####    #        #
 #  #  # #   #  #   #   #         #      #    # #   #    #   #    #
  ## ##  #    # #   #   ######    #      #    # #    #   #    ####
*/
void write_parts(char* filename, int num_verts, int* parts)
{
  ofstream outfile;
  outfile.open(filename);

  for (int i = 0; i < num_verts; ++i)
    outfile << parts[i] << endl;

  outfile.close();
}
