#!/bin/bash

# Sample of run commands

# Set the number of threads for OpenMP
export OMP_NUM_THREADS=4

l_weight_file=./data/interconnection_weight[cluster].txt
capacity_file=./data/processor_capacity[cluster].txt
graph_file=./data/task_graph.graph
result_file=./results.parts
cap_part_size=140000
num_part=46
vertex_balance=1.01

# **Bin Packing (Q-VSBPP) Homogeneous**
./pulp ${graph_file} $num_part -v ${vertex_balance}  -o ${result_file} -b -q -a $cap_part_size >> ./output.log

# **Bin Packing (Q-VSBPP) Heterogeneous**
./pulp ${graph_file} $num -v ${vertex_balance} -o ${result_file} -w ${l_weight_file} -p ${capacity_file} -b -q -a $cap_part_size >> ./${dir}/output.log

# **Assignment (BQSAP) heterogeneous**
./pulp ${graph_file} $num -v ${vertex_balance} -o ${result_file} -w ${l_weight_file} -p ${capacity_file}  >> ./output.log 2>&1

