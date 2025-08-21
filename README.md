
# Task Mapping for EMT Simulation

 Source code related to my Master's thesis at Polytechnique Montreal.

## Algorith `PuLP`-modified:

- Algorithm `PuLP-QSAP` : Modified version of `PuLP\pulp\0.2` of `https://github.com/HPCGraphAnalysis/PuLP`

  - Please refer to original license details in `scripts\run_PuLP\PuLP-QSAP\src`
  - Original reference :

> George M. Slota, Kamesh Madduri, and Sivasankaran Rajamanickam. PuLP : Scalable multi-objective multi-constraint partitioning for small-world networks. In 2014 IEEE International Conference on Big Data (Big Data), pages 481–490, October 2014.

- Modifications, for the task mapping problem:

  - weights between labels (communication cost partitions)
  - heterogenous label capacity (partition capacity)
  - option to allow empty label/partition (for bin packing)
  - maximal load constraint on label/partition (for bin packing)
  - gnuplot script for quick visualization of partitioning
  - fix of some allocations and seg fault problems.


## Publication

> Durette, J., Karabulut Kurt, G., & Lesage-Landry, A. (2025). *Task mapping strategies for electric power system simulations on heterogeneous clusters* (Les Cahiers Du GERAD No. G-2025-32; pp. 1–12). Groupe d’études et de recherche en analyse des décisions. [https://www.gerad.ca/en/papers/G-2025-32](https://www.gerad.ca/en/papers/G-2025-32)
>
