# Codes for the [VNC Matching Challenge](https://codex.flywire.ai/app/vnc_matching_challenge)

This repo contains the Team CIBR's solution for the VNC Challenge (members: Xin Zheng, Zuoyu Zhang, and Tatsuo Okubo at [Chinese Institute for Brain Research, Beijing](https://cibr.ac.cn/)). The [competition](https://codex.flywire.ai/app/vnc_matching_challenge) was organized by the FlyWire team at Princeton University, and we thank the organizers for this opportunity. 

## Access the data
To access the data, visit the [competition website](https://codex.flywire.ai/app/vnc_matching_challenge) and download the following:
- Two `.csv` files: one for the male connectome and one for the female connectome
- The benchmark solution

After downloading, please rename the columns labeled 'To Node Id' in both connectome files to 'To Node ID'. Place the connectome files in the `data/` directory and the benchmark solution in the `sol/` directory. 

## Files in the repository
- `random_swap_auto_detect.jl`: script for edge-based random swap
- `greedy_search_random.jl`: script for greedy refinement on male nodes
- `MatchingFunctions.jl`: utility functions for matching mapped nodes and calculating scores
- `Manifest.toml`: Julia environment file
- `Project.toml`: Julia package dependencies
- `sol_5816444.csv`: contains our solution with the score 5816444

## Julia version
This code is intended to run with `Julia 1.10.7`. Ensure that your Julia environment is set up accordingly to avoid compatibility issues.

## Methods
We used two algorithms to improve upon the benchmark solution: **edge-based random swap** and **greedy refinement on male nodes**. 

### Edge-based random swap
The algorithm randomly selects a male edge and swaps the corresponding female nodes. If the swap improves the score, it is retained; otherwise, it is discarded. This process continues until the score no longer improves significantly. The algorithm then switches to selecting female edges and swapping the corresponding male nodes. We continued alternating between male and female edge selection until no further significant improvement is achieved.

### Node-based greedy refinement
The algorithm begins by randomly selecting a male node. It then evaluates all female nodes to identify the one that yields the greatest score improvement. Once the best female node is found, the male node currently paired with this female is reassigned to the original female partner of the selected male node to maintain a 1:1 mapping. If the swap results in an improved score, it is accepted; otherwise, the original solution is retained.

## Instructions to run the code

To execute the code, you can use the following command in your terminal:

```bash
julia random_swap_auto_detect.jl --init_sol <initial_solution> --seed <random_seed> --tol <tolerance> 
```

```bash
julia greedy_search_random.jl --init_sol <initial_solution> --seed <random_seed> --max_iter <maximum_iteration> 
```
**Arguments**

Edge-based random swap
- `--init_sol`: name of the initial solution file (default: 'benchmark.csv')
- `--seed`: random seed for reproducibility (default: 392)
- `--tol`: tolerance to alternate between male edges and female edges

Greedy refinement
- `--init_sol`: name of the initial solution file (default: 'benchmark.csv')
- `--seed`: random seed for reproducibility (default: 392)
- `--max_iter`: maximum number of iterations

