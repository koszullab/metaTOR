#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Function to compute matrix of core communities Hamming distances in parallel. 
For a sample with ~300K core communities (with ~250K single contig communities),
it takes ~ 10-20 min to run.
Internally, for each core community, the function takes its array of 
N iterative Louvain communities and compares it to the arrays of N iterative 
Louvain communities for each other core community. Thus, for each pair of 
core community, it calculates a pairwise Hamming distance. The resulting 
matrix (N_cores x N_cores) is saved as a sparse matrix (COO). 
"""

import numpy as np
import pandas as pd
from sklearn import metrics
from scipy import stats
from scipy import sparse
import multiprocessing
import argparse
import glob
from functools import partial

def get_distances_splitmat (comm, core_communities_contigs):
    """
    This function takes a segment of the full iterative clustering matrix 
    and computes, for each index (i.e. contig), the hamming distance to each 
    of the other indices.
    """
    x = sparse.csr_matrix(1 - metrics.pairwise_distances(core_communities_contigs, comm.values, metric = "hamming"))
    return(x)

def main():
    
    # ------- Parse arguments ------- #
    parser = argparse.ArgumentParser(
        description="Generate matrix of Hamming distances between all pairs of core communities"
    )
    parser.add_argument(
        "-c", "--communities", help="Iterative clustering matrix, obtained by pasting all Louvain clusterings"
    )
    parser.add_argument(
        "-i", "--indices", help="Contig ids of core communities (last column of partition/core_size_indices_*.txt)"
    )
    parser.add_argument(
        "-o", "--output", help="Output file to store matrix of Hamming distances"
    )
    parser.add_argument(
        "-@", "--cores", help="Number of cores to parallelize computation", 
        default=5
    )
    args = parser.parse_args()
    communities = args.communities
    indices = args.indices
    output = args.output
    cores = int(args.cores)
    # communities = "/pasteur/homes/jaseriza/scratch/Projects/test/after0h/temp/after0h_iterations_100.txt"
    # indices = "/pasteur/homes/jaseriza/scratch/Projects/test/after0h/temp/after0h_indices_100.txt"
    # output = "/pasteur/homes/jaseriza/scratch/Projects/test/after0h/partition/partition/hamming_distance_100.npz"
    # cores = 1
    
    # ------- Import iterative clustering matrix ------- #
    communities = pd.read_csv(communities, sep='\t', header = None)
    n_contigs = communities.shape[0]
    n_iter = communities.shape[1]
    contig_names = ["contig_" + str(x+1) for x in range(0, n_contigs)]
    iter_names = ["iter_" + str(x+1) for x in range(0, n_iter)]
    communities.index = contig_names
    communities.columns = iter_names

    # ------- Import core communities ------- #
    core_communities = pd.read_csv(indices, sep='\t', header = None)
    core_communities = core_communities.apply(lambda x: [int(y) for y in x[0].split(' ')], axis = 1)
    
    # ------- Create core-community-level iterative clustering matrix ------- #
    core_communities_contigs = communities.loc[np.array(core_communities.apply(lambda x: "contig_" + str(x[0]))), ]
    core_communities_contigs.index = core_communities
    
    # ------- Compute Hamming distances in the core-community-level iterative clustering matrix, in parallel ------- #
    step = 1000
    steps = np.arange(step, len(core_communities_contigs.index)+step, step)
    split_core_communities = [core_communities_contigs[(k-step):k] for k in steps]
    pool = multiprocessing.Pool(processes = cores)
    res = pool.map(partial(get_distances_splitmat, core_communities_contigs=core_communities_contigs), split_core_communities)
    
    # ------- Save temporary matrices ------- #
    # for k in np.arange(len(res)):
    #     mat = res[k]
    #     file = os.path.dirname(output) + "/." + re.sub(".npz", "_" + str(k) + ".npz", os.path.basename(output))
    #     sparse.save_npz(file, mat)
    
    # ------- hStacking resulting list of arrays ------- #
    print('Merging sub-matrices...')
    res = sparse.hstack(res)
    # res_dense = res.todense()
    pool.close()
    # df = pd.DataFrame(res_dense)
    # df.index = [str(x) for x in range(1, len(df)+1)]
    # df.columns = [str(x) for x in range(1, len(df)+1)]
    
    # ------- Exporting matrix of Hamming distances ------- #
    print('Exporting full matrix of Hamming distances...')
    # df.to_csv(output, header=[str(x) for x in range(1, len(df)+1)], index=[str(x) for x in range(1, len(df)+1)], sep='\t', mode='a')
    # sparse_matrix = scipy.sparse.load_npz(output)
    sparse.save_npz(output, res)

if __name__ == "__main__":
    main()
