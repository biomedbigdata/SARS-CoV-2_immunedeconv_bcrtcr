from scirpy.ir_dist import *


def compute_distance_mat(sequences, cutoff):
  # create the AlignmentDistanceCalculator
    adc = metrics.AlignmentDistanceCalculator(cutoff, n_jobs = 1)
    # compute the weighted adjacencies
    dist_mat = adc.calc_dist_mat(sequences)
    return dist_mat


def get_alignment_distance_matrices(sequences, cutoffs, output = False):
  matrices = []
  for cutoff in cutoffs:
    dist_mat = compute_distance_mat(sequences, cutoff)
    # add matrix to list
    matrices.append(dist_mat)
  
  # return the computed matrices  
  return matrices

