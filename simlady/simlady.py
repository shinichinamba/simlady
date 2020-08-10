# -*- coding: utf-8 -*-
import sys
import datetime
import os
import tempfile
import shutil

from math import modf, ceil, floor, log10
from random import Random, randint
from scipy.stats import chi2, gamma, lognorm
import numpy as np
import dinopy as dp
import pysam

import dask.array as da # multi processing numpy
from concurrent import futures # multi processing

from simlord import simlord as sl

from . import parser as ps
from . import gtf2fasta

DEFAULT_LOGNORMAL_PARAMETERS = sl.DEFAULT_LOGNORMAL_PARAMETERS
BASES = sl.BASES  # byte values of bases
QUALITIES = sl.QUALITIES

# from sklearn.mixture import BayesianGaussianMixture
# def multi_lognorm_fit_prob_lady(length, tx_length, gamma_params):
#     """
#     Fit provided length with log baysian gaussian mixture model,
#     and predict probability of transcripts
#     """
#     gmm = BayesianGaussianMixture(n_components=12, max_iter=100)
#     loglength = np.log(np.asarray(length)).reshape(-1,1)
#     gmm.fit(loglength_2d)
#     ary_tx_len=np.asarray(tx_length)
#     ary_tx_len += gamma_params[0] * gamma_params[1] # Add expected 5'decay length. Expected mean of gamma distribution: shape * scale.
#     loglength_tx = np.log(ary_tx_len).reshape(-1,1)
#     np.exp(gmm.score_samples(loglength_tx))

def count_length(length, weights=None):
    """
    count provided length distribution
    """
    bc = np.bincount(length, weights)
    length_ary = np.nonzero(bc)[0]
    count_ary = bc[length_ary]
    return (length_ary, count_ary)

def count_length_no_deficiency(length, weights):
    """
    count provided length distribution.
    Even if counts of some lengths are 0 due to 0 weights, the result contains them. 
    """
    (length_ary, count_ary) = count_length(length, weights)
    length_ary_nodef = np.nonzero(np.bincount(length))[0]
    length_ary_def = length_ary_nodef[np.logical_not(np.in1d(length_ary_nodef, length_ary))]
    length_ary = np.append(length_ary, length_ary_def)
    count_ary = np.append(count_ary, np.zeros(len(length_ary_def), dtype = 'int'))
    count_ary = count_ary[np.argsort(length_ary)]
    length_ary.sort()
    return (length_ary, count_ary)


def lognorm_length(tx_len, ln_values):
    """
    probability mass function based on lognorm distribusion
    """
    max_len = max(tx_len) + 1
    length_ary = np.arange(max_len, dtype='float32')
    cdf_p1 = lognorm.cdf(length_ary + 1, s=ln_values[0], loc=ln_values[1], scale=ln_values[2]) 
    cdf_p0 = lognorm.cdf(length_ary, s=ln_values[0], loc=ln_values[1], scale=ln_values[2])
    dlognorm = np.array(cdf_p1 - cdf_p0, dtype='float32')
    return (length_ary, dlognorm)

def length_driven_weight(length_ary, count_ary, tx_length, fcs, gamma_params, min_readlength, threads):
    """
    calculate probabilities for which each transcript produces reads with provided length,
    multiply them with fold_change, and fit them into provided length distribution.
    """
    dgamma = gamma_pmf(gamma_params) # gamma distribution
    
    tx_length_ary = count_length(tx_length)[0]
    tx_order = tx_length_ary.searchsorted(tx_length)
    # probabilities for which each transcript produces reads with provided length
    n = len(tx_length_ary)
    with futures.ProcessPoolExecutor(max_workers=threads) as executor:
        results = executor.map(length_prob,
            [length_ary] * n, tx_length_ary, [dgamma] * n, [min_readlength] * n,
            chunksize=1000)
    length_prob_mtx = np.vstack(list(results))
    
    tx_weight_list = []
    for fc in fcs:
        tx_count_ary = count_length_no_deficiency(tx_length, fc)[1]
        da_length_prob_mtx = da.from_array(length_prob_mtx)
        da_tx_count_ary = da.from_array(tx_count_ary.astype("float32").reshape(-1, 1))
        da_count_ary = da.from_array(np.array(count_ary, dtype='float32'))
        da_fc_length_prob = da_tx_count_ary * da_length_prob_mtx
        da_sum_fclb = da_fc_length_prob.sum(axis=0)
        da_sum_fclb = da.where(da_sum_fclb == 0, np.float32('inf'), da_sum_fclb) # avoid zero devision
        da_fitted = da_fc_length_prob * (da_count_ary / da_sum_fclb) #fit into provided length distribution
        da_length_weight = da_fitted.sum(axis=1)
        length_weight = da_length_weight.compute(num_workers=threads)
        tx_weight = (length_weight[tx_order] * fc).astype("float64") # float32 may be not enough accurate for certification of total = 1
        tx_weight_list.append(tx_weight / tx_weight.sum()) # total = 1
    return np.vstack(tx_weight_list)

def gamma_pmf(gamma_params):
    """
    probability mass function based on gamma distribusion
    """
    max_decay = ceil(gamma.ppf(0.9925, a = gamma_params[0], scale = gamma_params[1]))
    cdf_p1 = gamma.cdf(np.arange(max_decay) + 1, a = gamma_params[0], scale = gamma_params[1]) 
    cdf_p0 = gamma.cdf(np.arange(max_decay), a = gamma_params[0], scale = gamma_params[1])
    dgamma = cdf_p1 - cdf_p0
    return dgamma.astype("float32") # for saving memory, weight related np objects have type 'float32'

def length_prob(length_ary, tx_len, dgamma, min_readlength):
    """
    probability for which paticular transcript produces reads with "length_ary" bases.
    """
    dgamma2 = np.append(dgamma, 0).astype("float32")
    len_d = len(dgamma)
    list_len = list(length_ary)
    idx = [int(tx_len - i) if min_readlength <= i and 0 <= tx_len - i < len_d else -1 for i in list_len]
    return dgamma2[idx]

def gamma_decay(tx_length, min_readlength, gamma_params):
    # gamma distribution -> prevents outlier
    my_decay = floor(gamma.rvs(a = gamma_params[0], scale = gamma_params[1]))
    new_readlength = tx_length - my_decay
    while my_decay > gamma.ppf(0.9925, a = gamma_params[0], scale = gamma_params[1]) or new_readlength < min_readlength:
        my_decay = floor(gamma.rvs(a = gamma_params[0], scale = gamma_params[1]))
        new_readlength = tx_length - my_decay
    return (new_readlength, my_decay)

def modify_cum_probs(cum_probs):
    """
    Insersion prolongs read length.
    Therefore, insersion, deletion, and mutation rate must be modified based on expected genarated read length.
    E(generated_readlength) = original_readlength * (1 + ci + ci^2 + ci^3 + ... + ci^inf) = original_readlength / (1 - ci)
    E(generated_ci) = (ci / (1 - ci)) / (1 / (1 - ci)) = ci
    E(generated_cd) = cd * (1 - ci)
    E(generated_cs) = cs * (1 - ci)
    """
    ci, cd, cs = cum_probs
    cs = cs - cd
    cd = cd - ci
    new_cs = cs / (1 - ci)
    new_cd = cd / (1 - ci)
    return ([ci, ci + new_cd, ci + new_cd + new_cs])

def read_fc(fold_change):
    """
    Read a fold_change file. fold_change file must have 1 numeric value in each line.
    """
    f = open(fold_change, 'r')
    chrlist = f.readlines()
    f.close()
    return np.float32(chrlist)
    
def read_fcs(fold_changes):
    """
    Read a fold_changes file. fold_changes file must be tab-delimited format and have the same number of values in each line.
    """
    return np.loadtxt(fold_changes, dtype="float32", delimiter="\t", converters=None, usecols=None, unpack=True, ndmin=2)


def sample_reads_lady(reference, num_reads, min_readlength,
        weights, gamma_params, chi2_params_n, chi2_params_s, max_passes,
        sqrt_params, norm_params, min_exp, prob_ins, prob_del, prob_subst,
        output_path, sam_output, no_sam, seed, thread_idx, sample_idx):
    """
    Sample 'num_reads' reads from the given 'reference' and write them to a file.
    For each read, 
        determine the read length, pick up a transcript, 
        remove 5' bases based on gamma distribution, 
        and draw a number of passes,
        depending on the length and the parameters of the chi^2 distribution.
        Adjust the base error probabilities according to the number of passes.
        Add errors to the read.
        Write the read at 'output_path' and  if 'no_sam' is False, the alignment
        at 'sam_output'.
    """
    if seed is None:
        seed = randint(0, 100000)
    seed += thread_idx * 1000 # avoid the same seeds
    seed += sample_idx * 111
    nprg = np.random.RandomState(seed)
    rg = Random(seed)

    (reference_names, reference_seqs, reference_lengths, n_ref) = reference
    
    # normal distributed noise for quality increase has to be adapted with a sigmoidal factor
    sigmoidal_factor = lambda x: 1 / (1 + 2**(-2.5/3*x+6.5/3))

    if not no_sam:
        sam_writer = pysam.AlignmentFile(sam_output, 'wh',
            reference_names=reference_names, 
            reference_lengths=reference_lengths)
    # with open(output_path, 'wb') as f:
    with sl.conditional_open(output_path) as f:
        # sample the reads
        for i in range(int(num_reads)):
            # choose read length and 
            length = -1

            # each tx has the same probability with fold change as weights
            index = nprg.choice(n_ref, p=weights)
            current_ref = reference_seqs[index]
            tx_name = reference_names[index]
            tx_length = reference_lengths[index]

            # get a random read from the chosen reference
            (current_readlength, start_pos) = gamma_decay(tx_length, min_readlength, gamma_params)
            read = get_raw_read_lady(bytearray(current_ref), start_pos, rg)
            
            # compute number of passes in prefix and suffix of the read
            (current_passes, cut_position, passes_left, passes_right, percentage_left) \
                = calculate_passes_lady(current_readlength, chi2_params_n, chi2_params_s, max_passes, rg)
            # modify 1-pass error probabilities according to number of passes
            (cum_probs_left, cum_probs_right) = sl.modify_probabilities(
                min_exp, sqrt_params, norm_params, sigmoidal_factor,
                passes_left, passes_right,
                prob_ins, prob_del, prob_subst)
            # modify cum_probs for read elongation by insertion
            cum_probs_left = modify_cum_probs(cum_probs_left)
            cum_probs_right = modify_cum_probs(cum_probs_right)
            # insert errors into read
            (read, insertions, deletions, substitutions, length_left, length_right, start_pos) = traverse_read_lady(
                read, current_ref, start_pos, cut_position, cum_probs_left, cum_probs_right, passes_left, nprg, rg)
            # update readlength
            generated_readlength = length_left + length_right
            # update percentage_left
            percentage_left = length_left / generated_readlength
            ### read = read[0:current_readlength]  # crop to readlength. no need for Isoseq
            # determine quality characters
            quality_left = QUALITIES[min(round(-10 * log10(cum_probs_left[2])), len(QUALITIES)-1)]
            quality_right = QUALITIES[min(round(-10 * log10(cum_probs_right[2])), len(QUALITIES)-1)]
            quality_values = quality_left * length_left + quality_right * length_right 
            # write the read
            num_total_errors = len(insertions) + len(deletions) + len(substitutions)
            total_error_prob = cum_probs_left[2]*percentage_left \
                               + cum_probs_right[2]*(1.0-percentage_left)
            name = ";".join(["Read="+str(i), "|read/f1p0", "length="+str(generated_readlength)+"bp",
                             "startpos="+str(start_pos),
                             "transcript="+tx_name,
                             "numberOfErrors="+str(num_total_errors),
                             "totalErrorProb={:.2f}".format(total_error_prob),
                             "passes="+str(current_passes), "passesLeft="+str(passes_left),
                             "passesRight="+str(passes_right), "cutPosition="+str(length_left),
                            ])
            # Create a mutable byte string then make a single write operation
            out_byte_str = bytearray()
            out_byte_str += "@{}\n".format(name).encode()
            # convert to string and reverse the read, maybe
            read = read.decode("ascii")
            out_byte_str += "{0}\n+\n{1}\n".format(read, quality_values).encode() # reversals is always false
            # write to the out_path in bytes, this method is faster than string writing when making large read sets
            f.write(out_byte_str)

            # write SAM file if desired
            if not no_sam:
                mapping_error_probability = 0.0000000001 # very small error probability, since we know the alignment
                write_sam_file_lady(sam_writer, name, read, quality_values.encode(),
                                    mapping_error_probability, tx_name, start_pos,
                                    current_readlength, insertions, deletions, substitutions)
            pass  # end for i
    if not no_sam:
        sam_writer.close()

def write_sam_file_lady(sam_writer, name, read, quality_values,
        mapping_error_probability, tx_name, start_pos, current_readlength,
        insertions, deletions, substitutions):
    """
    Write the sam file entry for a given read into an open sam_writer.
    name: name of the read (str)
    read: sequence of the read (bytearray)
    quality_values: (bytes)
    mapping_quality: error probability for correct mapping
    tx_name: (str)
    start_pos: start position of the read in the chromosome (int)
    current_readlength: (int)
    insertions: list of positions with insertions
    deletions: list of positions with deletions
    substitutions: list of positions with substitutions
    """
    assert type(quality_values) == bytes
    cigar = calculate_cigar_operations_lady(current_readlength, insertions, deletions, substitutions)

    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = read
    a.flag = 0 #always not reverse
    a.reference_id = sam_writer.gettid(tx_name)
    a.reference_start = start_pos
    a.mapping_quality = -10*log10(mapping_error_probability)
    a.cigar = cigar
    a.next_reference_id = -1
    a.next_reference_start = -1
    a.template_length = current_readlength
    a.query_qualities = pysam.qualitystring_to_array(quality_values)
    number_of_errors = len(insertions) + len(deletions) + len(substitutions)
    a.set_tag("NM", number_of_errors, "i")
    a.set_tag("XS", 0, "i")
    a.set_tag("XE", current_readlength, "i")
    a.set_tag("XQ", current_readlength, "i")
    sam_writer.write(a)

def get_raw_read_lady(current_ref, start_pos, random_generator, bases=BASES):
    """
    Given a reference 'current_ref' as bytearray,
    copy the whole read (omit 5'decayed bases).
    In the read, substitute all Ns with random bases.

    Return read as a bytearray.
    ADDITIONAL: reads must contain 3'. 
    start_pos is generated by 5'decay function.
    """
    assert type(current_ref) == bytearray
    read = current_ref[start_pos:]  # bytearray
    
    Ns = frozenset([ord('N'), ord('n')])
    # change Ns in reference to a random base in read
    for b, base in enumerate(read):
        if base in Ns:
            read[b] = random_generator.choice(bases)
    return (read)

def calculate_passes_lady(current_readlength, chi2_params_n, chi2_params_s, max_passes, random_generator):
    """
    Calculate the parameter n and s for the chi^2 distribution based on the
    'current_readlength' (int), 'chi2_params_n' (3-tuple float) and 'chi2_params_s'
    (5-tuple float) and draw the number of passes.
    The number of passes has max_passes (int) as upper bound.
    
    Calculate the cut position and the number of passes
    for the right and left side of the cut.
    Passes is leftside right against DNA, since Iso-seq read transcripts from 3'.

    Return (current_passes, cut_position, passes_left, passes_right, percentage_left):
    current_passes: passes for the current read (float)
    cut_position: position in the read where the qualities are split (int)
    passes_left: passes for the left part of the read (int)
    passes_right: passes for the right part of the read (int)
    percentage_left: percentage of the read on the left side of the cut (float)
    """
    # compute n parameter
    my_n = chi2_params_n[0] * min(current_readlength, chi2_params_n[2]) + chi2_params_n[1]
    current_n = max(0.001, my_n)

    # compute scale parameter
    if current_readlength <= chi2_params_s[2]:
        my_scale = chi2_params_s[0] * current_readlength - chi2_params_s[1]
        current_scale = max(0.001, my_scale)
    else:
        current_scale = chi2_params_s[3] / (current_readlength**chi2_params_s[4])

    # draw passes
    my_passes = chi2.rvs(current_n, scale=current_scale, loc=1)
    # draw new, when my_passes is greater than the 99.25 quantile of the current
    # chi^2 distribution -> prevents outlier over the a/x boundary
    while my_passes > chi2.ppf(0.9925, current_n, loc=1, scale=current_scale):
        my_passes = chi2.rvs(current_n, scale=current_scale, loc=1)
    # cut at maximum passes -> only for very small reads important
    current_passes = min(max_passes, my_passes)
    # find cut point
    (cut_factor, whole_passes) = modf(current_passes) # percentage of the read on one side of the cut, number of whole passes
    # right side of read has higher quality
    ### This is chosen at random.
    higher_side_flag = random_generator.random() < 0.5 # one of the fastest way to generate random boolean in python
    if higher_side_flag:
        cut_position = round(current_readlength * (1-cut_factor))
        passes_left = floor(current_passes)
        passes_right = ceil(current_passes)
        percentage_left = 1 - cut_factor
    # left side of read has higher quality
    else:
        cut_position = round(current_readlength*cut_factor)
        passes_left = ceil(current_passes)
        passes_right = floor(current_passes)
        percentage_left = cut_factor
    return (current_passes, cut_position, passes_left, passes_right, percentage_left)

def traverse_read_lady(read, current_ref, start_pos, cut_position, cum_probs_left, cum_probs_right, passes_left, np_random_generator, random_generator,
    bases=BASES):
    """
    Traverse the read and insert errors with the given probabilities.
    read: bytearray with ASCII codes of ACGT
    current_ref: bytearray of chromosome
    start_pos: start position of read in chromosome (current_ref)
    cut_position: use cum_probs_left in read[0:cut_position] 
        and cum_probs_right in read[cut_position:]
    cum_probs_left: triples of cumulative error probabilities (ins, ins+del, ins+del+subst)
    cum_probs_right: dito
    passes_left: ADDED. passes for the left part of the read (int). if passes_left == 0, left reads are removed.

    Return the modified read, the lists containing the errors, read length, and updated start_pos:
    (read, insertions, deletions, substitutions, length_left, length_right, start_pos)
    """
    insertions = []
    deletions = []
    substitutions = []
    current_readlength = len(read)
    Ns = frozenset([ord('N'), ord('n')])
    length_left = cut_position
    length_right = current_readlength - cut_position

    random_numbers = np_random_generator.random(current_readlength*2).tolist()  # *2 for possible insertions
    j = 0  # current read position
    m = 0  # position modifier. It is required for dynamic read update. equal to len(ins) - len(del)
    current_rand_pos = 0
    while j < current_readlength:
        assert 0 <= j < current_readlength, "error! j not in correct range"
        ci, cd, cs = cum_probs_left  if j < cut_position  else cum_probs_right
        r = random_numbers[current_rand_pos]
        current_rand_pos += 1
        jm = j + m

        if r < ci:
            # insert a base
            insertions.append(j)
            read[jm:] = bytearray([random_generator.choice(bases)]) + read[jm:]  # NOT preserve read length
            m += 1
            j -= 1 # for multiple insertion
            if j < cut_position:
                length_left += 1
            else:
                length_right += 1
        elif r < cd:
            # delete a base
            deletions.append(j)
            read = read[0:jm] + read[jm+1:]
            m -= 1
            if j < cut_position:
                length_left -= 1
            else:
                length_right -= 1
        elif r < cs:
            # substitute a base
            substitutions.append(j)
            base = random_generator.choice(bases)
            while base == read[jm]:
                base = random_generator.choice(bases)
            read[jm] = base
        #else:
            # no error at this position
        j += 1
    return (read, insertions, deletions, substitutions, length_left, length_right, start_pos)

def calculate_cigar_operations_lady(current_readlength, insertions, deletions, substitutions):
    """
    Given a read length, and three lists of positions
    with insertions, deletions, substitutions, respectively,
    calculate the cigar string for one read.
    Return list of pairs of (cigar operation codes, count) for pysam.
    Updated because current_readlength is not equal to generated read length any more.
    """
    MATCH, DELETION, INSERTION, SUBST = (7, 2, 1, 8)  # PySam CIGAR Operation Codes
    cigar = []
    count = 0
    last_op = MATCH
    point_ins = 0
    point_del = 0
    point_sub = 0
    for i in range(current_readlength):
        if point_ins < len(insertions) and i == insertions[point_ins]:
            # multiple insertions get the same index
            cigar.append((last_op, count))
            count = 1
            last_op = INSERTION
            point_ins += 1
            while point_ins < len(insertions) and i == insertions[point_ins]:
                count += 1
                point_ins += 1
                
        if point_del < len(deletions) and i == deletions[point_del]:
            # del
            point_del += 1
            if last_op == DELETION:
                count += 1
            else:
                cigar.append((last_op, count))
                count = 1
                last_op = DELETION
        elif point_sub < len(substitutions) and i == substitutions[point_sub]:
            point_sub += 1
            if last_op == SUBST:
                count += 1
            else:
                cigar.append((last_op, count))
                count = 1
                last_op = SUBST
        else:
            if last_op == MATCH:
                count += 1
            else:
                cigar.append((last_op, count))
                count = 1
                last_op = MATCH
    cigar.append((last_op, count))
    if cigar[0][1] == 0:
        cigar = cigar[1:]
    return cigar

def simulate_lady(args):
    """
    Simulate reads according to the given parameters.
    """
    t1 = datetime.datetime.now()
    # Obtain the reference.
    if args.gtf is None:
        reference = gtf2fasta.read_reference_lady(args.fasta)
    else:
        (tx_id, tx_seq) = gtf2fasta.gtf2fasta(args.fasta, args.gtf)
        print("extract transcriptomic fasta")
        if args.transcriptomic_fasta_output is not None:
            tx_fa_out = args.transcriptomic_fasta_output
        else:
            if args.output != "-":
                tx_fa_out = args.output + ".tx.fasta"
            else:
                tx_fa_out = "reads.tx.fasta"
        gtf2fasta.write_fasta(tx_id, tx_seq, tx_fa_out)
        del tx_id
        del tx_seq
        print("write out transcriptomic fasta")
        reference = gtf2fasta.read_reference_lady(tx_fa_out)
    t2 = datetime.datetime.now()
    print("Time for reading the reference: {} h".format(t2-t1), file=sys.stderr)

    num_reads = args.num_reads
    
    if args.fold_changes is not None:
        fold_changes = read_fcs(args.fold_changes)
        if fold_changes.shape[1] != reference[3]:
            raise ValueError("Fold change file must have the same row number as the number of reference transcripts.")
    else:
        fold_changes = np.ones(shape=(1,reference[3]), dtype='float32')
    
    if len(args.num_reps) != fold_changes.shape[0]:
        raise ValueError("Length of num-reps must be same as the column number of fold-changes file.")
    
    if args.prob_ins + args.prob_del + args.prob_sub > 1:
        raise ValueError("Sum of error probabilities must be < 1.")
    if args.probability_threshold <= 0.0 or args.probability_threshold >= 1.0:
        raise ValueError("Probability threshold t must be between 0 and 1 (0 < t < 1).")

    min_exp = sl.calculate_minimum_exponent(args.probability_threshold, args.prob_ins,
                                         args.prob_del, args.prob_sub)


    # Manage the read length sampling method.
    t3 = datetime.datetime.now()
    ln_values = DEFAULT_LOGNORMAL_PARAMETERS

    if args.sample_readlength_from_fastq is not None:
        lengths = sl.read_readlenghts_from_reads(args.sample_readlength_from_fastq, None)
        (length_ary, count_ary) = count_length(lengths)
    
    elif args.sample_readlength_from_text is not None:
        lengths = sl.read_readlenghts_from_file(args.sample_readlength_from_text)
        (length_ary, count_ary) = count_length(lengths)
        
    elif args.lognorm_readlength is not None:
        nparams = len(args.lognorm_readlength)
        if nparams == 0:  # -ln without parameter
            ln_values = DEFAULT_LOGNORMAL_PARAMETERS
        elif nparams == 3:
            ln_values = tuple(args.lognorm_readlength)
        else:
            raise ValueError("Wrong number of parameters for lognorm distribution. Values for sigma, loc and scale are required.")
        (length_ary, count_ary) = lognorm_length(reference[2], ln_values)
    
    elif args.asis_fold_changes:
        weights = fold_changes / np.sum(fold_changes, axis=1).reshape(-1,1)
    
    else:  # all read length args are None, use default of lognormal distribution
        ln_values = DEFAULT_LOGNORMAL_PARAMETERS
        (length_ary, count_ary) = lognorm_length(reference[2], ln_values)
    
    if not args.asis_fold_changes:
        weights = length_driven_weight(length_ary, count_ary, reference[2], fold_changes, args.gamma_params, args.min_readlength, args.threads)
    
    weights = weights.repeat(args.num_reps, axis=0) # repeat weights at num_reps times
    
    if args.write_weights:
        np.savetxt(args.output + '.weights.txt', weights.T, delimiter='\t')
    t4 = datetime.datetime.now()
    print("Time for calculation of weights: {} h.".format(t4-t3), file=sys.stderr)
    
    num_sample = weights.shape[0]
    if num_sample == 1:
        output_num = [""]
    else:
        output_num = ["_" + str(i).zfill(2) for i in range(1, num_sample + 1)]
    
    if args.bam:
        sam_ext = ".bam"
        sam_format = "BAM"
    else:
        sam_ext = ".sam"
        sam_format = "SAM"
    
    if args.output != "-":
        if args.gzip:
            fq_suffix = ".fastq.gz"
        else:
            fq_suffix = ".fastq"
        output = [args.output + n + fq_suffix for n in output_num]
        if args.sam_output is None:
            sam_output = [args.output + n + sam_ext for n in output_num]
    else:
        output = ["reads" + n + ".fastq" for n in output_num]
        if args.sam_output is None:
            sam_output = ["reads" + n + sam_ext for n in output_num]
            tempfile.mkstemp()
    # Actually sample the reads
    t5 = datetime.datetime.now()
    if args.threads == 1 or args.save_memory:
            (tmpfd, tmpsam) = tempfile.mkstemp()
    else:
        # touch tempfiles
        tmpfq = []
        tmpsam = []
        for i in range(args.threads):
            (tmpfd, tmppath) = tempfile.mkstemp(suffix=fq_suffix)
            # outsock = os.fdopen(tmpfd,'w')
            # outsock.close()
            tmpfq.append(tmppath)
            (tmpfd, tmppath) = tempfile.mkstemp() 
            # outsock = os.fdopen(tmpfd,'w')
            # outsock.close()
            tmpsam.append(tmppath)
        quotient = num_reads // args.threads
        remainder = num_reads % args.threads
        num_reads_threads = [quotient + 1 if i < remainder else quotient for i in range(args.threads)]

    for w, o, s, k in zip(weights, output, sam_output, range(num_sample)):
        if args.threads == 1 or args.save_memory:
            sample_reads_lady(
                reference=reference, 
                num_reads=num_reads, 
                min_readlength=args.min_readlength,
                weights=w,
                gamma_params=args.gamma_params,
                chi2_params_n=args.chi2_params_n, chi2_params_s=args.chi2_params_s,
                max_passes=args.max_passes, sqrt_params=args.sqrt_params,
                norm_params=args.norm_params, min_exp=min_exp, prob_ins=args.prob_ins,
                prob_del=args.prob_del, prob_subst=args.prob_sub,
                output_path=o, sam_output=tmpsam,
                no_sam=args.no_sam,
                seed=args.seed, thread_idx=1, sample_idx=k
            )
            if not args.no_sam:
                pysam.merge("-f", "-O", sam_format, "-@", str(args.threads), s, tmpsam)
        else:
            # multi threading
            future_list = []
            with futures.ProcessPoolExecutor() as executor:
                for i in range(args.threads):
                    future = executor.submit(
                        fn=sample_reads_lady, 
                        reference=reference, 
                        num_reads=num_reads_threads[i], 
                        min_readlength=args.min_readlength,
                        weights=w,
                        gamma_params=args.gamma_params,
                        chi2_params_n=args.chi2_params_n, chi2_params_s=args.chi2_params_s,
                        max_passes=args.max_passes, sqrt_params=args.sqrt_params,
                        norm_params=args.norm_params, min_exp=min_exp, prob_ins=args.prob_ins,
                        prob_del=args.prob_del, prob_subst=args.prob_sub,
                        output_path=tmpfq[i], sam_output=tmpsam[i],
                        no_sam=args.no_sam,
                        seed=args.seed, thread_idx=i, sample_idx=k)
                    future_list.append(future)
                _ = futures.as_completed(fs=future_list)
            # concatenate results
            with open(o, 'wb') as f:
                for indifq in tmpfq:
                    with open(indifq, 'rb') as fq:
                        shutil.copyfileobj(fq, f)
            if not args.no_sam:
                pysam.merge("-f", "-O", sam_format, "-@", str(args.threads), s, *tmpsam)
    # remove tempfiles
    if args.threads == 1 or args.save_memory:
        os.remove(tmpsam)
    else:
        for f in tmpfq:
            os.remove(f)
        for f in tmpsam:
            os.remove(f)
    t6 = datetime.datetime.now()
    print("Time for simulation of {} reads for {} samples: {} h.".format(num_reads, num_sample, t6-t5),
        file=sys.stderr)
        
def main():
    """
    (main function)
    Run simulate_lady() with given arguments.
    """
    parser = ps.get_argument_parser_lady()
    args = parser.parse_args()
    simulate_lady(args)

if __name__ == "__main__":
    main()
