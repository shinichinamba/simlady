import argparse
from . import __version__, DESCRIPTION

def get_argument_parser_lady():
    """
    Return the argument parser for this application
    removed from simlord:
        --generate-reference
        --save-reference
        --uniform-chromosome-probability
        --without-ns
        --fixed-readlength
    added to simlord: 
        --gtf
        --transcriptomic-fasta-output
        --fold-changes
        --gamma-params
        --asis-fold-changes
        --num-reps
        --seed
        --threads
        --bam
    """
    parser = argparse.ArgumentParser(prog="simlady",
        description="simlady v{} -- {}".format(__version__,DESCRIPTION))
    parser.add_argument("--version", action="version", version="simlady v"+__version__)

    parser.add_argument("fasta", metavar="FASTA_PATH",
        help="A transcriptomic fasta file to sample reads from, or a genomic fasta file if you specify '--gtf' flag. gzipped file is allowed.")
        
    parser.add_argument("--fold-changes", "-fc", metavar="PATH",
        help="Read fold changes from PATH. Each line must have the same number with tab-delimited format and be arranged in the same row order as a reference transcript fasta, empty for no fold change.")

    parser.add_argument("--num-reads", "-n", metavar="INT", type=int, default=1000,
        help="Number of reads to simulate (default= %(default)s).")
        
    parser.add_argument("--num-reps", "-r", metavar="INT", nargs="*", type=int, default=[1],
        help="Number of repeats par groups. If you specify '--fold-changes' option, number of groups must be same as the column number of fold-changes file. Otherwise 1 (default= %(default)s).")
    
    parser.add_argument("--gtf", "-g", metavar = "PATH",
        help="A transcriptomic gtf file with exonic features. Transcriptomic sequences will be extracted from this file and the genomic fasta file.")
    
    parser.add_argument("--gamma-params", "-ga", metavar = "PAR", type=float, nargs=2,
        default=(0.07611007, 771.80272892),
        help="Curve determining parameters for the gamma distribution: shape and scale.  (default= %(default)s)")

    parser.add_argument("--chi2-params-s", "-xs", metavar="PAR", type=float, nargs=5,
        default=(0.01214, -5.12, 675, 48303.0732881, 1.4691051212330266),
        help="Parameters for the curve determining the parameter scale for the chi^2 distribution: m,b, z, c, a for  'm*x + b' if x <= z and 'c * x^-a' if x > z (default= %(default)s)")

    parser.add_argument("--chi2-params-n", "-xn", metavar="PAR", type=float, nargs=3,
        default=(1.89237136e-03, 2.53944970e+00, 5500),
        help="Parameters for the function determining the parameter n for the chi^2 distribution: m, b, z  for 'm*x + b' if x < z and 'm*z + b' for x >=z (default= %(default)s).")

    parser.add_argument("--max-passes", "-mp", metavar="INT", type=int, default=40,
        help="Maximal number of passes for one molecule (default= %(default)s).")

    parser.add_argument("--sqrt-params", "-sq", metavar="PAR", type=float, nargs=2,
        default=(0.5, 0.2247),
        help="Parameters for the square root function for the quality increase: a, b for 'sqrt(x+a) - b' (default= %(default)s)")

    parser.add_argument("--norm-params", "-nd", metavar="PAR", type=float, nargs=2,
        default=(0, 0.2),
        help="Parameters for normal distributed noise added to quality increase sqare root function (default= %(default)s)")

    parser.add_argument("--probability-threshold", "-t", metavar="FLOAT", type=float,
        default=0.2,
        help="Upper bound for the modified total error probability (default= %(default)s)")

    parser.add_argument("--prob-ins", "-pi", metavar="FLOAT", type=float, default=0.11,
        help="Probability for insertions for reads with one pass. (default= %(default)s)")
    parser.add_argument("--prob-del", "-pd", metavar="FLOAT", type=float, default=0.04,
        help="Probability for deletions for reads with one pass. (default= %(default)s)", )
    parser.add_argument("--prob-sub", "-ps", metavar="FLOAT", type=float, default=0.01,
        help="Probability for substitutions for reads with one pass. (default= %(default)s)")

    parser.add_argument("--min-readlength", "-mr", metavar="INT", type=int,
        help="Minimum read length (default= %(default)s) for lognormal distribution", default=50)
    
    group_len = parser.add_mutually_exclusive_group()
    group_len.add_argument("--lognorm-readlength", "-ln", metavar="PARAMETER", nargs="*",
        type=float, 
        help="Parameters for lognormal read length distribution: (sigma, loc, scale), empty for defaults. (default = (0.200110276521, -10075.4363813, 17922.611306))")
    group_len.add_argument("--sample-readlength-from-fastq", "-sf", metavar="PATH", nargs="+",
        help="Sample read length from a fastq-file at PATH containing reads.")
    group_len.add_argument("--sample-readlength-from-text", "-st", metavar="PATH",
        help="Sample read length from a text file (one length per line).")
    group_len.add_argument("--asis-fold-changes", "-af", action="store_true",
        help="Do not consider read length but use fold change as is.")

    parser.add_argument("output", metavar="OUTPUT_PREFIX",
        help="Save the simulated reads as a fastq-file at OUTPUT_PREFIX.fastq")

    parser.add_argument("--sam-output", "-so", metavar="SAM_OUTPUT",
        help="Save the alignments in a transcript sam-file at SAM_OUTPUT. "
        "By default, use OUTPUT_PREFIX.sam. If '--bam' option is specified, use OUTPUT_PREFIX.bam")
    parser.add_argument("--no-sam", action="store_true",
        help="Do not calculate the alignment and write a sam file.")
    
    parser.add_argument("--transcriptomic-fasta-output", "-to", metavar = "PATH",
        help="A path to save a transcriptomic fasta file. By default, use OUTPUT_PREFIX.tx.fasta. This option is ignored unless you specify --gtf option.")
    
    parser.add_argument("--write-weights", action="store_true", help="Write out the weights matrix as a tab-delimited file. The file name will be OUTPUT_PREFIX.weights.txt.")

    parser.add_argument("--threads", "-T", metavar="INT", type=int, default=1,
        help="Number of threads. (default= %(default)s)")
        
    parser.add_argument("--save-memory", action="store_true", help="For saving memory, provoke multiprocessing only for calculating weight matrix and not for generating reads.")
        
    parser.add_argument("--seed", metavar="INT", type=int,
        help="Specify a seed of random generators. You should also match number of threads for reproducibility.")

    parser.add_argument("--gzip", action="store_true", help="Compress the simulated reads using gzip and save them at OUTPUT_PREFIX.fastq.gz")
    
    parser.add_argument("--bam", action="store_true", help="Save the simulated sam files in BAM format.")

    return parser
