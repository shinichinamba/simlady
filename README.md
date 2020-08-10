simlady
=======

SIMulator for Long read transcriptome Analysis with RNA DecaY model

Installation
------------

``` bash
git clone https://github.com/shinichinamba/simlady.git
cd simlady
pip3 install .
```

Usage
-----

``` bash
simlady --help
```

    usage: simlady [-h] [--version] [--fold-changes PATH] [--num-reads INT]
                   [--num-reps [INT [INT ...]]] [--gtf PATH]
                   [--gamma-params PAR PAR] [--chi2-params-s PAR PAR PAR PAR PAR]
                   [--chi2-params-n PAR PAR PAR] [--max-passes INT]
                   [--sqrt-params PAR PAR] [--norm-params PAR PAR]
                   [--probability-threshold FLOAT] [--prob-ins FLOAT]
                   [--prob-del FLOAT] [--prob-sub FLOAT] [--min-readlength INT]
                   [--lognorm-readlength [PARAMETER [PARAMETER ...]] |
                   --sample-readlength-from-fastq PATH [PATH ...] |
                   --sample-readlength-from-text PATH | --asis-fold-changes]
                   [--sam-output SAM_OUTPUT] [--no-sam]
                   [--transcriptomic-fasta-output PATH] [--write-weights]
                   [--threads INT] [--save-memory] [--seed INT] [--gzip] [--bam]
                   FASTA_PATH OUTPUT_PREFIX

    simlady v0.1.0 -- simlady: a SIMulator for Long read transcriptome Analysis
    with RNA DecaY model. Simlady uses the Pacific Biosciences SMRT error model
    which is implemented in SimLoRD.

    positional arguments:
      FASTA_PATH            A transcriptomic fasta file to sample reads from, or a
                            genomic fasta file if you specify '--gtf' flag.
                            gzipped file is allowed.
      OUTPUT_PREFIX         Save the simulated reads as a fastq-file at
                            OUTPUT_PREFIX.fastq

    optional arguments:
      -h, --help            show this help message and exit
      --version             show program's version number and exit
      --fold-changes PATH, -fc PATH
                            Read fold changes from PATH. Each line must have the
                            same number with tab-delimited format and be arranged
                            in the same row order as a reference transcript fasta,
                            empty for no fold change.
      --num-reads INT, -n INT
                            Number of reads to simulate (default= 1000).
      --num-reps [INT [INT ...]], -r [INT [INT ...]]
                            Number of repeats par groups. If you specify '--fold-
                            changes' option, number of groups must be same as the
                            column number of fold-changes file. Otherwise 1
                            (default= [1]).
      --gtf PATH, -g PATH   A transcriptomic gtf file with exonic features.
                            Transcriptomic sequences will be extracted from this
                            file and the genomic fasta file.
      --gamma-params PAR PAR, -ga PAR PAR
                            Curve determining parameters for the gamma
                            distribution: shape and scale. (default= (0.07611007,
                            771.80272892))
      --chi2-params-s PAR PAR PAR PAR PAR, -xs PAR PAR PAR PAR PAR
                            Parameters for the curve determining the parameter
                            scale for the chi^2 distribution: m,b, z, c, a for
                            'm*x + b' if x <= z and 'c * x^-a' if x > z (default=
                            (0.01214, -5.12, 675, 48303.0732881,
                            1.4691051212330266))
      --chi2-params-n PAR PAR PAR, -xn PAR PAR PAR
                            Parameters for the function determining the parameter
                            n for the chi^2 distribution: m, b, z for 'm*x + b' if
                            x < z and 'm*z + b' for x >=z (default=
                            (0.00189237136, 2.5394497, 5500)).
      --max-passes INT, -mp INT
                            Maximal number of passes for one molecule (default=
                            40).
      --sqrt-params PAR PAR, -sq PAR PAR
                            Parameters for the square root function for the
                            quality increase: a, b for 'sqrt(x+a) - b' (default=
                            (0.5, 0.2247))
      --norm-params PAR PAR, -nd PAR PAR
                            Parameters for normal distributed noise added to
                            quality increase sqare root function (default= (0,
                            0.2))
      --probability-threshold FLOAT, -t FLOAT
                            Upper bound for the modified total error probability
                            (default= 0.2)
      --prob-ins FLOAT, -pi FLOAT
                            Probability for insertions for reads with one pass.
                            (default= 0.11)
      --prob-del FLOAT, -pd FLOAT
                            Probability for deletions for reads with one pass.
                            (default= 0.04)
      --prob-sub FLOAT, -ps FLOAT
                            Probability for substitutions for reads with one pass.
                            (default= 0.01)
      --min-readlength INT, -mr INT
                            Minimum read length (default= 50) for lognormal
                            distribution
      --lognorm-readlength [PARAMETER [PARAMETER ...]], -ln [PARAMETER [PARAMETER ...]]
                            Parameters for lognormal read length distribution:
                            (sigma, loc, scale), empty for defaults. (default =
                            (0.200110276521, -10075.4363813, 17922.611306))
      --sample-readlength-from-fastq PATH [PATH ...], -sf PATH [PATH ...]
                            Sample read length from a fastq-file at PATH
                            containing reads.
      --sample-readlength-from-text PATH, -st PATH
                            Sample read length from a text file (one length per
                            line).
      --asis-fold-changes, -af
                            Do not consider read length but use fold change as is.
      --sam-output SAM_OUTPUT, -so SAM_OUTPUT
                            Save the alignments in a transcript sam-file at
                            SAM_OUTPUT. By default, use OUTPUT_PREFIX.sam. If '--
                            bam' option is specified, use OUTPUT_PREFIX.bam
      --no-sam              Do not calculate the alignment and write a sam file.
      --transcriptomic-fasta-output PATH, -to PATH
                            A path to save a transcriptomic fasta file. By
                            default, use OUTPUT_PREFIX.tx.fasta. This option is
                            ignored unless you specify --gtf option.
      --write-weights       Write out the weights matrix as a tab-delimited file.
                            The file name will be OUTPUT_PREFIX.weights.txt.
      --threads INT, -T INT
                            Number of threads. Increasing the number of threads
                            makes the processing faster at the expense of larger
                            memory requirement (default= 1)
      --save-memory         For saving memory, provoke multiprocessing only for
                            calculating weight matrix and not for generating
                            reads.
      --seed INT            Specify a seed of random generators. You should also
                            match number of threads for reproducibility.
      --gzip                Compress the simulated reads using gzip and save them
                            at OUTPUT_PREFIX.fastq.gz
      --bam                 Save the simulated sam files in BAM format.

Error model
-----------

Simlady uses the Pacific Biosciences SMRT error model which is
implemented in
[SimLoRD](https://bitbucket.org/genomeinformatics/simlord/src/master/).

Citation
--------

S Namba *et al*. Multi-sample Full-length Transcriptome Analysis of 22
Breast Cancer Clinical Specimens with Long-Read Sequencing.
***BioRxiv*** (2020)
<a href="https://doi.org/10.1101/2020.07.15.199851" class="uri">https://doi.org/10.1101/2020.07.15.199851</a>
