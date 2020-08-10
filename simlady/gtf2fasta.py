from gtfparse import read_gtf
import dinopy as dp
from dinopy.conversion import string_to_bytes
from dinopy.processors import reverse_complement
from pyfaidx import Fasta

def gtf2fasta(genome_fa_path, gtf_path):
    """
    read exonic gtf and reference genome fasta, and obtain transcript id and sequence
    """
    # read
    fa = Fasta(genome_fa_path, read_ahead=1000)
    df = read_gtf(gtf_path)
    # filter: exon
    df_exon = df[df["feature"] == "exon"]
    # get exon seq (time bottleneck) # TODO get seq only for redundant exon 
    seq = [string_to_bytes(fa.get_seq(c, s, e).seq) for c, s, e in zip(df_exon['seqname'].values.tolist(), df_exon['start'].values.tolist(), df_exon['end'].values.tolist())]
    # summarise key 
    key=df_exon['transcript_id'].values.tolist()
    key_lag = ['']
    key_lag.extend(key[:-1])
    tx_change = [k != kl for k, kl in zip(key, key_lag)]
    tx_change_idx = [i for i, val in enumerate(tx_change) if val]
    tx_change_idx.append(len(key) + 1)
    # summarise tx_id
    tx_id = [key[i] for i in tx_change_idx[:-1]]
    # summarise tx-level seq
    tx_seq_fwd = [b''.join(seq[s:e]) for s, e in zip(tx_change_idx[:-1], tx_change_idx[1:])]
    # summarise tx-level strand
    exon_strand = df_exon['strand'].values.tolist()
    tx_strand = [exon_strand[i] for i in tx_change_idx[:-1]]
    # reverse complement
    tx_seq = [seq if strand == '+' else reverse_complement(seq) for seq, strand in zip(tx_seq_fwd, tx_strand)]
    return (tx_id, tx_seq)

def write_fasta(tx_id, tx_seq, fa_out):
    """
    write out fasta file.
    """
    with open(fa_out, mode='w') as f:
        for txid, txseq in zip(tx_id, tx_seq):
            f.write(">" + txid + "\n" + txseq.decode("utf-8") + "\n")
    
def read_reference_lady(input_path):
    """
    Read a reference tx in fasta format by using pyfaidx.
    Return the reference as a list of tuples (pyfaidx_fasta_object, number_of_tx, reference_lengths),
    Reference_lengths has the same order with pyfaidx_fasta_object.keys().
    """
    fp = dp.FastaReader(input_path)
    reference_lengths, reference_names, reference_seqs = [], [], []
    for sequence, name, length, interval in fp.entries(dtype=bytearray):
        name = name.decode("utf-8").replace(" ", "_")
        reference_names.append(name)
        reference_lengths.append(length)
        reference_seqs.append(sequence)
    n_ref = len(reference_lengths)
    return (reference_names, reference_seqs, reference_lengths, n_ref)

# def fastarecord2bytes(fasta_record):
#     return string_to_bytes(fasta_record[:].seq)

# def fasta2ref(tx_id, tx_seq):
#     """
#     input: transcript id and sequence
#     output: same as read_reference() in simlord
#     """
#     transcripts, weights, reference_names, reference_lengths = [], [], [], []
#     max_chrom_length = 0
#     for name, seq in zip(tx_id, tx_seq):
#         name = name.replace(" ", "_")
#         length = len(seq)
#         reference_names.append(name)
#         reference_lengths.append(length)
#         if length > 1:
#             transcripts.append((seq, name, length, 0))
#             weights.append(length)
#             if length > max_chrom_length:
#                 max_chrom_length = length
#     sum_w = sum(weights)
#     weights = [x/sum_w for x in weights]
#     return (transcripts, reference_names, reference_lengths, max_chrom_length, weights)
