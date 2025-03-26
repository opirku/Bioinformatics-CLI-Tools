"""
Author: Orion Pirku
Date: 26.03.2025
Affiliation: Master Student/Research Assistant, 
Max Planck Institute for Evolutionary Biology, 
Ploen, Germany

 
"""
import sys
import argparse
from turtle import distance
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
from numba import njit
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np
import pandas as pd
def read_fasta_file(file_path: str) -> list[SeqRecord]:
    """
    Reads a FASTA amino acids file and returns a list of SeqRecord objects.

    :param file_path: Path to the FASTA file
    :return: List of SeqRecord objects
    """
    return list(SeqIO.parse(file_path, format="fasta"))


def find_SY_indeces(protein_records: list) -> list[int]:
        """
        Finds the index of the first occuring PY peptide meaning 2nd Zinc Finger
        :param protein_records: Input parsed .faa file
        :return: List of indeces where the first PY peptide occurs in each sequence 
        """    
        indeces: list = []
        
        for record in protein_records:
        
                index = record.seq.find("PY")
                
                if index != -1:
        
                        indeces.append(index)                
        
        return indeces

def trim_sequences(protein_records: list[SeqRecord], indeces: list[int]) -> list[SeqRecord]:
        
        """
        Trims protein sequences by removing any sequences before the first occurance of PY
        :param protein_records and indeces
        :return trimmed_sequences
        """
        
        trimmed_sequences: list[SeqRecord] = []
        
        for i, record in enumerate(protein_records):
                
                start_index = indeces[i] if i < len(indeces) else 0
                
                trimmed_sequence = record.seq[start_index:]
                
                trimmed_record = SeqRecord(id=record.id, seq=trimmed_sequence)
                
                trimmed_sequences.append(trimmed_record)
               
        return trimmed_sequences

def merge_identical_sequences(input_fasta: str, output_fasta: str):
    """Merges identical sequences in a FASTA file by combining their headers.

    Args:
        input_fasta (str): Path to input FASTA file.
        output_fasta (str): Path to output FASTA file.
    """
    sequence_dict = defaultdict(list)

    # Read sequences and group by sequence
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_str = str(record.seq)
        sequence_dict[seq_str].append(record.id)

    # Write unique sequences with merged headers
    with open(output_fasta, "w") as out_f:
        for seq, ids in sequence_dict.items():
            merged_header = ";".join(ids)  # Join headers with ";"
            out_f.write(f">{merged_header}\n{seq}\n") 
    
    return sequence_dict

@njit
def compute_edit_distance_njit(seq1_arr: np.ndarray, seq2_arr: np.ndarray, substitution_penalty: int, deletion_penalty: int, insertion_penalty: int) -> int:
    n = seq1_arr.shape[0]
    m = seq2_arr.shape[0]
    dp = np.zeros((n + 1, m + 1), dtype=np.int64)
    
    # Initialize DP table
    for i in range(n + 1):
        dp[i, 0] = i * deletion_penalty
    for j in range(m + 1):
        dp[0, j] = j * insertion_penalty
    
    # Compute the DP table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1_arr[i - 1] == seq2_arr[j - 1]:
                cost = 0
            else:
                cost = substitution_penalty
            dp[i, j] = min(dp[i - 1, j] + deletion_penalty,
                           dp[i, j - 1] + insertion_penalty,
                           dp[i - 1, j - 1] + cost)
    return dp[n, m]


def compute_edit_distance(seq1: str, seq2: str, substitution_penalty: int, deletion_penalty: int, insertion_penalty: int) -> int:
    """
    Converts sequences to numpy arrays of ASCII codes and computes the edit distance.
    """
    seq1_arr = np.frombuffer(seq1.encode('ascii'), dtype=np.uint8)
    seq2_arr = np.frombuffer(seq2.encode('ascii'), dtype=np.uint8)
    return compute_edit_distance_njit(seq1_arr, seq2_arr, substitution_penalty, deletion_penalty, insertion_penalty)

def compute_edit_distance_table(fasta_file: str, substitution_penalty: int, deletion_penalty: int, insertion_penalty: int) -> pd.DataFrame:
    """
    Computes a table of Levenshtein distances between all sequences in a FASTA file.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    seqs = [str(record.seq) for record in records]
    ids = [record.id for record in records]

    distance_data = []
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            dist = compute_edit_distance(seqs[i], seqs[j], substitution_penalty, deletion_penalty, insertion_penalty)
            distance_data.append({"Sequence_1_ID": ids[i], "Sequence_2_ID": ids[j], "Edit_Distance": dist})
    return pd.DataFrame(distance_data)

def cluster_sequences(edit_distance_df: pd.DataFrame, method: str = "average", threshold: float = None) -> pd.DataFrame:
    """
    Clusters sequences based on edit distances using hierarchical clustering and plots a pretty dendrogram.
    Also saves the dendrogram as a PNG file and the clusters as a CSV.
    """
    # Extract unique sequence IDs
    unique_ids = sorted(set(edit_distance_df["Sequence_1_ID"]) | set(edit_distance_df["Sequence_2_ID"]))

    # Create a distance matrix
    distance_matrix = pd.DataFrame(np.zeros((len(unique_ids), len(unique_ids))), index=unique_ids, columns=unique_ids)
    for _, row in edit_distance_df.iterrows():
        seq1, seq2, dist = row["Sequence_1_ID"], row["Sequence_2_ID"], row["Edit_Distance"]
        distance_matrix.loc[seq1, seq2] = dist
        distance_matrix.loc[seq2, seq1] = dist  # Symmetric

    # Convert distance matrix to condensed form
    condensed_dist_matrix = sch.distance.squareform(distance_matrix)

    # Perform hierarchical clustering
    linkage_matrix = sch.linkage(condensed_dist_matrix, method=method)

    # Set a color palette for clusters
    sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(12, 6))
    dendro = sch.dendrogram(
        linkage_matrix,
        labels=unique_ids,
        leaf_rotation=90,
        leaf_font_size=10,
        above_threshold_color="gray",
        color_threshold=threshold,
        ax=ax
    )

    ax.set_title("Hierarchical Clustering of Sequences", fontsize=14, fontweight="bold", pad=15)
    ax.set_xlabel("Sequence ID", fontsize=12)
    ax.set_ylabel("Edit Distance", fontsize=12)
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    plt.tight_layout()
    # Save dendrogram as PNG
    plt.savefig("dendrogram.png", format="png")
    plt.close()

    # Assign clusters based on threshold
    if threshold is not None:
        cluster_labels = sch.fcluster(linkage_matrix, threshold, criterion="distance")
        cluster_df = pd.DataFrame({"Sequence_ID": unique_ids, "Cluster": cluster_labels})
        # Save clusters to CSV
        cluster_df.to_csv("clusters.csv", index=False)
        return cluster_df
    return None

def main():
    # Set up the command line argument parser
    parser = argparse.ArgumentParser(description="Cluster sequences based on Levenshtein distance.")
    parser.add_argument("--fasta_file", required=True, help="Input FASTA file with sequences")
    parser.add_argument("--substitution_penalty", type=float, default=1, help="Substitution penalty for Levenshtein distance")
    parser.add_argument("--deletion_penalty", type=float, default=1, help="Deletion penalty for Levenshtein distance")
    parser.add_argument("--insertion_penalty", type=float, default=1, help="Insertion penalty for Levenshtein distance")
    parser.add_argument("--out_csv_dendogram", required=True, help="Output CSV filename for cluster data")
    parser.add_argument("--clustering_threshold", type=float, default = 4, help="input threshold for hirearchical clustering of edit distances")
    
    
    # Print help and exit if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()

    # Compute the edit distance table
    edit_distance_df = compute_edit_distance_table(args.fasta_file, args.substitution_penalty, args.deletion_penalty, args.insertion_penalty)
    # Perform clustering and save results
    clusters_df = cluster_sequences(edit_distance_df, method="average", threshold=args.clustering_threshold)
    if clusters_df is not None:
        clusters_df.to_csv(args.out_csv_dendogram, index=False)
        print(f"Clustering results saved to {args.out_csv_dendogram}")
        print("Dendrogram saved as dendrogram.png")

if __name__ == "__main__":
    main()
