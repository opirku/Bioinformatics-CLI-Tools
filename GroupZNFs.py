"""
Author: Orion Pirku
Date: 26.03.2025
Affiliation: Master Student/Research Assistant, 
Max Planck Institute for Evolutionary Biology, 
Ploen, Germany

 
"""
import argparse
import dis
from turtle import distance
import cython
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse
from numba import njit

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

def group_exact_sequences(filtered_sequences: SeqRecord) -> dict[str, List[str]]: 
    
    """The function does the following:
    1. Make a set of unique sequences
    2. initialize a empty dictionary with unique sequences as keys
    3. add sequence id that have that sequence as a list value to that dictionary
    """ 
    
    unique_sequences: set[str] = {str(record.seq) for record in filtered_sequences}
    
    sequence_groups: dict[str, list[str]] ={sequence : [] for sequence in unique_sequences}
    
    for record in filtered_sequences:
            sequence = str(record.seq)
            
            sequence_groups[sequence].append(record.id)
    
    return sequence_groups

def aminoacid_edit_distance(sequence_1: str, sequence_2: str, deletion_cost: float, insertion_cost: float, substitution_cost: float) -> int:
        
        n: int = len(sequence_1)
        m: int = len(sequence_2)
        
        distance_matrix = [[0] * (m + 1) for _ in range(n + 1)]
        
        for i in range(n+1):
                distance_matrix[i][0] = i
        for j in range(m+1):
                distance_matrix[0][j] = j 
        
        for i range(1, n+1):
                
                for j in range(1, m+1):
                        
                        if sequence_1[i-1] == sequence_2[j-1]:
                                cost: int = 0
                        else:
                                cost: float = substitution_cost
                        
                        distance_matrix[i][j] = min(
                                distance_matrix[i-1][j] + deletion_cost, 
                                distance_matrix[i, j-1] + insertion_cost, 
                                distance_matrix[i-1][j-1] + cost
			)
        
        return distance_matrix[n+1][m+1]