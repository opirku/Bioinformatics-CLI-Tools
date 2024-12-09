#!/opt/anaconda3/envs/Phylo/bin/python
from ete3 import Tree
import sys
from typing import Dict

class getSpeciesNames:
    def __init__(self, fastaFile: str):
        self.fasta_file = fastaFile
        self.accessionSpeciesDictionary: Dict[str, str] = {}

    def readFastaFile(self):
        """Reads the FASTA file line by line."""
        with open(self.fasta_file, "r") as openedFastaFile:
            for fastaFileLine in openedFastaFile:
                if fastaFileLine.startswith(">"):
                    yield fastaFileLine.strip()[1:]

    def makeAccessionSpeciesDictionary(self):
        """Populates the accession-to-species name dictionary."""
        for fastaHeader in self.readFastaFile():
            headerParts = fastaHeader.split()
            if len(headerParts) >= 3:  # Ensure there are enough parts
                accessions = headerParts[0]
                speciesNames = headerParts[2]
                self.accessionSpeciesDictionary[accessions] = speciesNames

    def getDictionary(self) -> Dict[str, str]:
        """Returns the populated dictionary."""
        return self.accessionSpeciesDictionary


class ReplaceNewickLeafNames:
    def __init__(self, newick_tree_file: str, accessionSpeciesDictionary: Dict[str, str]):
        self.newick_tree_file = newick_tree_file
        self.accessionSpeciesDictionary = accessionSpeciesDictionary

    def openNewickTree(self):
        """Parses the Newick tree file."""
        with open(self.newick_tree_file, "r") as newickTreeFile:
            newickTreeString = newickTreeFile.read().strip()  # Read and strip the content

        try:
            newickTree = Tree(newickTreeString, format=1)
            return newickTree
        except Exception as e:
            print(f"Error parsing Newick tree: {e}")
            return None

    def replaceLeafNames(self):
        """Replaces leaf names in the tree with species names."""
        tree = self.openNewickTree()

        if tree:
            for leaf in tree.get_leaves():
                accession = leaf.name
                speciesNames = self.accessionSpeciesDictionary.get(accession)

                if speciesNames:
                    leaf.name = speciesNames

            return tree
        else:
            print("Tree parsing failed.")
            return None


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: changeLeafNames.py <newick_file> <fasta_file> <output_file>")
        sys.exit(1)

    newick_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    # Step 1: Create an instance of getSpeciesNames and populate the dictionary
    species_names = getSpeciesNames(fasta_file)
    species_names.makeAccessionSpeciesDictionary()

    # Step 2: Create an instance of ReplaceNewickLeafNames
    replace_names = ReplaceNewickLeafNames(newick_file, species_names.getDictionary())

    # Step 3: Replace leaf names in the tree
    updated_tree = replace_names.replaceLeafNames()

    if updated_tree:
        with open(output_file, "w") as out_file:
            out_file.write(updated_tree.write(format=1))
            print(f"Updated Newick tree written to {output_file}")

