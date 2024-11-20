import pandas as pd
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
from Levenshtein import distance
import random

class PhyloProfiling:
    """
    Class to perform phylogenetic profiling-based interaction prediction analysis between host and pathogen proteins.
    """

    def __init__(self):
        """
        Initialize PhyloProfiling class.
        """
        
        return


    # Extract protein IDs
    def extractIDs(self, hostProtFasta, pathogenProtFasta):
        """
        Extract input sequence IDs and initialize patterns for host and pathogen proteins.
        
        :param hostProtFasta: File path to the FASTA file containing host protein sequences.
        :param pathogenProtFasta: File path to the FASTA file containing pathogen protein sequences.

        :type hostProtFasta: str
        :type pathogenProtFasta: str

        :return: A tuple containing:
             - List of host protein sequence IDs.
             - List of pathogen protein sequence IDs.
             - Number of host protein sequences.
             - Number of pathogen protein sequences.
             - Dictionary to store patterns for host protein sequences.
             - Dictionary to store patterns for pathogen protein sequences.
        :rtype: tuple
        """

        # Host sequences
        host_seq = list(SeqIO.parse(f"{hostProtFasta}", 'fasta'))
        hostIDs = [host.id for host in host_seq]
        numberHost = len(hostIDs)
        pattern_host = {i: "" for i in range(numberHost)}

        # Pathogen sequences
        pathogen_seq = list(SeqIO.parse(f"{pathogenProtFasta}", 'fasta'))
        pathogenIDs = [pathogen.id for pathogen in pathogen_seq]
        numberPathogen = len(pathogenIDs)
        pattern_pathogen = {k: "" for k in range(numberPathogen)}

        return hostIDs, pathogenIDs, numberHost, numberPathogen, pattern_host, pattern_pathogen
    

    # Extract pattern
    def process_files(self, files, ids, phyloIdent, phyloCov):
        """
        Process BLAST output files to generate patterns based on sequence identity and coverage thresholds.

        :param files: List of file paths to BLAST output files.
        :param ids: List of sequence IDs to be checked against the BLAST results.
        :param phyloIdent: Identity threshold for filtering BLAST results.
        :param phyloCov: Coverage threshold for filtering BLAST results.

        :type files: list
        :type ids: list
        :type phyloIdent: float
        :type phyloCov: float

        :return: A list of patterns representing whether each sequence ID meets the specified thresholds across all files.
        :rtype: list
        """

        patterns = [''] * len(ids)
        for file in files:
            with open(file) as f:
                blast_ids = set()
                for line in f:
                    parts = line.split("\t")
                    if len(parts) > 5:  # Ensure there are enough columns
                        identity, coverage = parts[2], parts[5]
                        try:
                            identity = float(identity)
                            coverage = float(coverage)
                            if identity > phyloIdent and coverage > phyloCov:
                                blast_ids.add(parts[0])
                        except ValueError:
                            # Handle cases where conversion to float fails
                            continue
                
            for i, id in enumerate(ids):
                patterns[i] += '1' if id in blast_ids else '0'
    
        return patterns


    # Calculate phylogenetic distance
    def compute_similarity(self, args):
        """
        Compute similarity scores between host and pathogen protein patterns based on a threshold.
        
        :param args: A tuple containing the following arguments:
            - i (int): Index of the host pattern in the overall list.
            - host_pattern (str): The binary pattern representing the presence/absence of features in the host protein.
            - host (str): ID of the host protein.
            - nullPool (str): A pattern used to represent the absence of any features.
            - pattern_pathogen (list): A list of binary patterns representing the presence/absence of features in the pathogen proteins.
            - pathogenIDs (list): A list of pathogen protein IDs corresponding to the patterns in `pattern_pathogen`.
            - ngenome (int): Total number of genomes in the pool.
            - threshold (float): The minimum similarity score required to include a host-pathogen pair in the results.

        :type args: tuple

        :return: A list of lists that contains host protein ID, pathogen protein ID and similarity score.
        :rtype: list
        """

        i, host_pattern, host, nullPool, pattern_pathogen, pathogenIDs, ngenome, threshold = args
        result = []
        
        if host_pattern != nullPool:
            for k, pathogen_pattern in enumerate(pattern_pathogen):
                if pathogen_pattern != nullPool:
                    score = round((ngenome - distance(host_pattern, pathogen_pattern)) / ngenome, 2)
                    if score > threshold:
                        result.append([host, pathogenIDs[k], score])

        return result


    # Predicting PPIs
    def phylo_ppis(self, numberHost, pattern_host, pattern_pathogen, nullPool, hostIDs, pathogenIDs, ngenome, threshold, chunk_size=10000):
        """
        Perform phylogenetic profiling to predict protein-protein interactions between host and pathogen proteins.

        :param numberHost: Number of host proteins to process.
        :param pattern_host: A list of binary patterns representing the presence/absence of features in host proteins.
        :param pattern_pathogen: A list of binary patterns representing the presence/absence of features in pathogen proteins.
        :param nullPool: A pattern used to represent the absence of any features.
        :param hostIDs: A list of host protein IDs corresponding to the patterns in `pattern_host`.
        :param pathogenIDs: A list of pathogen protein IDs corresponding to the patterns in `pattern_pathogen`.
        :param ngenome: Total number of genomes in the pool used for pattern generation.
        :param threshold: The minimum similarity score required to include a host-pathogen pair in the results.
        :param chunk_size: The number of host proteins to process in each chunk for parallel computation. Default is 10,000.

        :type numberHost: int
        :type pattern_host: list
        :type pattern_pathogen: list
        :type nullPool: str
        :type hostIDs: list
        :type pathogenIDs: list
        :type ngenome: int
        :type threshold: float
        :type chunk_size: int

        :return: A DataFrame containing unique host-pathogen pairs with their computed similarity scores.
        :rtype: pandas.DataFrame
        """

        random.seed(42)
        
        pool = Pool(cpu_count())
        args = [(i, pattern_host[i], hostIDs[i], nullPool, pattern_pathogen, pathogenIDs, ngenome, threshold) for i in range(numberHost)]
        
        all_results = []
        chunks = [args[i:i + chunk_size] for i in range(0, len(args), chunk_size)]
        
        for count, chunk in enumerate(chunks, 1):
            print(f" >> Processing chunk: {count}")

            results = pool.map(self.compute_similarity, chunk)
            results = [item for sublist in results for item in sublist]
            results_df = pd.DataFrame(results, columns=['Host', 'Pathogen', 'Score'])
            results_df = results_df.drop_duplicates(subset=['Host', 'Pathogen'])
            all_results.append(results_df)
        
        pool.close()
        pool.join()
        
        combined_results = pd.concat(all_results, ignore_index=True)
        
        return combined_results