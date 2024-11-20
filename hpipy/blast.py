import os
import subprocess
import requests
import tarfile
import pandas as pd
from . import utility as util
from pathlib import Path
from tqdm import tqdm

class BLAST:
    """
    Class for performing BLAST alignments.
    """

    def __init__(self, model = None, hostFile = None, pathogenFile = None, use_slurm = False, num_threads = None, log = None, genome_pool = None, phyloEvalue = None, threshold = None):
        """
        Initialize BLAST object.

        :param model: Name of the host-pathogen model to be used for analysis.
        :param hostFile: The path to the host FASTA file.
        :param pathogenFile: The path to the pathogen FASTA file.
        :param use_slurm: Whether to use SLURM for parallel processing.
        :param num_threads: Number of threads to be used for processing.
        :param log: Logger object for logging messages.
        :param genome_pool:
        :param phyloEvalue:
        :param threshold:
        
        :type model: str
        :type hostFile: str
        :type pathogenFile: str
        :type use_slurm: bool
        :type num_threads: int
        :type log: str
        :type genome_pool:
        :type phyloEvalue:
        :type threshold:

        Attributes:
            databases (dict): A dictionary mapping model names to lists of associated databases.
        """

        self.model = model
        self.hostFile = hostFile
        hostFastaPath = Path(self.hostFile)
        self.hostFasta = hostFastaPath.name
        self.pathogenFile = pathogenFile
        pathogenFastaPath = Path(self.pathogenFile)
        self.pathogenFasta = pathogenFastaPath.name
        self.genome_pool = genome_pool
        self.phyloEvalue = phyloEvalue
        self.threshold = threshold
        self.use_slurm = use_slurm
        self.num_threads = num_threads
        self.log = log
        self.databases = {"humanVirus": ["biogrid", "dip", "hpidb", "intact", "mint", "virhostnet"], 
                          "plantPathogen": ["biogrid", "dip", "hpidb", "intact", "mint"],
                          "animalPathogen": ["biogrid", "dip", "hpidb", "intact", "mint", "virhostnet"],
                          "humanBacteria": ["biogrid", "dip", "hpidb", "intact", "mint"]}
    
    '''
    Interolog BLAST
    '''
    # Index BLAST databases
    def indexBlastDB(self, data_directory, outputdir):
        """
        Index BLAST databases.

        :param data_directory: Directory containing input data for indexing.
        :param outputdir: Results directory to save the output.

        :type data_directory: str
        :type outputdir: str

        :return: Tuple containing the subprocess objects for host and pathogen database creation.
        :rtype: tuple
        """

        dbs = self.databases[self.model]
        dataDir = f"{data_directory}"
        logDir = f"{outputdir}/logs"
        hostDB = None
        pathogenDB = None

        try:
            for db in dbs:
                self.log.info(f"Indexing '{db}' database")
                print(f"HPIpy: Indexing '{db}' database")
                hostDBcmd = f"makeblastdb -in {dataDir}/{self.model}/blast_dbs/{db}_host.fasta -dbtype prot > {logDir}/{self.hostFasta.split('.')[0]}_blastdb_{db}.log"
                hostDB = subprocess.run(hostDBcmd, shell=True, check=True)
                pathogenDBcmd = f"makeblastdb -in {dataDir}/{self.model}/blast_dbs/{db}_pathogen.fasta -dbtype prot > {logDir}/{self.pathogenFasta.split('.')[0]}_blastdb_{db}.log"
                pathogenDB = subprocess.run(pathogenDBcmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            self.log.error(f"An error occurred while generating BLAST databases: {e}")
            print(f"HPIpy: An error occurred while generating BLAST databases: {e}")

        return hostDB, pathogenDB


    # Execute BLAST
    def executeBLAST(self, data_directory, outputdir):
        """
        Execute BLAST alignments for host and pathogen input fasta files.

        :param data_directory: Directory containing indexed database files.
        :param outputdir: Directory to save the BLAST results.
        
        :type data_directory: str
        :type outputdir: str

        :return List of job IDs for submitted SLURM jobs, if any.
        :rtype: list
        """

        dbs = self.databases[self.model]
        blastJobsIDs = []
        dataDir = f"{data_directory}"
        inputClusters = f"{outputdir}/Clustering"
        blast_output = f"{outputdir}/Alignment/Interolog"
        abs_blast_outdir = os.path.abspath(blast_output)

        if not os.path.exists(blast_output):
            os.makedirs(blast_output)
            self.log.info(f"Interolog BLAST results directory created successfully :: {abs_blast_outdir}")
            print(f"HPIpy: Interolog BLAST results directory created successfully :: {abs_blast_outdir}")
        
        for db in dbs:
            # Host
            hostBLASTcmd = f"blastp -db {dataDir}/{self.model}/blast_dbs/{db}_host.fasta -query {inputClusters}/{self.hostFasta} -out {blast_output}/{self.hostFasta.split('.')[0]}_{db}_blast.txt -num_threads {self.num_threads} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'"
            if self.use_slurm:
                job_id, error = util.submit_slurm_job("BLAST_host", hostBLASTcmd, use_slurm=True, outputdir = outputdir)
                if job_id:
                    blastJobsIDs.append(job_id)
                else:
                    self.log.error(f"Host BLAST failed with error: {error}")
                    print(f"HPIpy: Host BLAST failed with error: {error}")
            else:
                self.log.info(f"Running BLASTp alignments for '{self.hostFasta.split('.')[0]}' sequences against '{db}' database")
                print(f"HPIpy: Running BLASTp alignments for '{self.hostFasta.split('.')[0]}' sequences against '{db}' database")
                result = subprocess.run(hostBLASTcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                if result.returncode != 0:
                    self.log.error(f"Host BLAST failed with stderr: {result.stderr.decode('utf-8')}")
                    print(f"HPIpy: Host BLAST failed with stderr: {result.stderr.decode('utf-8')}")

            # Pathogen
            pathogenBLASTcmd = f"blastp -db {dataDir}/{self.model}/blast_dbs/{db}_pathogen.fasta -query {inputClusters}/{self.pathogenFasta} -out {blast_output}/{self.pathogenFasta.split('.')[0]}_{db}_blast.txt -num_threads {self.num_threads} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'"
            if self.use_slurm:
                job_id, error = util.submit_slurm_job("BLAST_pathogen", pathogenBLASTcmd, use_slurm=True, outputdir = outputdir)
                if job_id:
                    blastJobsIDs.append(job_id)
                else:
                    self.log.info(f"Pathogen BLAST failed with error: {error}")
                    print(f"HPIpy: Pathogen BLAST failed with error: {error}")
            else:
                self.log.info(f"Running BLASTp alignments for '{self.pathogenFasta.split('.')[0]}' sequences against '{db}' database")
                print(f"HPIpy: Running BLASTp alignments for '{self.pathogenFasta.split('.')[0]}' sequences against '{db}' database")
                result = subprocess.run(pathogenBLASTcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                if result.returncode != 0:
                    self.log.error(f"Pathogen BLAST failed with stderr: {result.stderr.decode('utf-8')}")
                    print(f"HPIpy: Pathogen BLAST failed with stderr: {result.stderr.decode('utf-8')}")

        return blastJobsIDs


    '''
    Phylogenetic profiling BLAST
    '''
    # Index databases
    def indexPhyloDB(self, data_directory, ngenome, outputdir):
        """
        Index BLAST databases.

        :param data_directory: Directory containing input data in FASTA format.
        :param ngenome: Number of genomes in the selected pool.
        :param outputdir: Results directory to save the output.

        :type data_directory: str
        :type ngenome: int
        :type outputdir: str

        :return: Tuple containing the subprocess objects for host and pathogen database creation and a list of genome names processed.
        :rtype: tuple
        """

        dataDir = f"{data_directory}"
        phyloDB = None
        genomeListPath = f"{dataDir}/{self.genome_pool}/{self.genome_pool}.txt"
        genomeList = pd.read_csv(genomeListPath, names=['genome_name'])
        poolList = genomeList['genome_name'].tolist()

        try:
            with tqdm(total=ngenome, desc="HPIpy: Indexing genome pools", unit="genome") as pbar:
                for i in range(1, ngenome+1):
                    self.log.info(f"Indexing (DIAMOND BLAST) genome {i} ({poolList[i-1]})")                   
                    phyloDBindex = f"diamond makedb --in {dataDir}/{self.genome_pool}/{poolList[i-1]}.fasta --db {dataDir}/{self.genome_pool}/{poolList[i-1]} --quiet"
                    phyloDB = subprocess.run(phyloDBindex, shell=True, check=True)
                    pbar.update(1)
        except subprocess.CalledProcessError as e:
            self.log.error(f"An error occurred while indexing DIAMOND databases: {e}")
            print(f"HPIpy: An error occurred while indexing DIAMOND databases: {e}")

        return phyloDB, poolList
        

    # Execute Diamond BLAST
    def phylo_blast(self, ngenome, phyloPoolList, data_directory, outputdir):
        """
        Perform DIAMOND BLAST for host and pathogen protein sequences.

        :param ngenome: Number of genomes to perform BLAST searches against from the pool.
        :param phyloPoolList: List of genome identifiers in the phylogenetic pool to be used for BLAST searches.
        :param data_directory: Directory containing genome databases.
        :param outputdir: Directory to save the generated BLAST output files.

        :type ngenome: int
        :type phyloPoolList: list
        :type data_directory: str
        :type outputdir: str

        :return: A list of job IDs for submitted SLURM jobs if SLURM is used, or an empty list if SLURM is not used.
        :rtype: list
        """

        dataDir = f"{data_directory}"
        inputClusters = f"{outputdir}/Clustering"
        phylo_blast_output = f"{outputdir}/Alignment/PhyloProfiling"
        abs_phylo_blast_outdir = os.path.abspath(phylo_blast_output)
        phyloBlastJobsIDs = []

        if not os.path.exists(phylo_blast_output):
            os.makedirs(phylo_blast_output)
            self.log.info(f"Phylogenetic profiling BLAST results directory created successfully :: {abs_phylo_blast_outdir}")
            print(f"HPIpy: Phylogenetic profiling BLAST results directory created successfully :: {abs_phylo_blast_outdir}")
        
        host_files = [f"{phyloPoolList[i-1]}_{self.hostFasta.split('.')[0]}_blastOut.txt" for i in range(1, ngenome+1)]
        pathogen_files = [f"{phyloPoolList[i-1]}_{self.pathogenFasta.split('.')[0]}_blastOut.txt" for i in range(1, ngenome+1)]

        with tqdm(total=ngenome, desc="HPIpy: BLAST searching against genomes", unit="genome") as pbar:
            for i in range(1, ngenome+1):
                self.log.info(f"BLAST searching against genome {i} ({phyloPoolList[i-1]}) of {ngenome}")
                
                # Host
                hostPhyloBlast = f"diamond blastp --db {dataDir}/{self.genome_pool}/{phyloPoolList[i-1]} --query {inputClusters}/{self.hostFasta} --evalue {self.phyloEvalue} --out {phylo_blast_output}/{host_files[i-1]} --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp --max-target-seqs 1 --threads {self.num_threads}"
                if self.use_slurm:
                    job_id, error = util.submit_slurm_job("Phylo_BLAST_host", hostPhyloBlast, use_slurm=True, outputdir = outputdir)
                    if job_id:
                        phyloBlastJobsIDs.append(job_id)
                    else:
                        self.log.error(f"Host DIAMOND BLAST failed with error: {error}")
                        print(f"HPIpy: Host DIAMOND BLAST failed with error: {error}")
                else:
                    result = subprocess.run(hostPhyloBlast, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    if result.returncode != 0:
                        self.log.error(f"Host DIAMOND BLAST failed with stderr: {result.stderr.decode('utf-8')}")
                        print(f"HPIpy: Host DIAMOND BLAST failed with stderr: {result.stderr.decode('utf-8')}")

                # Pathogen
                pathogenPhyloBlast = f"diamond blastp --db {dataDir}/{self.genome_pool}/{phyloPoolList[i-1]} --query {inputClusters}/{self.pathogenFasta} --evalue {self.phyloEvalue} --out {phylo_blast_output}/{pathogen_files[i-1]} --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp --max-target-seqs 1 --threads {self.num_threads}"
                if self.use_slurm:
                    job_id, error = util.submit_slurm_job("Phylo_BLAST_pathogen", pathogenPhyloBlast, use_slurm=True, outputdir = outputdir)
                    if job_id:
                        phyloBlastJobsIDs.append(job_id)
                    else:
                        self.log.error(f"Pathogen DIAMOND BLAST failed with error: {error}")
                        print(f"HPIpy: Pathogen DIAMOND BLAST failed with error: {error}")
                else:
                    result = subprocess.run(pathogenPhyloBlast, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    if result.returncode != 0:
                        self.log.error(f"Pathogen DIAMOND BLAST failed with stderr: {result.stderr.decode('utf-8')}")
                        print(f"HPIpy: Pathogen DIAMOND BLAST failed with stderr: {result.stderr.decode('utf-8')}")
                pbar.update(1)

        return phyloBlastJobsIDs