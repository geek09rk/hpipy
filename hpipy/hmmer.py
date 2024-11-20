import os
import requests
import subprocess
from pathlib import Path
from . import utility as util

class HMMER:
    """
    A class for performing HMMER (profile hidden Markov model) analysis.
    """

    def __init__(self, model = None, hostFile = None, pathogenFile = None,  use_slurm = False, num_threads = None, log = None):
        """
        Initialize HMMER object.

        :param model: Name of the host-pathogen model to be used for analysis.
        :param hostFile: The path to the host FASTA file.
        :param pathogenFile: The path to the pathogen FASTA file.
        :param use_slurm: Whether to use SLURM for parallel processing.
        :param num_threads: Number of threads to be used for processing.
        :param log: Logger object for logging messages.

        :type model: str
        :type hostFile: str
        :type pathogenFile: str
        :type use_slurm: bool
        :type num_threads: int
        :type log: str
        """
        
        self.model = model
        self.hostFile = hostFile
        hostFastaPath = Path(self.hostFile)
        self.hostFasta = hostFastaPath.name
        self.pathogenFile = pathogenFile
        pathogenFastaPath = Path(self.pathogenFile)
        self.pathogenFasta = pathogenFastaPath.name
        self.use_slurm = use_slurm
        self.num_threads = str(num_threads)
        self.log = log
        self.ppiModel = {"humanVirus": [0.2, 10], 
                         'plantPathogen': [0.2, 0.45],
                         "animalPathogen": [0.2, 0.45],
                         "humanBacteria": [0.2, 0.35]}


    # hmmpress
    def pfam(self, data_directory, outputdir):
        """
        Executes Pfam analysis.

        :param data_directory: Directory to save the downloaded Pfam database.
        :param outputdir: Results directory to save the output.

        :type data_directory: str
        :type outputdir: str

        :return The directory containing Pfam database files.
        :rtype: str
        """

        pfamOutdir = f"{data_directory}/PfamDB"
        pfam_files = [os.path.join(pfamOutdir, f"Pfam-A{ext}") for ext in ['.hmm', '.hmm.h3m', '.hmm.h3i', '.hmm.h3f', '.hmm.h3p']]
        pfam_files_exist = os.path.exists(pfamOutdir) and all(os.path.exists(file) for file in pfam_files)

        if pfam_files_exist:
            self.log.info("Required Pfam files exist, proceeding to next step")
            print("HPIpy: Required Pfam files exist, proceeding to next step...")
        else:
            if not os.path.exists(pfamOutdir):
                os.mkdir(pfamOutdir)
                self.log.info("Pfam directory created successfully")
                print("HPIpy: Pfam directory created successfully")

            # Pfam from InterPro
            self.log.info("Downloading the latest Pfam database from InterPro (https://www.ebi.ac.uk/interpro/)")
            print("HPIpy: Downloading latest Pfam database from InterPro (https://www.ebi.ac.uk/interpro/)")

            pfam_url = 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
            pfam_path = os.path.join(pfamOutdir, 'Pfam-A.hmm.gz')

            try:
                pfam_db = requests.get(pfam_url, allow_redirects=True)
                pfam_db.raise_for_status()
                with open(pfam_path, 'wb') as f:
                    f.write(pfam_db.content)
                self.log.info("Pfam database download complete")
                print("HPIpy: Pfam database downloaded successfully")
            except requests.RequestException as e:
                self.log.error(f"Failed to download Pfam database: {e}")
                print(f"HPIpy: Failed to download Pfam database: {e}")
                return None

            # Unzip Pfam database
            self.log.info("Processing Pfam database files")
            print("HPIpy: Processing Pfam database files")

            try:
                subprocess.run(f"gunzip -c {pfam_path} > {os.path.join(pfamOutdir, 'Pfam-A.hmm')}", shell=True, check=True)
            except subprocess.CalledProcessError as e:
                self.log.error(f"Failed to unzip Pfam database: {e}")
                print(f"HPIpy: Failed to unzip Pfam database: {e}")
                return None

            # Run hmmpress
            self.log.info("Executing 'hmmpress'")
            print("HPIpy: Running 'hmmpress'")

            hmmpress_cmd = f"hmmpress {os.path.join(pfamOutdir, 'Pfam-A.hmm')} > {os.path.join(outputdir, 'logs', 'hmmpress.log')}"
            try:
                subprocess.run(hmmpress_cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                self.log.error(f"Failed to run hmmpress: {e}")
                print(f"HPIpy: Failed to run hmmpress: {e}")
                return None

        return pfamOutdir
    
    
    # hmmscan
    def hmmer(self, PfamDir = None, outputdir = None):
        """
        Executes HMMER analysis.

        :param PfamDir: Directory containing Pfam files.
        :param outputdir: Results directory to save the output.
        
        :type PfamDir: str
        :type outputdir: str

        :return: A list of SLURM job IDs for submitted jobs, if any.
        :rtype: list
        """

        logDir = f"{outputdir}/logs"
        inputClusters = f"{outputdir}/Clustering"
        hmmer_output = f"{outputdir}/Domains"
        abs_hmm_outdir = os.path.abspath(hmmer_output)

        if not os.path.exists(hmmer_output):
            os.mkdir(hmmer_output)
            self.log.info(f"HMMER results directory created successfully :: {abs_hmm_outdir}")
            print(f"HPIpy: HMMER results directory created successfully :: {abs_hmm_outdir}")

        self.log.info("HMMER analysis started. Time duration depends on the number of sequences")
        print("HPIpy: Running 'hmmscan'. Time duration depends on the number of sequences")

        hmmerJobsIDs = []

        def execute_hmmscan(command, job_name):
            try:
                subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
                self.log.info(f"{job_name} completed successfully.")
                print(f"HPIpy: {job_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                self.log.error(f"{job_name} failed with error: {e}")
                print(f"HPIpy: {job_name} failed with error: {e}")
                return None, str(e)
            return job_name, None

        # Host
        hostHMMSCANcmd = f"hmmscan --domE {self.ppiModel[self.model][0]} --cpu {self.num_threads} --tblout {hmmer_output}/{self.hostFasta.split('.')[0]}_domains.txt {PfamDir}/Pfam-A.hmm {inputClusters}/{self.hostFasta} > {logDir}/{self.hostFasta.split('.')[0]}_hmmscan.log"
        if self.use_slurm:
            job_id, error = util.submit_slurm_job("HMM_host", hostHMMSCANcmd, use_slurm=True, outputdir = outputdir)
            if job_id:
                hmmerJobsIDs.append(job_id)
            else:
                self.log.info(f"Host HMMER failed with error: {error}")
                print(f"HPIpy: Host HMMER failed with error: {error}")
        else:
            job_id, error = execute_hmmscan(hostHMMSCANcmd, "Host HMMER")

        # Pathogen
        pathogenHMMSCANcmd = f"hmmscan --domE {self.ppiModel[self.model][1]} --cpu {self.num_threads} --tblout {hmmer_output}/{self.pathogenFasta.split('.')[0]}_domains.txt {PfamDir}/Pfam-A.hmm {inputClusters}/{self.pathogenFasta} > {logDir}/{self.pathogenFasta.split('.')[0]}_hmmscan.log"
        if self.use_slurm:
            job_id, error = util.submit_slurm_job("HMM_pathogen", pathogenHMMSCANcmd, use_slurm=True, outputdir = outputdir)
            if job_id:
                hmmerJobsIDs.append(job_id)
            else:
                self.log.info(f"Pathogen HMMER failed with error: {error}")
                print(f"HPIpy: Pathogen HMMER failed with error: {error}")
        else:
            job_id, error = execute_hmmscan(pathogenHMMSCANcmd, "Pathogen HMMER")

        return hmmerJobsIDs