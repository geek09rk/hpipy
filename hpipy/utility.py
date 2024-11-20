import os
import subprocess
import requests
import tarfile
from Bio import SeqIO
from .args_parse import *


## Create new directory
def make_directory(directory):
    """
    Create a new directory or modify the name if it already exists.

    :param directory: The name of the directory to be created or modified name if it already exists.
    :type directory: str

    :return: The name of the created or modified directory.
    :rtype: str
    """
    
    if not os.path.exists(directory):
        os.mkdir(directory)
    else:
        i = 1
        while os.path.exists(f'{directory}_{i}'):
            i += 1
        os.mkdir(f'{directory}_{i}')
        return f'{directory}_{i}'
    return directory


## Compressed file check
def decompress_file(file_path):
    """
    Decompresses a compressed file, if necessary.

    :param file_path: The path to the file to be decompressed.
    :type file_path: str

    :return: Path to the decompressed file.
    :rtype: str

    :raises ValueError: If no files are found in the extracted folder (for .zip files)
    """

    _, file_extension = os.path.splitext(file_path)
    
    if file_extension == '.gz':
        decompressed_path = file_path[:-3]
        os.system(f"gzip -d -c {file_path} > {decompressed_path}")
    
    elif file_extension == '.zip':
        decompressed_folder = file_path[:-4]
        os.system(f"unzip {file_path} -d {decompressed_folder}")
        
        extracted_files = os.listdir(decompressed_folder)
        if extracted_files:
            decompressed_path = os.path.join(decompressed_folder, extracted_files[0])
        else:
            raise ValueError("No files found in the extracted folder")
        
    else:
        decompressed_path = file_path
    
    return decompressed_path


## Fasta format check
def is_fasta(filename):
    """
    Checks if the input file is in FASTA format

    :param filename: The name of the file to be checked.
    :type filename: str

    :return: True if the file is in FASTA format, FALSE otherwise.
    :rtype: bool
    """
    
    file_extension = os.path.splitext(filename)[1]
    if file_extension.lower() in ['.fasta', '.fa', '.faa']:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)
    else:
        return False


## Validate input protein sequences
def is_protein(fasta_file):
    """
    Checks if the input FASTA file contains protein sequences.

    :param fasta_file: The path to the FASTA file to be validated.
    :type fasta_file: str

    :return: True if the file contains only valid protein sequences, False otherwise.
    :rtype: bool
    """
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            for char in line.strip():
                if char not in "ACDEFGHIKLMNOPQRSTUVWXY":
                    return False
    return True


## Sequence homology
def cdhit(inputFasta, outputdir, outfile):
    """
    Performs sequence clustering using CD-HIT.

    :param inputFasta: The path to the input FASTA file containing protein sequences.
    :param outputdir: The directory where the output files will be stored.
    :param outfile: The name of the output file containing clustered sequences.

    :type inputFasta: str
    :type outputdir: str
    :type outfile: str

    :return: The number of clusters formed/identified.
    :rtype: int
    """

    cdhitDir = f"{outputdir}/Clustering"
    if not os.path.exists(cdhitDir):
        os.mkdir(cdhitDir)

    cmd = f"cd-hit -i {inputFasta} -o {cdhitDir}/{outfile} -c {args.seq_homology} > {outputdir}/logs/{outfile.split('.fasta')[0]}_cdhit_out.log"
    os.system(cmd)

    clusters_cmd = f'grep -c ">" {cdhitDir}/{outfile}'
    clusters_output = subprocess.run(clusters_cmd, shell=True, capture_output=True, text=True)

    if clusters_output.returncode == 0:
        num_clusters = int(clusters_output.stdout.strip())
        print(f"HPIpy: Number of clusters identified in '{outfile}': {num_clusters}")
    else:
        print("Error occurred while counting clusters.")

    return num_clusters


## Obtain data
def downloadData(data_directory):
        """
        Downlaod required data files.

        :param data_directory: Directory to save the downloaded files.
        :type data_directory: str
        """

        print(f"HPIpy: Downloading required data files...")

        dataDir = f"{data_directory}"

        if isinstance(args.computation, list):
            args.computation = args.computation[0] if len(args.computation) == 1 and isinstance(args.computation[0], list) else args.computation

        urls_to_download = []

        for method in args.computation:
            if method == 'interolog':
                url = f'https://kaabil.net/hpipy/data/{args.model}.tar.gz'
                urls_to_download.append(url)

            elif method == 'domain':
                url = f'https://kaabil.net/hpipy/data/domainDBs.tar.gz'
                urls_to_download.append(url)

            elif method == 'phyloProfiling':
                url = f'https://kaabil.net/hpipy/data/{args.genome_pool}.tar.gz'
                urls_to_download.append(url)

        # Download data file
        for url in urls_to_download:
            filename = url.split('/')[-1]

            try:
                response = requests.get(url)
                response.raise_for_status()  # bad responses
            except requests.exceptions.RequestException as e:
                print(f"HPIpy: Failed to download the file from {url}. Error: {e}")
                print("HPIpy stopped.")
                exit()
        
            file_path = os.path.join(dataDir, filename)
            with open(file_path, 'wb') as f:
                f.write(response.content)
            
            # Decompress file
            try:
                with tarfile.open(file_path, 'r:gz') as tar:
                    tar.extractall(path=dataDir)
            except tarfile.TarError as e:
                print(f"HPIpy: Failed to decompress the file: {e}")
                print("HPIpy stopped.")
                exit()
            finally:
                if os.path.exists(file_path):
                    os.remove(file_path)

        return


## Obtain annotations
def downloadAnnot(data_directory):
        """
        Downlaod required annotation files.

        :param data_directory: Directory to save the downloaded files.
        :type data_directory: str
        """

        print(f"HPIpy: Downloading annotations...")

        dataDir = f"{data_directory}"
        url = None

        if args.model == "humanVirus" or args.model == "humanBacteria":
            url = f'https://kaabil.net/hpipy/data/annotations.tar.gz'
        else:
            pass
        
        # URL
        filename = url.split('/')[-1]

        try:
            response = requests.get(url)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f"HPIpy: Failed to download annotations from {url}. Error: {e}")
            print("HPIpy stopped.")
            exit()

        file_path = os.path.join(dataDir, filename)
        with open(file_path, 'wb') as f:
            f.write(response.content)
        
        # Decompress file
        try:
            with tarfile.open(file_path, 'r:gz') as tar:
                tar.extractall(path=dataDir)
        except tarfile.TarError as e:
            print(f"HPIpy: Failed to decompress the file: {e}")
            print("HPIpy stopped.")
            exit()
        finally:
            if os.path.exists(file_path):
                os.remove(file_path)

        print("HPIpy: Required annotations downloaded successfully")

        return


## Files generation by makeblastdb
def check_blastDB_files(data_directory):
    """
    Checks if all the required files generated by "makeblastdb" are present.

    :param data_directory: The directory containing the indexed BLAST databases.
    :type data_directory: str

    :return: True if all required files are present, False otherwise.
    :rtype: bool
    """
    
    dataDir = f"{data_directory}"

    databases = {"humanVirus": ["biogrid", "dip", "hpidb", "intact", "mint", "virhostnet"], 
                 "plantPathogen": ["biogrid", "dip", "hpidb", "intact", "mint"], 
                 "animalPathogen": ["biogrid", "dip", "hpidb", "intact", "mint", "virhostnet"], 
                 "humanBacteria": ["biogrid", "dip", "hpidb", "intact", "mint"]}
    
    dbs = databases[args.model]

    for db in dbs:
        required_extensions = [".pdb", ".phr", ".pin", ".pjs", ".pot", ".psq", ".ptf", ".pto"]
        for extension in required_extensions:
            if not (os.path.exists(f"{dataDir}/{args.model}/blast_dbs/{db}_host.fasta{extension}") and
                    os.path.exists(f"{dataDir}/{args.model}/blast_dbs/{db}_pathogen.fasta{extension}")):
                return False
    return True


## Files generation by hmmpress
def check_hmm_files(directory):
    """
    Checks if the required files generated by "hmmpress" are present.

    :param directory: The directory containing the HMM files.
    :type directory: str

    :return: True if all the required files are present, False otherwise.
    :rtype: bool
    """
    
    required_extensions = [".h3m", ".h3i", ".h3f", ".h3p"]
    for extension in required_extensions:
        if not os.path.exists(f"{directory}/Pfam-A.hmm{extension}"):
            return False
    return True


## Submit job using SLURM
def submit_slurm_job(job_name, command, use_slurm=False, dep='', outputdir = None, account=None):
    """
    Submits a job using SLURM manager.

    :param job_name: The name of the job.
    :param command: The command to be executed by the job.
    :param use_slurm: Flag indicating whether to use SLURM for job submission.
    :param dep: Job dependencies.
    :param outputdir: The directory where job output will be stored.
    :param account: SLURM account name.
    
    :type job_name: str
    :type command: str
    :type use_slurm: bool
    :type dep: str
    :type outputdir: str
    :type account: str

    :return: A tuple containing the job ID if the submission was successful and an error message if any occurred.
    :rtype: tuple
    """
    
    job_id = None
    error = None
    try:
        if dep != '':
            dep = '--dependency=afterok:{} --kill-on-invalid-dep=yes '.format(dep)

        command = command.replace("'", "\"")
        sbatch_command = f"sbatch -J {job_name} -o {outputdir}/logs/{job_name}.out -e {outputdir}/logs/{job_name}.err -t 4-00:00 --cpus-per-task={args.num_threads} --ntasks=1 --wrap='{command}' {dep}"

        if account:
            sbatch_command += f" --account={account}"
        
        sbatch_command += f" --wrap='{command}' {dep}"
        sbatch_response = subprocess.getoutput(sbatch_command)
        job_id = sbatch_response.split(' ')[-1].strip()

    except Exception as e:  
       error = str(e)

    return job_id, error


## Check running SLURM jobs
def is_job_running(job_id):
    """
    Checks if a SLURM job is currently running using its ID.

    :param job_id: The ID of the job to check.
    :type job_id: str

    :return: True if the job is running, False otherwise.
    :rtype: bool
    """
    
    result = subprocess.run(["squeue", "-j", job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if str(result.stderr.decode("utf-8")).startswith("slurm_load_jobs error: Invalid job id specified"):
        return False
    else:
        return True
    

## Check DB files
def check_db_files(dbDirectory):
    """
    Checks if all required database files are present.

    :param dbDirectory: The directory containing the database files.
    :type dbDirectory: str

    :return: True if all required files are present, False otherwise.
    :rtype: bool
    """
    
    if not os.path.exists(f"{dbDirectory}/interactions.db") and \
        not os.path.exists(f"{dbDirectory}/{args.host.split('.')[0]}.db") and \
        not os.path.exists(f"{dbDirectory}/{args.pathogen.split('.')[0]}.db"):
        return False
    return True


## Extract protein IDs & sequences
def extract_sequences(input_fasta, protein_ids, proteinCol, output_dir, computation):   # organism
    """
    Extracts sequences from an input FASTA file based on a list of protein IDs.

    :param input_fasta: Path to the input FASTA file.
    :param protein_ids: Dataframe of protein IDs to extract.
    :param proteinCol: Column name in the DataFrame that contains the protein IDs.
    :param output_dir: Directory to save the output FASTA file.
    :param organism: Organism name to include in the output file.

    :type input_fasta: str
    :type protein_ids: pandas.DataFrame
    :type proteinCol: str
    :type output_dir: str
    :type organism: str

    :return: None   
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    sequences_dir = os.path.join(output_dir, "extracted_sequences")
    os.makedirs(sequences_dir, exist_ok=True)
    output_fasta = os.path.join(sequences_dir, f"{computation}_{proteinCol}_sequences.fasta")
    protein_ids_set = set(protein_ids[proteinCol])

    with open(f"{sequences_dir}/{computation}_{proteinCol}_IDs.txt", "w") as f:
        for protein_id in protein_ids_set:
            f.write(f"{protein_id}\n")
    
    sequences_to_write = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        for pid in protein_ids_set:
            if pid in record.id:
                sequences_to_write.append(record)
                break

    SeqIO.write(sequences_to_write, output_fasta, "fasta")

    return


## Remove file
def delFile(filePath):
    """
    Removes a file.

    :param filePath: The path to the file to be removed.
    :type filePath: str

    :return: Return code of the removal operation (0 if successful, 1 if an error occurred).
    :rtype: int
    """

    try:
        os.remove(filePath)
        return 0
    except Exception as e:
        print(f"Error removing file {filePath}: {e}")
        return 1


## Remove directory
def del_directory(dirPath):
    """
    Recursively removes a directory and all its contents.

    :param dirPath: The path to the directory to be removed.
    :type dirPath: str

    :return: Return code of the removal operation (0 if successful, 1 if an error occurred).
    :rtype: int
    """
    
    dir_cmd = os.system(f"rm -rf {dirPath}")

    return dir_cmd
