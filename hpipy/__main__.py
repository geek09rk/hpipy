import os
import sys
import time
import re
from .args_parse import *
from . import utility as util
from . import logger
from . import blast
from . import hmmer
from . import goSimilarity
from . import tables
from . import predict
from pathlib import Path


def main():

    if not args.resume_ppis:
        print("\n -------------------------------------------------------------------------")
        print("  HPIpy: A tool for host-pathogen protein-protein interactions prediction ")
        print(" -------------------------------------------------------------------------\n")

        # Results directory
        outputdir = util.make_directory(args.outputdir)
        abs_outputdir = os.path.abspath(outputdir)
        pkgData = os.path.join(os.path.dirname(outputdir), "HPIpy_data")
        if not os.path.exists(pkgData):
            os.makedirs(pkgData)

        # Log directory
        logsDir = util.make_directory(os.path.join(outputdir, "logs"))
        abs_logdir = os.path.abspath(logsDir)
        log = logger.logHPIpy(logdir = outputdir, mode = "a")

        # Executed command
        command = ""
        for i, arg in enumerate(sys.argv):
            command += arg + " "
        log.info(f"Command executed: {command.strip()}")

        log.info(f"Successfully created results directory :: {abs_outputdir}")
        print(f"HPIpy: Successfully created output directory :: {abs_outputdir}")
        log.info(f"Successfully created logs directory :: {abs_logdir}")
        print(f"HPIpy: Successfully created logs directory :: {abs_logdir}")
        
        # Processing input files
        hostInputFile = args.host
        hostDecomp = util.decompress_file(hostInputFile)
        hostFastaPath = Path(hostDecomp)
        hostFastaName = hostFastaPath.name
        pathogenInputFile = args.pathogen
        pathogenDecomp = util.decompress_file(pathogenInputFile)
        pathogenFastaPath = Path(pathogenDecomp)
        pathogenFastaName = pathogenFastaPath.name

        # GO terms files
        hostGOPath = None
        pathogenGOPath = None
        if args.hostGOFile:
            hostGOFile = args.hostGOFile
            hostGOPath = Path(hostGOFile)
        if args.pathogenGOFile:
            pathogenGOFile = args.pathogenGOFile
            pathogenGOPath = Path(pathogenGOFile)

        print("HPIpy: Parameters provided by the user for analysis:")
        log.info("Parameters provided by the user for analysis:")
        
        # Parameters
        related_args = {
            'interIdentity': ['interCoverage', 'interEvalue', 'outputdir', 'network', 'num_threads'],  # if 'input_file' is provided, also show these args
            'domHostEvalue': ['domPathogenEvalue', 'outputdir', 'network', 'num_threads'],
            'phyloThreshold': ['genome_pool', 'phyloIdentity', 'phyloCoverage', 'phyloEvalue', 'outputdir', 'network', 'num_threads'],
            'pathogenGOFile': ['goSimThreshold', 'interproscan', 'outputdir', 'network', 'num_threads'],
            'interproscan': ['hostGOFile', 'pathogenGOFile', 'goSimThreshold', 'outputdir', 'network', 'num_threads']
        }

        printed_args = set()
        for arg in vars(args):
            # Check if the argument is set to a non-default value
            if getattr(args, arg) != parser.get_default(arg):
                print(f"    > {arg.replace('_', ' ').title()}: '{getattr(args, arg)}'")
                log.info(f"    > {arg.replace('_', ' ').title()}: '{getattr(args, arg)}'")
                printed_args.add(arg)
                if arg in related_args:
                    for related_arg in related_args[arg]:
                        if related_arg not in printed_args:
                            print(f"    > {related_arg.replace('_', ' ').title()}: '{getattr(args, related_arg)}'")
                            log.info(f"    > {related_arg.replace('_', ' ').title()}: '{getattr(args, related_arg)}'")
                            printed_args.add(related_arg)


        ## Validate fasta format
        log.info("Input FASTA files validation started")

        # Host
        if util.is_fasta(hostDecomp):
            log.info("Host fasta file loaded")
            print("HPIpy: Host fasta file loaded")
        else:
            log.error("Host input file is not in fasta format, check input file format")
            log.info("HPIpy stopped!")
            print("HPIpy: Host input file is not in fasta format, check input file format")
            print("HPIpy stopped.")
            exit()

        # Pathogen
        if util.is_fasta(pathogenDecomp):
            log.info("Pathogen fasta file loaded")
            print("HPIpy: Pathogen fasta file loaded")
        else:
            log.error("Pathogen input file is not in fasta format, check input file format.")
            log.info("HPIpy stopped!")
            print("HPIpy: Pathogen input file is not in fasta format, check input file format")
            print("HPIpy stopped.")
            exit()


        ## Validate protein sequences
        # Host sequences
        if util.is_protein(hostDecomp):
            log.info("Host fasta file passed protein sequence check")
            print("HPIpy: Host FASTA file passed protein sequence check")
        else:
            log.error("Host file does not contain protein sequences. Check input sequences")
            log.info("HPIpy stopped!")
            print("HPIpy: Host FASTA file does not contain protein sequences")
            print("HPIpy stopped.")
            exit()

        # Pathogen sequences
        if util.is_protein(pathogenDecomp):
            log.info("Pathogen fasta file passed protein sequence check")
            print("HPIpy: Pathogen FASTA file passed protein sequence check")
        else:
            log.error("Pathogen file does not contain protein sequences. Check input sequences")
            log.info("HPIpy stopped!")
            print("HPIpy: Pathogen FASTA file does not contain protein sequences")
            print("HPIpy stopped.")
            exit()

        log.info("Input FASTA files successfully validated")
        print("HPIpy: Input FASTA files successfully validated")


        ## Sequence homology
        log.info(f'Sequence clustering directory created successfully :: {abs_outputdir}/Clustering')
        print(f'HPIpy: Sequence clustering directory created successfully :: {abs_outputdir}/Clustering')

        # Host
        hostClusters = util.cdhit(hostDecomp, outputdir, hostFastaName)
        log.info(f"Number of clusters identified in '{hostFastaName}': {hostClusters}")
        
        # Pathogen
        pathogeClusters = util.cdhit(pathogenDecomp, outputdir, pathogenFastaName)
        log.info(f"Number of clusters identified in '{pathogenFastaName}': {pathogeClusters}")


        ## Extract computational method from args
        if isinstance(args.computation, list):
            args.computation = args.computation[0] if len(args.computation) == 1 and isinstance(args.computation[0], list) else args.computation

        for method in args.computation:
            ## Interolog model
            if method == 'interolog':
                log.info(f"Performing BLASTp alignments for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")
                print(f"HPIpy: Performing BLASTp alignments for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")
                
                base_blast = blast.BLAST(model = args.model, hostFile = hostDecomp, pathogenFile = pathogenDecomp, use_slurm = args.use_slurm, 
                                            num_threads = args.num_threads, log=log, genome_pool = args.genome_pool, phyloEvalue = args.phyloEvalue,
                                            threshold = args.phyloThreshold)
                util.downloadData(data_directory = pkgData)
                base_blast.indexBlastDB(data_directory = pkgData, outputdir = outputdir)
                
                # check blast index files
                if not util.check_blastDB_files(data_directory = pkgData):
                    log.error("makeblastdb failed to generate all required files")
                    print("HPIpy Error: makeblastdb failed to generate all required files")
                    exit()
                else:
                    log.info("BLAST index files generated successfully. Proceeding...")
                    print("HPIpy: BLAST index files generated successfully. Proceeding...")
                
                # Run blast
                if args.use_slurm:
                    blastJobs = base_blast.executeBLAST(data_directory = pkgData, outputdir = outputdir)
                    log.info(f"BLASTp job for host and pathogen sequences submitted with SLURM job IDs: {blastJobs}")
                    print(f"HPIpy: BLAST job for host and pathogen submitted with SLURM job IDs: {blastJobs}")

                    log.info("Waiting for all the SLURM jobs to finish")
                    print("HPIpy: Waiting for all the SLURM jobs to finish...")

                    while any(util.is_job_running(job_id) for job_id in blastJobs):
                        time.sleep(60)

                    log.info("BLAST alignments completed successfully")
                    print("HPIpy: BLAST alignments completed successfully")
                else:
                    base_blast.executeBLAST(data_directory = pkgData, outputdir = outputdir)
                    log.info("BLAST alignments completed successfully")
                    print("HPIpy: BLAST alignments completed successfully")


            ## Domain model
            if method == 'domain':
                log.info(f"Extracting significant domains from '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")
                print(f"HPIpy: Extracting significant domains from '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")

                base_hmmer = hmmer.HMMER(model = args.model, hostFile = hostDecomp, pathogenFile = pathogenDecomp, use_slurm = args.use_slurm, 
                                            num_threads = args.num_threads, log=log)
                pfam_dir = base_hmmer.pfam(data_directory = pkgData, outputdir = outputdir)
                util.downloadData(data_directory = pkgData)

                # check hmmpress file
                if not util.check_hmm_files(pfam_dir):
                    log.error("'hmmpress' failed to generate the required files. Delete the existing 'PfamDB' directory from 'HPIpy_data' directory and run the program again")
                    print("HPIpy: Error: 'hmmpress' failed to generate the required files. Delete the existing 'PfamDB' directory from 'HPIpy_data' directory and run the program again")
                    exit()

                # execute hmmscan
                if args.use_slurm:
                    hmmerJobs = base_hmmer.hmmer(PfamDir = pfam_dir, outputdir = outputdir)
                    log.info(f"HMMER job for host and pathogen sequences submitted with SLURM job IDs: {hmmerJobs}")
                    print(f"HPIpy: HMMER job for host and pathogen submitted with SLURM job IDs: {hmmerJobs}")

                    log.info("Waiting for all the SLURM jobs to finish")
                    print("HPIpy: Waiting for all the SLURM jobs to finish...")

                    while any(util.is_job_running(job_id) for job_id in hmmerJobs):
                        time.sleep(60)                
                else:
                    base_hmmer.hmmer(PfamDir = pfam_dir, outputdir = outputdir)

                # Modifying hmmscan results
                os.system("awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18}' " + outputdir + "/Domains/" + hostFastaName.split(".")[0] + "_domains.txt > " + outputdir + "/Domains/" + hostFastaName.split(".")[0] + "_iddi_domains.txt")
                os.system(f"cp {outputdir}/Domains/{hostFastaName.split('.')[0]}_iddi_domains.txt {outputdir}/Domains/{hostFastaName.split('.')[0]}_domine_domains.txt")
                os.system(f"cp {outputdir}/Domains/{hostFastaName.split('.')[0]}_iddi_domains.txt {outputdir}/Domains/{hostFastaName.split('.')[0]}_did3_domains.txt")
                os.system("awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18}' " + outputdir + "/Domains/" + pathogenFastaName.split(".")[0] + "_domains.txt > " + outputdir + "/Domains/" + pathogenFastaName.split(".")[0] + "_iddi_domains.txt")
                os.system(f"cp {outputdir}/Domains/{pathogenFastaName.split('.')[0]}_iddi_domains.txt {outputdir}/Domains/{pathogenFastaName.split('.')[0]}_domine_domains.txt")
                os.system(f"cp {outputdir}/Domains/{pathogenFastaName.split('.')[0]}_iddi_domains.txt {outputdir}/Domains/{pathogenFastaName.split('.')[0]}_did3_domains.txt")

                log.info(f"Domain extraction from {hostFastaName.split('.')[0]} and {pathogenFastaName.split('.')[0]} sequences completed successfully")
                print(f"HPIpy: Domain extraction from {hostFastaName.split('.')[0]} and {pathogenFastaName.split('.')[0]} sequences completed successfully")
                        

            ## Phylogenetic profiling model
            if method == 'phyloProfiling':
                log.info(f"Performing DIAMOND BLASTp alignments for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")
                print(f"HPIpy: Performing DIAMOND BLASTp alignments for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")

                base_blast = blast.BLAST(model = args.model, hostFile = hostDecomp, pathogenFile = pathogenDecomp, use_slurm = args.use_slurm, 
                                            num_threads = args.num_threads, log=log, genome_pool = args.genome_pool, phyloEvalue = args.phyloEvalue,
                                            threshold = args.phyloThreshold)
                util.downloadData(data_directory = pkgData)
                ngenome_number = int(re.search(r'\d+', args.genome_pool).group())
                _, poolList = base_blast.indexPhyloDB(data_directory = pkgData, ngenome = ngenome_number, outputdir = outputdir)
                                
                # Run DIAMOND BLASTp
                if args.use_slurm:
                    phyloBlastJobs = base_blast.phylo_blast(ngenome = ngenome_number, phyloPoolList = poolList, data_directory = pkgData, outputdir = outputdir)
                    log.info(f"DIAMOND BLASTp job for host and pathogen sequences submitted with SLURM job IDs: {phyloBlastJobs}")
                    print(f"HPIpy: DIAMOND BLAST job for host and pathogen submitted with SLURM job IDs: {phyloBlastJobs}")
                else:
                    base_blast.phylo_blast(ngenome = ngenome_number, phyloPoolList = poolList, data_directory = pkgData, outputdir = outputdir)
                    log.info("DIAMOND BLAST alignments completed successfully")
                    print("HPIpy: DIAMOND BLAST alignments completed successfully")


            ## GO similarity
            if args.interproscan:
                log.info(f"Running Interproscan for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")
                print(f"HPIpy: Running Interproscan for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")

                base_go = goSimilarity.GOSimilarity(log = log)

                # Host interproscan
                host_interpro_file = base_go.run_interproscan(inputFasta = hostFastaName, outputdir = outputdir, data_directory = pkgData)
                log.info(f"Extracting GO terms from Interproscan output for '{hostFastaName.split('.')[0]}'")
                print(f"HPIpy: Extracting GO terms from Interproscan output for '{hostFastaName.split('.')[0]}'")
                base_go.extractGOTerms(host_interpro_file)

                # Pathogen interproscan
                pathogen_interpro_file = base_go.run_interproscan(inputFasta = pathogenFastaName, outputdir = outputdir, data_directory = pkgData)
                log.info(f"Extracting GO terms from Interproscan output for '{pathogenFastaName.split('.')[0]}'")
                print(f"HPIpy: Extracting GO terms from Interproscan output for '{pathogenFastaName.split('.')[0]}'")
                base_go.extractGOTerms(pathogen_interpro_file)
            else:
                pass
        

        ## Generate tables
        for method in args.computation:
            base_tables = tables.TABLES(model = args.model, hostFile = hostFastaName, pathogenFile = pathogenFastaName)

            # Interolog tables
            if method == 'interolog':
                print(f"HPIpy: Generating local SQL tables for {method}-based predictions")
                log.info(f"Generating local SQL tables for {method}-based predictions")

                databases_Interolog = {"humanVirus": ["biogrid", "dip", "hpidb", "intact", "mint", "virhostnet"], 
                                    "plantPathogen": ["biogrid", "dip", "hpidb", "intact", "mint"], 
                                    "animalPathogen": ["biogrid", "dip", "hpidb", "intact", "mint", "virhostnet"], 
                                    "humanBacteria": ["biogrid", "dip", "hpidb", "intact", "mint"]}
                dbsPPI = databases_Interolog[args.model]
                
                for i in dbsPPI:
                    base_tables.create_interaction_table(file = f"{pkgData}/{args.model}/dbs/{i}_interactions.txt", 
                                                                sep = "\t", table = i, db = f"{pkgData}/interactions.db")
                    base_tables.create_table_blast(file = f"{outputdir}/Alignment/Interolog/{hostFastaName.split('.')[0]}_{i}_blast.txt", 
                                                sep = "\t", table = i, db = f"{pkgData}/{hostFastaName.split('.')[0]}.db")
                    base_tables.create_table_blast(file = f"{outputdir}/Alignment/Interolog/{pathogenFastaName.split('.')[0]}_{i}_blast.txt", 
                                                sep = "\t", table = i, db = f"{pkgData}/{pathogenFastaName.split('.')[0]}.db")

                    log.info(f"Table for '{i}' database created")
                    print(f"HPIpy: Table for '{i}' database created")

            # Domain tables
            if method == 'domain':
                print(f"HPIpy: Generating local SQL tables for {method}-based predictions")
                log.info(f"Generating local SQL tables for {method}-based predictions")
            
                databases_Domain = {"humanVirus": ["did3", "domine", "iddi"], 
                                    "plantPathogen": ["did3", "domine", "iddi"], 
                                    "animalPathogen": ["did3", "domine", "iddi"], 
                                    "humanBacteria": ["did3", "domine", "iddi"]}
                dbsDDI = databases_Domain[args.model]
                
                for j in dbsDDI:
                    base_tables.create_interaction_table(file = f"{pkgData}/domainDBs/{j}_interactions.txt", 
                                                                sep = "\t", table = j, db = f"{pkgData}/interactions.db")
                    base_tables.create_table_hmm(file = f"{outputdir}/Domains/{hostFastaName.split('.')[0]}_{j}_domains.txt", 
                                                sep = " ", table = j, db = f"{pkgData}/{hostFastaName.split('.')[0]}.db")
                    base_tables.create_table_hmm(file = f"{outputdir}/Domains/{pathogenFastaName.split('.')[0]}_{j}_domains.txt", 
                                                sep = " ", table = j, db = f"{pkgData}/{pathogenFastaName.split('.')[0]}.db")
                    
                    log.info(f"Table for '{j}' database created")
                    print(f"HPIpy: Table for '{j}' database created")

            # SQL databases generation check
            if method == 'interolog' or method == 'domain':
                if not util.check_db_files(pkgData):
                    log.error(f"Local SQL databases failed to generate for '{method}' method. Delete the existing '.db' files and run the program again")
                    print(f"HPIpy: Error: Local SQL databases failed to generate '{method}' method. Delete the existing '.db' files and run the program again")
                    exit()
                else:
                    log.info(f"Local tables generated successfully for '{method}' method")
                    print(f"HPIpy: Tables generated successfully for '{method}' method")

        # Human annotations tables
        if args.model == "humanVirus" or args.model == "humanBacteria":
            util.downloadAnnot(data_directory = pkgData)
            annots = ['human_localization', 'human_pathways', 'human_drugs', 'human_chembl']

            for k in annots:
                try:
                    base_tables.table_human_annotations(file=f"{pkgData}/annotations/{k}.csv", sep=";", table=k, db=f"{pkgData}/annotations.db")
                except Exception as err:
                    print(f"HPIpy: Error generating annotation database: {err}")
                    log.error(f"Error generating annotation database: {err}")
            
            util.del_directory(f'{pkgData}/annotations')
        else:
            pass
        

        ## Predict interactions between host and pathogen
        log.info(f"Predicting interactions between '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}'")    
        print(f"HPIpy: Predicting interactions between '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}'")
        
        base_ppi = predict.PREDICT(model = args.model, hostFile = hostFastaName, hostFilePath = hostInputFile, pathogenFile = pathogenFastaName, 
                                   pathogenFilePath = pathogenInputFile, hostGOFile = hostGOPath if hostGOPath else None, 
                                   pathogenGOFile = pathogenGOPath if pathogenGOPath else None, log=log)
        base_ppi.predictInteractions(outputdir = outputdir, data_directory = pkgData, computation_methods = args.computation)

        log.info(f"Interactions between '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' successfully predicted")
        print(f"HPIpy: Interactions between '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' successfully predicted")


    ## Resume predictions
    if args.resume_ppis:
        print("***** HPIpy: A tool for host-pathogen interaction prediction *****")
        print("\n -------------------------------------------------------------------------")
        print("  HPIpy: A tool for host-pathogen protein-protein interactions prediction ")
        print(" -------------------------------------------------------------------------\n")

        # Results directory
        outputdir = util.make_directory(args.outputdir)
        abs_outputdir = os.path.abspath(outputdir)
        pkgData = os.path.join(os.path.dirname(outputdir), "HPIpy_data")

        # Log directory
        logsDir = util.make_directory(os.path.join(outputdir, "logs"))
        abs_logdir = os.path.abspath(logsDir)
        log = logger.logHPIpy(logdir = outputdir, mode = "a")

        # Executed command
        command = ""
        for i, arg in enumerate(sys.argv):
            command += arg + " "
        log.info(f"Command executed: {command.strip()}")

        log.info("Resuming pipeline for protein-protein interactions prediction with new parameters")
        print("HPIpy: Resuming pipeline for protein-protein interactions prediction with new parameters")
        log.info(f"Successfully created results directory :: {abs_outputdir}")
        print(f"HPIpy: Successfully created output directory :: {abs_outputdir}")
        log.info(f"Successfully created logs directory :: {abs_logdir}")
        print(f"HPIpy: Successfully created logs directory :: {abs_logdir}")
        
        # Processing input files
        hostInputFile = args.host
        hostDecomp = util.decompress_file(hostInputFile)
        hostFastaPath = Path(hostDecomp)
        hostFastaName = hostFastaPath.name
        pathogenInputFile = args.pathogen
        pathogenDecomp = util.decompress_file(pathogenInputFile)
        pathogenFastaPath = Path(pathogenDecomp)
        pathogenFastaName = pathogenFastaPath.name

        # GO terms files
        hostGOPath = None
        pathogenGOPath = None
        if args.hostGOFile:
            hostGOFile = args.hostGOFile
            hostGOPath = Path(hostGOFile)
        if args.pathogenGOFile:
            pathogenGOFile = args.pathogenGOFile
            pathogenGOPath = Path(pathogenGOFile)
        
        # Parameters
        print("HPIpy: Parameters provided by the user for analysis:")
        log.info("Parameters provided by the user for analysis:")
        
        related_args = {
            'interIdentity': ['interCoverage', 'interEvalue', 'outputdir', 'network', 'num_threads'],
            'domHostEvalue': ['domPathogenEvalue', 'outputdir', 'network', 'num_threads'],
            'phyloThreshold': ['genome_pool', 'phyloIdentity', 'phyloCoverage', 'phyloEvalue', 'outputdir', 'network', 'num_threads'],
            'pathogenGOFile': ['goSimThreshold', 'interproscan', 'outputdir', 'network', 'num_threads'],
            'interproscan': ['hostGOFile', 'pathogenGOFile', 'goSimThreshold', 'outputdir', 'network', 'num_threads']
        }

        printed_args = set()
        for arg in vars(args):
            if getattr(args, arg) != parser.get_default(arg):
                print(f"    > {arg.replace('_', ' ').title()}: '{getattr(args, arg)}'")
                log.info(f"    > {arg.replace('_', ' ').title()}: '{getattr(args, arg)}'")
                printed_args.add(arg)
                if arg in related_args:
                    for related_arg in related_args[arg]:
                        if related_arg not in printed_args:
                            print(f"    > {related_arg.replace('_', ' ').title()}: '{getattr(args, related_arg)}'")
                            log.info(f"    > {related_arg.replace('_', ' ').title()}: '{getattr(args, related_arg)}'")
                            printed_args.add(related_arg)

        
        ## Validate fasta format
        log.info("Input FASTA files validation started")
        
        # Host
        if util.is_fasta(hostDecomp):
            log.info("Host fasta file loaded")
            print("HPIpy: Host fasta file loaded")
        else:
            log.error("Host input file is not in fasta format, check input file format")
            log.info("HPIpy stopped!")
            print("HPIpy: Host input file is not in fasta format, check input file format")
            print("HPIpy stopped.")
            exit()

        # Pathogen
        if util.is_fasta(pathogenDecomp): # pathogenFasta
            log.info("Pathogen fasta file loaded")
            print("HPIpy: Pathogen fasta file loaded")
        else:
            log.error("Pathogen input file is not in fasta format, check input file format.")
            log.info("HPIpy stopped!")
            print("HPIpy: Pathogen input file is not in fasta format, check input file format")
            print("HPIpy stopped.")
            exit()
    

        ## Validate protein sequences
        # Host sequences
        if util.is_protein(hostDecomp):
            log.info("Host fasta file passed protein sequence check")
            print("HPIpy: Host FASTA file passed protein sequence check")
        else:
            log.error("Host file does not contain protein sequences. Check input sequences")
            log.info("HPIpy stopped!")
            print("HPIpy: Host FASTA file does not contain protein sequences")
            print("HPIpy stopped.")
            exit()

        # Pathogen sequences
        if util.is_protein(pathogenDecomp):
            log.info("Pathogen fasta file passed protein sequence check")
            print("HPIpy: Pathogen FASTA file passed protein sequence check")
        else:
            log.error("Pathogen file does not contain protein sequences. Check input sequences")
            log.info("HPIpy stopped!")
            print("HPIpy: Pathogen FASTA file does not contain protein sequences")
            print("HPIpy stopped.")
            exit()

        log.info("Input FASTA files successfully validated")
        print("HPIpy: Input FASTA files successfully validated")


        ## Sequence homology
        log.info(f'Sequence clustering directory created successfully :: {abs_outputdir}/Clustering')
        print(f'HPIpy: Sequence clustering directory created successfully :: {abs_outputdir}/Clustering')

        # Host
        hostClusters = util.cdhit(hostDecomp, outputdir, hostFastaName)
        log.info(f"Number of clusters identified in '{hostFastaName}': {hostClusters}")
        
        # Pathogen
        pathogeClusters = util.cdhit(pathogenDecomp, outputdir, pathogenFastaName)
        log.info(f"Number of clusters identified in '{pathogenFastaName}': {pathogeClusters}")

        ## Extract computational method from args
        if isinstance(args.computation, list):
            args.computation = args.computation[0] if len(args.computation) == 1 and isinstance(args.computation[0], list) else args.computation

        for method in args.computation:
            ## Phylogenetic profiling model
            if method == 'phyloProfiling':
                log.info(f"Performing DIAMOND BLASTp alignments for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")
                print(f"HPIpy: Performing DIAMOND BLASTp alignments for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")
                
                base_blast = blast.BLAST(model = args.model, hostFile = hostDecomp, pathogenFile = pathogenDecomp, use_slurm = args.use_slurm, 
                                            num_threads = args.num_threads, log=log, genome_pool = args.genome_pool, phyloEvalue = args.phyloEvalue,
                                            threshold = args.phyloThreshold)
                util.downloadData(data_directory = pkgData)
                ngenome_number = int(re.search(r'\d+', args.genome_pool).group())
                _, poolList = base_blast.indexPhyloDB(data_directory = pkgData, ngenome = ngenome_number, outputdir = outputdir)
                        
                # Run DIAMOND BLASTp
                if args.use_slurm:
                    phyloBlastJobs = base_blast.phylo_blast(ngenome = ngenome_number, phyloPoolList = poolList, data_directory = pkgData, outputdir = outputdir)
                    log.info(f"DIAMOND BLASTp job for host and pathogen sequences submitted with SLURM job IDs: {phyloBlastJobs}")
                    print(f"HPIpy: DIAMOND BLAST job for host and pathogen submitted with SLURM job IDs: {phyloBlastJobs}")
                else:
                    base_blast.phylo_blast(ngenome = ngenome_number, phyloPoolList = poolList, data_directory = pkgData, outputdir = outputdir)
                    log.info("DIAMOND BLAST alignments completed successfully")
                    print("HPIpy: DIAMOND BLAST alignments completed successfully")


            ## GO similarity
            if args.interproscan:
                if not os.path.exists(f"{pkgData}/Interproscan_output/{hostFastaName.split('.')[0]}_go.csv") and \
                    not os.path.exists(f"{pkgData}/Interproscan_output/{pathogenFastaName.split('.')[0]}_go.csv"):

                    log.info(f"Running Interproscan for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")
                    print(f"HPIpy: Running Interproscan for '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' protein sequences")

                    base_go = goSimilarity.GOSimilarity(log = log)

                    # Host interproscan
                    host_interpro_file = base_go.run_interproscan(inputFasta = hostFastaName, outputdir = outputdir, data_directory = pkgData)
                    log.info(f"Extracting GO terms from Interproscan output for '{hostFastaName.split('.')[0]}'")
                    print(f"HPIpy: Extracting GO terms from Interproscan output for '{hostFastaName.split('.')[0]}'")
                    base_go.extractGOTerms(host_interpro_file)

                    # Pathogen interproscan
                    pathogen_interpro_file = base_go.run_interproscan(inputFasta = pathogenFastaName, outputdir = outputdir, data_directory = pkgData)
                    log.info(f"Extracting GO terms from Interproscan output for '{pathogenFastaName.split('.')[0]}'")
                    print(f"HPIpy: Extracting GO terms from Interproscan output for '{pathogenFastaName.split('.')[0]}'")
                    base_go.extractGOTerms(pathogen_interpro_file)
                else:
                    pass
            else:
                pass

            # SQL databases generation check
            if method == 'interolog' or method == 'domain':
                if not util.check_db_files(pkgData):
                    log.error("Local SQL databases failed to generate. Delete the existing '.db' files and re-run the program from start")
                    print("HPIpy: Error: Local SQL databases failed to generate. Delete the existing '.db' files and re-run the program from start")
                    exit()
                else:
                    log.info(f"Required SQL databases exist for '{method}' method, proceeding...")
                    print(f"HPIpy: Required SQL databases exist for '{method}' method, proceeding...")

        ## Predict interactions between host and pathogen
        log.info(f"Predicting interactions between '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}'")    
        print(f"HPIpy: Predicting interactions between '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}'")
        
        base_ppi = predict.PREDICT(model = args.model, hostFile = hostFastaName, hostFilePath = hostInputFile, pathogenFile = pathogenFastaName, 
                                   pathogenFilePath = pathogenInputFile, hostGOFile = hostGOPath if hostGOPath else None, 
                                   pathogenGOFile = pathogenGOPath if pathogenGOPath else None, log=log)
        base_ppi.predictInteractions(outputdir = outputdir, data_directory = pkgData, computation_methods = args.computation)

        log.info(f"Protein-protein interactions between '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' successfully predicted")
        print(f"HPIpy: Protein-protein interactions between '{hostFastaName.split('.')[0]}' and '{pathogenFastaName.split('.')[0]}' successfully predicted")

if __name__ == '__main__':
    main()