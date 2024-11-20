from .args_parse import *
from .blast import *
from .hmmer import *
from .interolog import *
from .domain import *
from .logger import *
from . import network
from pathlib import Path
from . import utility as util
import itertools
from . import tables
from .phyloProfile import *
from .goSimilarity import *
import glob
import re

class PREDICT:
    """
    Class for predicting protein-protein interactions
    """

    def __init__(self, model = None, hostFile = None, hostFilePath = None, pathogenFile = None, pathogenFilePath = None, hostGOFile = None, pathogenGOFile = None, log = None):
        """
        Initialize the PREDICT class.

        :param model: Model used for prediction.
        :param hostFile: File containing host protein sequences.
        :param hostFilePath: Complete path to host FASTA file.
        :param pathogenFile: File containing pathogen protein sequences.
        :param pathogenFilePath: Complete path to pathogen FASTA file.
        :param log: Log file for logging messages.

        :type model: str
        :type hostFile: str
        :type hostFilePath: str
        :type pathogenFile: str
        :type pathogenFilePath: str
        :type log: str
        """
        
        self.model = model
        self.hostFile = hostFile
        hostFastaPath = Path(self.hostFile)
        self.hostFasta = hostFastaPath.name
        self.hPath = hostFilePath
        self.pathogenFile = pathogenFile
        pathogenFastaPath = Path(self.pathogenFile)
        self.pathogenFasta = pathogenFastaPath.name
        self.pPath = pathogenFilePath
        self.interologModel = Interolog()
        self.domainModel = Domain()
        self.phyloModel = PhyloProfiling()
        self.goSimModel = GOSimilarity()
        self.log = log
        self.phyloIdentity = args.phyloIdentity
        self.phyloCoverage = args.phyloCoverage
        self.phyloThreshold = args.phyloThreshold
        self.hostGOFile = hostGOFile
        self.pathogenGOFile = pathogenGOFile
        self.go_similarity = 'Wang'
        self.go_combine = args.go_combine
        self.goSimThreshold = args.goSimThreshold
        self.computation = args.computation
        self.identity = args.interIdentity
        self.coverage = args.interCoverage
        self.evalue = args.interEvalue
        self.databasesInterolog = {"humanVirus": ["biogrid", "dip", "hpidb", "intact", "mint", "virhostnet"], 
                                   'plantPathogen': ["biogrid", "dip", "hpidb", "intact", "mint"], 
                                   'animalPathogen': ["biogrid", "dip", "hpidb", "intact", "mint", "virhostnet"],
                                   "humanBacteria": ["biogrid", "dip", "hpidb", "intact", "mint"]}
        self.databasesDomain = {"humanVirus": ["did3", "domine", "iddi"], 
                                'plantPathogen': ["did3", "domine", "iddi"],
                                'animalPathogen': ["did3", "domine", "iddi"],
                                'humanBacteria': ["did3", "domine", "iddi"]}
        self.ddiModel = {"humanVirus": [1e-23, 1e-01], 
                         'plantPathogen': [1e-23, 1e-17],
                         'animalPathogen': [1e-23, 1e-17],
                         'humanBacteria': [1e-23, 1e-18]}
        self.hostEvalue = self.ddiModel[self.model][0]
        self.pathogenEvalue = self.ddiModel[self.model][1]

        if args.domHostEvalue:
            self.hostEvalue = args.domHostEvalue
        if args.domPathogenEvalue:
            self.pathogenEvalue = args.domPathogenEvalue      

    
    ## Interactions prediction
    def predictInteractions(self, outputdir, data_directory, computation_methods):
        """
        This function predicts protein-protein interactions (PPIs) based on various computational methods.

        :param outputdir: Directory to save the output files.
        :param data_directory: Directory containing data files necessary for each prediction method.
        :param computation_methods: List of computational methods for PPI prediction.

        :type outputdir: str
        :type data_directory: str
        :type computation_methods: list of str

        :return: Summary statistics for PPIs predictions.
        :rtype: None
        """

        ppi_output = os.path.join(outputdir, "Predictions")
        os.makedirs(ppi_output, exist_ok=True)
        abs_ppis_outdir = os.path.abspath(ppi_output)
        self.log.info(f"Output directory for storing predicted interactions created successfully :: {abs_ppis_outdir}")
        print(f"HPIpy: Output directory for storing predicted interactions created successfully :: {abs_ppis_outdir}")

        if isinstance(self.computation, list):
            self.computation = self.computation[0] if len(self.computation) == 1 and isinstance(self.computation[0], list) else self.computation
        
        # Prediction methods
        method_to_function = {
            'domain': self.predict_domain,
            'phyloProfiling': self.predict_phylo,
            'gosim': self.predict_go,
            'interolog': self.predict_interolog
        }
        
        selected_results = {}
        concatenated_results = []

        with open(f"{ppi_output}/Prediction_stats.txt", "a") as stat:
            stat.write(f">> Protein-protein interactions prediction results <<\n\n")
        
        for method in computation_methods:
            if method in method_to_function:
                result = method_to_function[method](outputdir, ppi_output, data_directory)
                ppis, methodStats, methodDir = result
                selected_results[method] = ppis     # consensus
                concatenated_results.append(result) # combined
                
                # Write stats
                with open(f"{ppi_output}/Prediction_stats.txt", "a") as stat:
                    stat.write(f"Method: {method.upper()}\n")
                    for col in methodStats.columns:
                        stat.write(f"   {col:<20}|  {methodStats[col].iloc[0]:>5}\n")
                    stat.write("\n")

                # Extract proteins IDs & their sequences
                host = ppis[['Host']].drop_duplicates()
                util.extract_sequences(self.hPath, host, 'Host', methodDir, method)
                pathogen = ppis[['Pathogen']].drop_duplicates()
                util.extract_sequences(self.pPath, pathogen, 'Pathogen', methodDir, method)
                
                # Network
                self.network(ppis, methodDir)

                # Human annotations
                if self.model == "humanVirus" or args.model == "humanBacteria":
                    self.human_annot(ppis[['Host']], data_directory, methodDir)
                    
                    drugs_file = f"{methodDir}/human_annotations/human_drugs.txt"
                    if os.path.exists(drugs_file):
                        with open(drugs_file, "a") as stat:
                            stat.write(f"\nNote: Drug IDs and their common names are obtained from DrugBank (https://go.drugbank.com/).\n")
                    else:
                        pass

        # Combined PPIs
        if len(computation_methods) > 1:
            combined_stats = self.combined_predictions(ppi_output, concatenated_results)
            with open(f"{ppi_output}/Prediction_stats.txt", "a") as stat:
                stat.write(f"COMBINED INTERACTIONS:\n")

                for col in combined_stats.columns:
                    stat.write(f"   {col:<20}|  {combined_stats[col].iloc[0]:>5}\n")
                stat.write("\n")
        
        # Consensus PPIs
        if len(selected_results) > 1:
            consensus_stats = self.consensus_predictions(out_dir=ppi_output, **selected_results)
            with open(f"{ppi_output}/Prediction_stats.txt", "a") as stat:
                stat.write(f"CONSENSUS INTERACTIONS:\n")
                for _, row in consensus_stats.iterrows():
                    stat.write(f"   Consensus_Pair      |  {row['Consensus_Pair']}\n")
                    stat.write(f"   Interactions        |  {row['Interactions']}\n")
                    stat.write(f"   Host_Proteins       |  {row['Host_Proteins']}\n")
                    stat.write(f"   Pathogen_Proteins   |  {row['Pathogen_Proteins']}\n\n")
        else:
            pass
        
        return
    

    ## Combined PPIs
    def combined_predictions(self, out_dir, results):
        """
        This function merges interactions and statistics from multiple prediction methods.

        :param out_dir: Directory to save the combined output files.
        :param results: List of tuples containing results from individual computation methods.

        :type out_dir: str
        :type results: list of tuples

        :return: Summary statistics for the combined protein-protein interactions.
        :rtype: pd.DataFrame
        """

        self.log.info("Predicting combined interactions from the provided computational methods")
        print("HPIpy: Predicting combined interactions from the provided computational methods")

        out_dir_combined = os.makedirs(os.path.join(out_dir, "Combined_PPIs"), exist_ok=True) or os.path.join(out_dir, "Combined_PPIs")
        statCombo = []
        abs_combined_outdir = os.path.abspath(out_dir_combined)
        interaction_dfs = [result[0] for result in results]
        combined_interactions = pd.concat(interaction_dfs, ignore_index=True).drop_duplicates(subset = ["Host", "Pathogen"], keep='first')
        combined_interactions = combined_interactions[['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'Interaction_Database', 'InteractionType']]
        combined_interactions.to_csv(f"{out_dir_combined}/Combined_PPIs.csv", index=False)

        # stat
        statCombinedDict = {"Interactions": [combined_interactions.shape[0]], 
                            "Host_Proteins": [combined_interactions['Host'].nunique()],
                            "Pathogen_Proteins": [combined_interactions['Pathogen'].nunique()]}
        statCombined_df = pd.DataFrame(statCombinedDict, index=None)
        statCombo.append(statCombined_df)
        statCombinedFinal = pd.concat(statCombo)

        self.log.info(f"Combined interactions successfully saved to :: {abs_combined_outdir}")
        print(f"HPIpy: Combined interactions successfully saved to :: {abs_combined_outdir}")
        
        return statCombinedFinal
    

    ## Consensus PPIs
    def consensus_predictions(self, out_dir, **model_PPIs):
        """
        Function to find consensus protein-protein interactions (PPIs) across multiple models.

        :param out_dir: Directory to save the consensus output files.
        :param model_PPIs: Keyword arguments where each key is a model name and each value is a DataFrame of predicted PPIs.

        :type out_dir: str
        :type model_PPIs: dict

        :return: Summary statistics for consensus predictions.
        :rtype: pd.DataFrame
        """

        self.log.info("Predicting consensus interactions between the provided computational methods")
        print("HPIpy: Predicting consensus interactions between the provided computational methods")

        out_dir_consensus = os.makedirs(os.path.join(out_dir, "Consensus_PPIs"), exist_ok=True) or os.path.join(out_dir, "Consensus_PPIs")
        statConsensus = []
        abs_consensus_outdir = os.path.abspath(out_dir_consensus)
        consensus_dfs = {}

        model_columns = {
            'interolog': ['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'DetectionMethod', 'ConfidenceScore', 'PubMedID', 'Interaction_Database'],
            'domain': ['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'Interaction_Database'],
            'phyloProfiling': ['Host', 'Pathogen', 'Score', 'Interaction_Database'],
            'gosim': ['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'Score', 'Interaction_Database']
        }
                
        model_names = list(model_PPIs.keys())
        for i, model1 in enumerate(model_names):
            for model2 in model_names[i+1:]:
                # Merge models
                consensus = model_PPIs[model1].merge(model_PPIs[model2], on=['Host', 'Pathogen'], how='outer', indicator=True)
                consensus = consensus[consensus['_merge'] == 'both'].drop_duplicates(subset=['Host', 'Pathogen']).drop(columns='_merge')
                columns1 = [col for col in model_columns[model1] if col in consensus.columns]
                columns2 = [col for col in model_columns[model2] if col in consensus.columns and col not in ['Host', 'Pathogen']]
                consensus_selected = consensus[columns1 + columns2].copy()
                renamed_columns = {col: f"{col}_{model1}" for col in columns1 if col not in ['Host', 'Pathogen']}
                renamed_columns.update({col: f"{col}_{model2}" for col in columns2})
                consensus_selected = consensus_selected.rename(columns=renamed_columns)

                # Column for model source
                consensus_selected['Method_Source'] = f"{model1}_{model2}"
                file_name = f"{model1}_{model2}_consensus_PPIs.csv"
                consensus_selected.to_csv(f"{out_dir_consensus}/{file_name}", sep="\t", index=False)
                consensus_dfs[file_name] = consensus_selected
                
                # stats
                statConsensus.append({
                    "Consensus_Pair": f"{model1}_{model2}",
                    "Interactions": consensus_selected.shape[0],
                    "Host_Proteins": consensus_selected['Host'].nunique(),
                    "Pathogen_Proteins": consensus_selected['Pathogen'].nunique()
                })

        # Consensus from all models
        if len(consensus_dfs) > 2:
            consensus_all = pd.concat(consensus_dfs.values(), ignore_index=True).drop_duplicates(subset=['Host', 'Pathogen'])
            consensus_all = consensus_all[['Host', 'Pathogen', 'Method_Source']]
            consensus_all.to_csv(f"{out_dir_consensus}/Consensus_PPIs_all_methods.csv", sep="\t", index=False)

            statConsensus.append({
                "Consensus_Pair": "All_Consensus",
                "Interactions": consensus_all.shape[0],
                "Host_Proteins": consensus_all['Host'].nunique(),
                "Pathogen_Proteins": consensus_all['Pathogen'].nunique()
                })
        statConsensus_df = pd.DataFrame(statConsensus)

        self.log.info(f"Consensus interactions successfully saved to :: {abs_consensus_outdir}")
        print(f"HPIpy: Consensus interactions successfully saved to :: {abs_consensus_outdir}")
        
        return statConsensus_df


    ## Domain-based PPIs
    def predict_domain(self, resultDir, outputdir, data_directory):
        """
        Predict domain-based protein-protein interactions (PPIs)

        :param outputdir: Directory to save the domain-based prediction output files.
        :param data_directory: Directory containing interaction data files.

        :type outputdir: str
        :type data_directory: str

        :return: A tuple containing the domain-based PPIs, statistics, and the path to domain-based results directory.
        :rtype: tuple(pd.DataFrame, pd.DataFrame, str)
        """
        
        self.log.info("Predicting domain-based interactions")
        print("HPIpy: Predicting domain-based interactions")

        dbsDDI = self.databasesDomain[self.model]
        domain_dir = os.path.join(outputdir, "Domain-based")
        os.makedirs(domain_dir, exist_ok=True)
        abs_domain_outdir = os.path.abspath(domain_dir)
        self.log.info(f"Domain-based PPIs directory created successfully :: {abs_domain_outdir}")
        print(f"HPIpy: Domain-based PPIs directory created successfully :: {abs_domain_outdir}")
        
        domainPredicted = []
        domainAnnotAB = []
        statDomain = []
        
        for j in dbsDDI:
            host_ddi, ddi_df_host = self.domainModel.filter_hmm_domain(f"{data_directory}/{self.hostFasta.split('.')[0]}.db", j, self.hostEvalue, host = True)
            pathogen_ddi, ddi_df_pathogen = self.domainModel.filter_hmm_domain(f"{data_directory}/{self.pathogenFasta.split('.')[0]}.db", j, self.pathogenEvalue, host = False)
            try:
                ddiAB = self.domainModel.search_domain_id(f"{data_directory}/interactions.db", host_ddi, pathogen_ddi, j)
                ddi_db = self.domainModel.df_merge_domain(ddiAB, ddi_df_pathogen, ddi_df_host).drop_duplicates(subset = ["Host", "Pathogen"])
                ddi_db['Interaction_Database'] = j
                ddi_db_ppis = ddi_db[['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'Interaction_Database']]
                domainPredicted.append(ddi_db_ppis)
                ddi_db_ppis.to_csv(f"{domain_dir}/{j}_PPI.csv", sep = "\t", index = False)

                ddi_annotationA = ddi_db[['Interactor_A', 'PfamName_A', 'InterproID_A', 'InterproName_A', 'GO_IDs_A', 'GO_name_A', 'PDB_A']].drop_duplicates()
                ddi_annotationA.columns = ['InteractorID', 'PfamName', 'InterproID', 'InterproName', 'GO_IDs', 'GO_name', 'PDB']
                domainAnnotAB.append(ddi_annotationA)

                ddi_annotationB = ddi_db[['Interactor_B', 'PfamName_B', 'InterproID_B', 'InterproName_B', 'GO_IDs_B', 'GO_name_B', 'PDB_B']].drop_duplicates()
                ddi_annotationB.columns = ['InteractorID', 'PfamName', 'InterproID', 'InterproName', 'GO_IDs', 'GO_name', 'PDB']
                domainAnnotAB.append(ddi_annotationB)
            except Exception:
                ddi_db = pd.DataFrame()

        # Concatenating Domain results
        domainPPIs = pd.concat(domainPredicted)
        domainPPIs['InteractionType'] = 'Domain'
        domainUniquePPI = domainPPIs.drop_duplicates(subset = ["Host", "Pathogen"])
        domainUniquePPI.to_csv(f"{domain_dir}/Domain_PPIs.csv", sep = "\t", index = False)
        domainAnnotations = pd.concat(domainAnnotAB, ignore_index=True).drop_duplicates()
        domainAnnotations.to_csv(f"{domain_dir}/Domain_Annotations.txt", sep = "\t", index = False)

        # Domain-based statistics
        statDomainDict = {"Interactions": [domainUniquePPI.shape[0]], 
                        "Host_Proteins": [domainUniquePPI['Host'].nunique()], 
                        "Pathogen_Proteins": [domainUniquePPI['Pathogen'].nunique()]}
        statDomain_df = pd.DataFrame(statDomainDict, index=None)
        statDomain.append(statDomain_df)
        statDomainFinal = pd.concat(statDomain)
        util.del_directory(f'{data_directory}/domainDBs')
        
        return domainUniquePPI, statDomainFinal, domain_dir


    ## Phylogenetic profiling-based PPIs
    def predict_phylo(self, resultDir, outputdir, data_directory):
        """
        Predict phylogenetic profile-based protein-protein interactions

        :param outputdir: Directory where the phylogenetic profiling-based output files will be saved.
        :param data_directory: Directory containing necessary input data files.

        :type outputdir: str
        :type data_directory: str

        :return: A tuple containing the phylogenetic profiling-based PPIs, statistics, and the path to phylogenetic profiling-based results directory.
        :rtype: tuple(pd.DataFrame, pd.DataFrame, str)
        """

        self.log.info("Predicting phylogenetic profiling-based interactions")
        print("HPIpy: Predicting phylogenetic profiling-based interactions")

        inputClusters = f"{resultDir}/Clustering"
        phylo_dir = os.path.join(outputdir, "PhyloProfiling-based")
        os.makedirs(phylo_dir, exist_ok=True)
        abs_phylo_outdir = os.path.abspath(phylo_dir)
        self.log.info(f"Phylogenetic profiling-based PPIs directory created successfully :: {abs_phylo_outdir}")
        print(f"HPIpy: Phylogenetic profiling-based PPIs directory created successfully :: {abs_phylo_outdir}")

        phyloPredicted = []
        statPhylo = []
        ngenome_number = int(re.search(r'\d+', args.genome_pool).group())
        nullPool = '0' * ngenome_number
        host_blast_files = glob.glob(f"{resultDir}/Alignment/PhyloProfiling/*_{self.hostFasta.split('.')[0]}_blastOut.txt")
        pathogen_blast_files = glob.glob(f"{resultDir}/Alignment/PhyloProfiling/*_{self.pathogenFasta.split('.')[0]}_blastOut.txt")
        
        try:
            hostIDs, pathogenIDs, numberHost, _, _, _ = self.phyloModel.extractIDs(f"{inputClusters}/{self.hostFasta}", f"{inputClusters}/{self.pathogenFasta}")
            pattern_host = self.phyloModel.process_files(host_blast_files, hostIDs, self.phyloIdentity, self.phyloCoverage)
            pattern_pathogen = self.phyloModel.process_files(pathogen_blast_files, pathogenIDs, self.phyloIdentity, self.phyloCoverage)
            phylo_ppis = self.phyloModel.phylo_ppis(numberHost, pattern_host, pattern_pathogen, nullPool, hostIDs, pathogenIDs, ngenome_number, self.phyloThreshold)
            phyloPredicted.append(phylo_ppis)
        except Exception:
            phylo_ppis = pd.DataFrame()
        
        if phyloPredicted:
            phyloPPIs = pd.concat(phyloPredicted)
            phyloPPIs['GenomePool'] = args.genome_pool
            phyloPPIs['InteractionType'] = 'PhyloProfiling'
            phyloUniquePPI = phyloPPIs.drop_duplicates(subset = ["Host", "Pathogen"])
            phyloUniquePPI.to_csv(f"{phylo_dir}/PhyloProfiling_PPIs.csv", sep = "\t", index = False)
            phyloUniquePPI['Interactor_A'] = '-'
            phyloUniquePPI['Interactor_B'] = '-'
            phyloUniquePPI = phyloUniquePPI[['Host', 'Pathogen', 'Score', 'Interactor_A', 'Interactor_B', 'GenomePool', 'InteractionType']]
            phyloUniquePPI.columns = ['Host', 'Pathogen', 'Score', 'Interactor_A', 'Interactor_B', 'Interaction_Database', 'InteractionType']
        else:
            phyloUniquePPI = pd.DataFrame(columns = ['Host', 'Pathogen', 'Score', 'Interactor_A', 'Interactor_B', 'Interaction_Database', 'InteractionType'])
            phyloUniquePPI.to_csv(f"{phylo_dir}/PhyloProfiling_PPIs.csv", sep = "\t", index = False)
        
        # Phylogenetic profiling-based PPIs statistics
        statPhyloDict = {"Interactions": [phyloUniquePPI.shape[0]], 
                        "Host_Proteins": [phyloUniquePPI['Host'].nunique()], 
                        "Pathogen_Proteins": [phyloUniquePPI['Pathogen'].nunique()]}
        statPhylo_df = pd.DataFrame(statPhyloDict, index=None)
        statPhylo.append(statPhylo_df)
        statPhyloFinal = pd.concat(statPhylo)
        util.del_directory(f'{data_directory}/{args.genome_pool}')

        return phyloUniquePPI, statPhyloFinal, phylo_dir


    ## GO similarity-based PPIs
    def predict_go(self, resultDir, outputdir, data_directory):
        """
        Predict Gene Ontology-based protein-protein interactions using semantic similarity

        :param outputdir: Directory where the GO similarity-based output files will be saved.
        :param data_directory: Directory containing necessary input data files.

        :type outputdir: str
        :type data_directory: str

        :return: A tuple containing the GO similarity-based PPIs, statistics, and the path to GO similarity-based results directory.
        :rtype: tuple(pd.DataFrame, pd.DataFrame, str)
        """

        self.log.info("Predicting GO similarity-based interactions")
        print("HPIpy: Predicting GO similarity-based interactions")
        
        gosim_dir = os.path.join(outputdir, "GOSimilarity-based")
        os.makedirs(gosim_dir, exist_ok=True)
        abs_go_outdir = os.path.abspath(gosim_dir)
        self.log.info(f"GO similarity-based PPIs directory created successfully :: {abs_go_outdir}")
        print(f"HPIpy: GO similarity-based PPIs directory created successfully :: {abs_go_outdir}")

        statGOSim = []

        # Interproscan
        if args.interproscan:
            host_go = f"{data_directory}/Interproscan_output/{self.hostFasta.split('.')[0]}_go.csv"
            pathogen_go = f"{data_directory}/Interproscan_output/{self.pathogenFasta.split('.')[0]}_go.csv"
            if host_go is None or pathogen_go is None:
                self.log.error(f"GO terms could not be extracted from InterProScan output")
                print(f"HPIpy: GO terms could not be extracted from InterProScan output")
                goUniquePPI = pd.DataFrame(columns=['Host', 'Pathogen', 'Score', 'Interactor_A', 'Interactor_B', 'Interaction_Database', 'InteractionType'])
        
        # GO terms files
        else:
            host_go = None
            pathogen_go = None
            if self.hostGOFile and self.pathogenGOFile:
                try:
                    host_go = self.goSimModel.readGOFile(self.hostGOFile)
                    pathogen_go = self.goSimModel.readGOFile(self.pathogenGOFile)
                except Exception as e:
                    self.log.error(f"One or both of the host and pathogen GO term files are missing.")
                    print(f"HPIpy: One or both of the host and pathogen GO term files are missing.")
                    goUniquePPI = pd.DataFrame(columns=['Host', 'Pathogen', 'Score', 'Interactor_A', 'Interactor_B', 'Interaction_Database', 'InteractionType'])
            else:
                pass

        # Run GOSemSim
        goUniquePPI = pd.DataFrame(columns=['Host', 'Pathogen', 'Score', 'Interactor_A', 'Interactor_B', 'Interaction_Database', 'InteractionType'])

        if host_go is not None and pathogen_go is not None:
            semData = self.goSimModel.initializeSemData(ontology='BP')
            try:
                if semData is not None:
                    go_ppis = self.goSimModel.predictGOPPIs(host_go, pathogen_go, semData, self.go_similarity, self.go_combine, self.goSimThreshold)
                else:
                    self.log.info("GO Semantic data could not be initialized. Proceeding without GO similarity prediction.")
                    print(f"HPIpy: Warning: GO Semantic data could not be initialized. Proceeding without GO similarity prediction.")
            except Exception:
                go_ppis = pd.DataFrame(columns=['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'Score'])

            go_ppis['GOSimMethod'] = self.go_similarity
            go_ppis['InteractionType'] = 'GOSimilarity'
            goUniquePPI = go_ppis.drop_duplicates(subset = ["Host", "Pathogen"])
            goUniquePPI.to_csv(f"{gosim_dir}/GOSimilarity_PPIs.csv", sep = "\t", index = False)
            goUniquePPI.columns = ['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'Score', 'Interaction_Database', 'InteractionType']
            goUniquePPI = goUniquePPI[['Host', 'Pathogen', 'Score', 'Interactor_A', 'Interactor_B', 'Interaction_Database', 'InteractionType']]

        # GO similarity-based PPIs statistics
        statGOSimDict = {"Interactions": [goUniquePPI.shape[0]], 
                        "Host_Proteins": [goUniquePPI['Host'].nunique()], 
                        "Pathogen_Proteins": [goUniquePPI['Pathogen'].nunique()]}
        statGOSim_df = pd.DataFrame(statGOSimDict, index=None)
        statGOSim.append(statGOSim_df)
        statGOSimFinal = pd.concat(statGOSim)

        return goUniquePPI, statGOSimFinal, gosim_dir


    ## Interolog-based PPIs  
    def predict_interolog(self, resultDir, outputdir, data_directory):
        """
        Predict protein-protein interactions using the interolog-based approach

        :param outputdir: Directory where the interolog-based output files will be saved.
        :param data_directory: Directory containing necessary input data files.

        :type outputdir: str
        :type data_directory: str

        :return: A tuple containing the interolog-based unique PPIs, statistics, and the path to interolog-based results directory.
        :rtype: tuple(pd.DataFrame, pd.DataFrame, str)
        """

        self.log.info("Predicting interolog-based interactions")
        print("HPIpy: Predicting interolog-based interactions")

        dbsPPI = self.databasesInterolog[args.model]
        count = 0
        interolog_dir = os.path.join(outputdir, "Interolog-based")
        os.makedirs(interolog_dir, exist_ok=True)
        abs_interolog_outdir = os.path.abspath(interolog_dir)
        self.log.info(f"Interolog-based PPIs directory created successfully :: {abs_interolog_outdir}")
        print(f"HPIpy: Interolog-based PPIs directory created successfully :: {abs_interolog_outdir}")

        statInterolog = []
        interologPredicted = []
        interologAnnotAB = []
        count += 1

        for i in dbsPPI:
            host_ppidb, ppidb_df_host = self.interologModel.filter_blast(f"{data_directory}/{self.hostFasta.split('.')[0]}.db", i, self.identity, self.evalue, self.coverage, host = True)
            pathogen_ppidb, ppidb_df_pathogen = self.interologModel.filter_blast(f"{data_directory}/{self.pathogenFasta.split('.')[0]}.db", i, self.identity, self.evalue, self.coverage, host = False)
            try:
                ppiAB = self.interologModel.search_id(f"{data_directory}/interactions.db", host_ppidb, pathogen_ppidb, i)
                ppi_db = self.interologModel.df_merge_interolog(ppiAB, ppidb_df_pathogen, ppidb_df_host).drop_duplicates(subset = ["Host", "Pathogen"])
                ppi_db['Interaction_Database'] = i
                ppi_db_ppis = ppi_db[['Host', 'Pathogen', 'Interactor_A', 'Interactor_B', 'DetectionMethod', 'ConfidenceScore', 'PubMedID', 'Interaction_Database']]
                interologPredicted.append(ppi_db_ppis)
                ppi_db_ppis.to_csv(f"{interolog_dir}/{i}_PPI.csv", sep = "\t", index = False)
                
                ppi_annotationA = ppi_db[['Interactor_A', 'EntrezGeneID_A', 'GO_IDs_A', 'GO_name_A', 'PDB_A']].drop_duplicates()
                ppi_annotationA.columns = ['InteractorID', 'EntrezGeneID', 'GO_IDs', 'GO_name', 'PDB']
                interologAnnotAB.append(ppi_annotationA)

                ppi_annotationB = ppi_db[['Interactor_B', 'EntrezGeneID_B', 'GO_IDs_B', 'GO_name_B', 'PDB_B']].drop_duplicates()
                ppi_annotationB.columns = ['InteractorID', 'EntrezGeneID', 'GO_IDs', 'GO_name', 'PDB']
                interologAnnotAB.append(ppi_annotationB)
            except Exception:
                ppi_db = pd.DataFrame()

        # Concatenating Interolog results
        interologPPIs = pd.concat(interologPredicted, ignore_index=True)
        interologPPIs['InteractionType'] = 'Interolog'
        interologUniquePPI = interologPPIs.drop_duplicates(subset = ["Host", "Pathogen"])
        interologUniquePPI.to_csv(f"{interolog_dir}/Interolog_PPIs.csv", sep = "\t", index = False)
        interologAnnotations = pd.concat(interologAnnotAB, ignore_index=True).drop_duplicates()
        interologAnnotations.to_csv(f"{interolog_dir}/Interolog_Annotations.txt", sep = "\t", index = False)

        # Interolog-based PPIs statistics
        statInterologDict = {"Interactions": [interologUniquePPI.shape[0]],
                             "Host_Proteins": [interologUniquePPI['Host'].nunique()],
                             "Pathogen_Proteins": [interologUniquePPI['Pathogen'].nunique()]}
        statInterolog_df = pd.DataFrame(statInterologDict, index=None)
        statInterolog.append(statInterolog_df)
        statInterologFinal = pd.concat(statInterolog)
        util.del_directory(f'{data_directory}/{args.model}')
        
        return interologUniquePPI, statInterologFinal, interolog_dir    

    
    ## Human annotations
    def human_annot(self, hostList, data_directory, out_dir):
        """
        Retrieve human protein annotations based on specified protein IDs.

        :param hostList: List of human protein IDs.
        :param data_directory: Directory containing the annotations database.
        :param out_dir: Directory where the human annotations will be saved.

        :type hostList: list[str]
        :type data_directory: str
        :type out_dir: str

        :return: None
        """

        out_dir_annot = os.makedirs(os.path.join(out_dir, "human_annotations"), exist_ok=True) or os.path.join(out_dir, "human_annotations")
        base_tables = tables.TABLES(model = self.model, hostFile = self.hostFasta.split('.')[0], pathogenFile = self.pathogenFasta.split('.')[0])
        
        if self.model == "humanVirus" or self.model == "humanBacteria":
            try:
                base_tables.query_human_annotations(proteinIDs = hostList, db = f"{data_directory}/annotations.db", output_directory = out_dir_annot)
            except Exception as err:
                print(f"HPIpy: No annotations found for human proteins: {err}")
                self.log.info(f"No annotations found for human proteins: {err}")
        else:
            pass


    ## Network analysis
    def network(self, ppis, outputdir):
        """
        Perform network analysis on protein-protein interaction (PPI) data and save the results.

        :param ppis: DataFrame containing protein-protein interactions to analyze.
        :param outputdir: Directory where the network analysis results will be saved.

        :type ppis: pd.DataFrame
        :type outputdir: str

        :return: None
        """
        
        if args.network:
            self.log.info(f"   Network analysis started")
            print(f"   Network analysis started")
            
            networkDir = f"{outputdir}/network_analysis/"
            os.makedirs(networkDir, exist_ok=True)

            base_network = network.ProteinInteractionNetwork(log=self.log)
            graph_obj = base_network.initiate_graph(interactions = ppis)
            base_network.calculate_protein_hubs(graph_object = graph_obj, networkOutDir = networkDir)
            base_network.calculate_betweenness_centrality(graph_object = graph_obj, networkOutDir = networkDir)
            base_network.calculate_degree_centrality(graph_object = graph_obj, networkOutDir = networkDir)
            base_network.calculate_closeness_centrality(graph_object = graph_obj, networkOutDir = networkDir)
            base_network.save_network(graph_object = graph_obj, networkOutDir = networkDir)

            self.log.info("   Network analysis completed")
            print("   Network analysis completed")
        else:
            pass