import argparse
from colorama import init, Fore, Style
init(autoreset=True)

parser = argparse.ArgumentParser(prog='hpipy', 
                                description=f'''{Fore.GREEN + Style.BRIGHT}
  hpipy: A package to predict host-microbe protein-protein interactions
  ---------------------------------------------------------------------

  To obtain more information about the package, visit: https://kaabil.net/hpipy/
                                {Style.RESET_ALL}''',
                                usage='python3 -m %(prog)s [options]',
                                epilog=f'''

{Fore.GREEN + Style.BRIGHT}Examples:{Style.RESET_ALL}
    {Fore.GREEN + Style.BRIGHT}## Single computational method:{Style.RESET_ALL}
    {Fore.WHITE}python3 -m hpipy 
        --host hostProteins.fasta 
        --pathogen pathogenProteins.fasta 
        --computation interolog 
        --model humanVirus 
        --num_threads 10 
        --network {Style.RESET_ALL}

    {Fore.GREEN + Style.BRIGHT}## Multiple computational methods:{Style.RESET_ALL}
    {Fore.WHITE}python3 -m hpipy 
        --host hostProteins.fasta 
        --pathogen pathogenProteins.fasta 
        --computation interolog domain phyloProfiling gosim 
        --model humanVirus 
        --num_threads 10 
        --hostGOFile hostGOTerms.csv 
        --pathogenGOFile pathogenGOTerms.csv 
        --network {Style.RESET_ALL}
{Fore.BLUE + Style.BRIGHT}
--------------------------------------------------------------------------------------
| Developed by Raghav Kataria | KAABiL (https://kaabil.net/) | Utah State University |
--------------------------------------------------------------------------------------
                                {Style.RESET_ALL}''', 

                                formatter_class = argparse.RawDescriptionHelpFormatter)

# Required arguments
requiredArgs = parser.add_argument_group(f"{Fore.RED + Style.BRIGHT}Required arguments{Style.RESET_ALL}")

requiredArgs.add_argument('--host', type=str, required=True, 
                    help='Protein sequences of host species (formats accepted: .fasta, .fasta.gz, .fasta.zip, .fa, .fa.gz, .fa.zip, .faa, .faa.gz, .faa.zip)', 
                    metavar = '<hostFastaFile>')

requiredArgs.add_argument('--pathogen', type=str, required=True, 
                    help='Protein sequences of pathogen species (formats accepted: .fasta, .fasta.gz, .fasta.zip, .fa, .fa.gz, .fa.zip, .faa, .faa.gz, .faa.zip)', 
                    metavar = '<pathogenFastaFile>')

requiredArgs.add_argument('--computation', type=str, required=True, nargs='+', 
                    help= "Computational method(s) to be implemented for the analysis. Provide a space-separated list; Available methods: interolog, domain, phyloProfiling, gosim",
                    metavar = '<method>')

requiredArgs.add_argument('--model', type=str, required=True, 
                    help= "Host-pathogen model to be implemented for the analysis. Available models: plantPathogen, animalPathogen, humanVirus, humanBacteria",
                    metavar = '<choose model>')


# Optional arguments
optional_args = parser.add_argument_group(f"{Fore.RED + Style.BRIGHT}Optional Arguments{Style.RESET_ALL}")

pkgVersion = 1.0
optional_args.add_argument("--version", action="version", version = f"hpipy v{pkgVersion}", 
                        help='Display package version and exit')

optional_args.add_argument("--outputdir", type=str, default="HPIpy_results", 
                        help='Directory where output files will be written; default: "HPIpy_results"',
                        metavar = '<results directory>')

optional_args.add_argument("--use_slurm", action="store_true", 
                        help="To run jobs using SLURM job scheduler")

optional_args.add_argument("--slurm_account", action="store_true", 
                        help="To run SLURM jobs on a specific account on the cluster")

optional_args.add_argument("--network", action="store_true", 
                        help="To perform network analysis for the predicted interactions. This step will take more time based on the number of predicted interactions.")

optional_args.add_argument("--num_threads", type=int, default=4, 
                        help="Number of threads to be used; default: 4",
                        metavar = '<int>')

optional_args.add_argument("--seq_homology", type=float, default=1.0, 
                        help="Sequence identity for CD-HIT (0.1 to 1.0); default: 1.0 (100 percent)",
                        metavar = '<float>')

optional_args.add_argument("--resume_ppis", action="store_true",
                        help="To predict interactions using different parameters without running the whole pipeline. Provide 'interproscan' option\
                            if it was used before, although it will no be executed again if interproscan's output files already exist")


# Interolog model
interolog_args = parser.add_argument_group(f"{Fore.RED + Style.BRIGHT}Interolog model prediction arguments (optional){Style.RESET_ALL}")

interolog_args.add_argument("--interIdentity", type=int, default=50, 
                        help='Sequence identity to filter BLAST alignments; default: 50',
                        metavar = '<int>')

interolog_args.add_argument("--interCoverage", type=int, default=50, 
                        help='Sequence coverage to filter BLAST alignments; default: 50',
                        metavar = '<int>')

interolog_args.add_argument("--interEvalue", type=str, default="1e-05", 
                        help='e-value to filter BLAST alignments; default: 1e-05',
                        metavar = '<str>')


# Domain model
domain_args = parser.add_argument_group(f"{Fore.RED + Style.BRIGHT}Domain model prediction arguments (optional){Style.RESET_ALL}")

domain_args.add_argument("--domHostEvalue", type=str, 
                        help='e-value to filter host HMMER output; default is based on the selected model',
                        metavar = '<str>')

domain_args.add_argument("--domPathogenEvalue", type=str, 
                        help='e-value to filter pathogen HMMER output; default is based on the selected model',
                        metavar = '<str>')


# Phylogenetic profiling model
phylo_args = parser.add_argument_group(f"{Fore.RED + Style.BRIGHT}Phylogenetic profiling model prediction arguments (optional){Style.RESET_ALL}")

phylo_args.add_argument('--genome_pool', type=str, default='BC20', 
                    help= "Genome pool to be used for phylogenetic profiling model. Available pools: UP82, BC20, protPhylo490; default: BC20",
                    metavar = '<choose pool>')

phylo_args.add_argument("--phyloEvalue", type=str, default="1e-05", nargs='+', 
                        help='e-value to filter DIAMOND BLAST alignments; default: 1e-05',
                        metavar = '<str>')

phylo_args.add_argument("--phyloIdentity", type=int, default=50, nargs='+', 
                        help='Sequence identity to filter DIAMOND BLAST alignments; default: 50',
                        metavar = '<int>')

phylo_args.add_argument("--phyloCoverage", type=int, default=50, nargs='+', 
                        help='Sequence coverage to filter DIAMOND BLAST alignments; default: 50',
                        metavar = '<int>')

phylo_args.add_argument("--phyloThreshold", type=float, default=0.9, 
                        help="Threshold value to filter predicted interactions based on phylogenetic distance (0.1 to 1.0); default: 0.9",
                        metavar = '<float>')


# GO similarity model
gosem_args = parser.add_argument_group(f"{Fore.RED + Style.BRIGHT}GO semantic similarity model prediction arguments (optional){Style.RESET_ALL}")

gosem_args.add_argument("--interproscan", action="store_true", 
                        help='To run InterProScan locally to obtain GO terms for GO similarity model. InterProScan should be installed locally')

gosem_args.add_argument('--hostGOFile', type=str, 
                    help= 'Comma- or tab-separated file containing GO terms of host proteins',
                    metavar = '<Host GO terms file>')

gosem_args.add_argument('--pathogenGOFile', type=str, 
                    help= 'Comma- or tab-separated file containing GO terms of pathogen proteins',
                    metavar = '<Pathogen GO terms file>')

gosem_args.add_argument('--go_combine', type=str, default='BMA', 
                    help= 'Method to combine GO similarity scores. Available methods: max, avg, rcmax and BMA; default: BMA',
                    metavar = '<choose combine method>')

gosem_args.add_argument("--goSimThreshold", type=float, default=0.9, 
                        help="Threshold value (0.1 to 1.0) for GO semantic similarity scores to filter predicted interactions; default: 0.9",
                        metavar = '<float>')

args = parser.parse_args()