# HPIpy  

<img src="https://kaabil.net/hpipy/downloads/homepage.png" alt="HPIpy" width="1000" height="200">

HPIpy is a Python-based standalone package that leverages various sequence-based computational models for the prediction of protein-protein interactions (PPIs) between host and pathogen species.  

<details open><summary><b>Table of Contents</b></summary>  

  - [Computational approaches](#computation)
  - [Host-pathogen models](#ppi-models)  
  - [Download Miniconda](#miniconda)  
  - [Download the package](#package)  
  - [Usage](#usage)  
  - [Contact Us](#contact)

</details>  


## Computational approaches <a name="computation"></a>  
Various sequence-based computational approaches implemented in HPIpy:  
- Interolog mapping  
- Domain-based  
- Phylogenetic profiling  
- Gene ontology (GO) semantic similarity


## Host-pathogen interaction models <a name="ppi-models"></a>  
HPIpy supports the following pathosystems for protein-protein interactions prediction.  

| Model | Host | Pathogen(s) |
|------------|-------|--------------|
|`humanVirus`| Human | Human-related virus |
|`humanBacteria`| Human | Human-related bacteria |
|`animalPathogen`| All animals | Virus, fungi, bacteria |
|`plantPathogen`| All plants | Virus, fungi, bacteria |


## Download and Install Miniconda <a name="miniconda"></a>    
If not installed, download and install Miniconda:  
- Download the latest version from [Miniconda](https://docs.anaconda.com/free/miniconda/)  
- To install, execute the command:
```
bash Miniconda3-latest-Linux-x86_64.sh
```  


## Download HPIpy <a name="package"></a>  
Download the HPIpy package using one of the below options (1 or 2):  
1. Clone this repository:
```
git clone https://github.com/usubioinfo/hpipy.git
```  
2. Obtain package using link:
```
wget https://kaabil.net/hpipy/downloads/hpipy.tar.gz
tar -xvzf hpipy.tar.gz
```  

To create the conda environment for HPIpy:  
```
cd hpipy
conda env create -f environment.yml
conda activate hpipy
```


## Usage <a name="usage"></a>  
To view the "Help" section of HPIpy, run:
```
python3 -m hpipy --help
```

Basic usage of HPIpy:  
```
python3 -m hpipy --host exampleData/hostProteins.fasta --pathogen exampleData/pathogenProteins.fasta --computation interolog --model humanVirus
```

Apart from required arguments (as mentioned above), HPIpy contains several options (see package help), which can be used based on the requirement. Example for advanced analysis:  
```
python3 -m hpipy --host exampleData/hostProteins.fasta --pathogen exampleData/pathogenProteins.fasta \
  --computation interolog \
  --model humanVirus \
  --seq_homology 0.8 \
  --num_threads 20 \
  --network \
  --outputdir samplePPIs \
  --interIdentity 40 60 \
  --interEvalue 1e-15 1e-20 \
  --interCoverage 60
```

If you want to predict the interactions based on different parameters for prediction, use the `--resume_ppis` option. The program will not run all the steps but only execute PPIs prediction step.
```
python3 -m hpipy --host exampleData/hostProteins.fasta --pathogen exampleData/pathogenProteins.fasta \
  --computation interolog \
  --model humanVirus \
  --seq_homology 0.8 \
  --num_threads 20 \
  --network \
  --outputdir samplePPIs \
  --interIdentity 70 80 \
  --interEvalue 1e-25 \
  --resume_ppis
```


## Contact Us <a name="contact"></a>  
For any queries, contact us at [bioinfo@kaabil.net](mailto:raghav.kataria@usu.edu).  
