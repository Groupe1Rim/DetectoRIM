

  # DetectoRIM(Detection of Recombination In streptoMyces)

Tool to detect recombination in Streptomyces species very close genetically

## Contents

## Installation

### Required dependencies
      
Python version-->3.0 or more

Python modules:

      -pandas
      -statistics
      -openpyxl
      -dna_features_viewer
      -Biopython
      
      
### Installation procedure

Use the git command ``git clone https://github.com/Groupe1Rim/DetectoRIM.git``
You can then use each script independently


## Create chimeria
### Input files

Chimerator.py needs .backbone output file from a Mauve alignment and the two aligned sequences

.backbone output file must be in text format and sequences files must be in Fasta format

### Usage

To create chimeria use chimerator.py :

``python chimerator.py``

### Output files

A successful Chimerator.py run will generate differents files :

- A number of chimeras equivalent to the number of homologous zones between each parent
- A .csv file containing for each chimera the percentage of membership in each parent







## Authors

* **Alizée Bueche** 
* **Clémentine Isembart** 
* **Julie Charles** 
* **Louis L'Hôte** 
* **Pauline Gascht**




    
     


