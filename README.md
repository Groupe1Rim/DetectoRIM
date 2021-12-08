

  # DetectoRIM(Detection of Recombination In streptoMyces)

![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)

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
To create chimeria you have to use Chimerator.py
First you have to aligned the 2 parent sequences on the Mauve tool (http://darlinglab.org/mauve/mauve.html)

### Input files

Chimerator.py needs .backbone output file from a Mauve alignment and the two aligned sequences

.backbone output file must be in text format and sequences files must be in Fasta format

### Usage

To create chimeria use Chimerator.py :

``python Chimerator.py``

### Output files

A successful Chimerator.py run will generate differents files :

- A number of chimeras in fasta file equivalent to the number of homologous zones between each parent
- A text file with the percentage of membership for each parent for each chimera
- A .csv file containing for each chimera the percentage of membership in each parent

### Warning

Be careful, if the two sequences studied are genetically close, this can result in the creation of a large number of chimeras. So make sure you have enough memory on your hard drive. 

## Create final sequene from consensus sequence

To create a final mix sequence you have to use Consensusator.py
First you have to aligned the 2 consensus sequences on the Mauve tool (http://darlinglab.org/mauve/mauve.html)

### Input files

Consensusator.py needs .backbone output file from a Mauve alignment and the two aligned sequences

.backbone output file must be in text format and sequences files must be in Fasta format

### Usage

To create final sequence use Consensusator.py :

``python Consensusator.py``

### Output files

A successful Consensusator.py run will generate :
- A final sequence in fasta format

## Calculates and represents the recombination

To calculates and represents the recombination you have to use Recombinator.py
First you have to aligned the 2 prents sequences and the recombinant sequence on the Mauve tool (http://darlinglab.org/mauve/mauve.html)

### Input files
Recombinator.py needs SNP output file from a Mauve alignment and the three aligned sequences

SNP output file must be in text format and sequences files must be in Fasta format

### Usage

To create calculates recombination in the recombinant sequence and have a graphic representation of the recombinant genome use Recombinator.py :

``python Recombinator.py``

Recombinator.py will ask which are the two parent sequences. For the first parent sequence you have to write the name of the first sequence aligned in Mauve and for the second parent sequence you have to write the name of the second sequence aligned in Mauve. You must then define a window in which SNPs blocks will be created in order to define from which parent each part of the genome of the recombinant sequence comes from. The higher the density of SNPs on the genome, the better it is to use a small window(like 5 SNPs.)


### Output files

A successful Recombinator.py run will generate :

- An excel document containing all the SNPs that concern each of the 3 sequences named "total_SNPs.xlsx".With for each SNP its position on the genome of the recombinant, its origin (To which parent it belongs) as well as the distance which separates it from each SNP (upstream and downstream of the genome)
- An excel document containing only the SNPs blocks that are created following the user defined window named "Sorted_SNPs.xlsx". .With for each SNPs block  its position on the genome of the recombinant (start position and end position), its origin (To which parent it belongs) and the size of each block (Distance between each block).
- A png image showing which parent each part of the recombinant sequence genome comes from. As shown below: 

![alt text](https://github.com/Groupe1Rim/DetectoRIM/blob/main/2_S2102finalpetit.png?raw=true)


## Authors

* **Alizée Bueche** 
* **Clémentine Isembart** 
* **Julie Charles** 
* **Louis L'Hôte** 
* **Pauline Gascht**




    
     


