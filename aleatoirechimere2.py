from Bio import SeqIO
import pandas as pd
import statistics
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from pathlib import Path
import os.path
from random import *




seq1 = input('Entrez le chemin de la seq1: ')
seq2 = input('Entrez le chemin de la seq2: ')

path1 = Path(seq1)
nom1=path1.name
print(nom1)


path2 = Path(seq2)
nom2=path2.name
print(nom2)


#on commence par définir les séquences d'interet
seqRecord1 = SeqIO.read(seq1, 'fasta')
seqRecord2 = SeqIO.read(seq2, 'fasta')



#on verifie que ce sont bien elles 
print(seqRecord1)
print(seqRecord2)
longueur1=len(seqRecord1)
longueur2=len(seqRecord2)
print(longueur1)
print(longueur2)

#voici le doc contenant toutes les séqunces backbone
#df=pd.read_csv('chimere_backbone.txt',sep="\t")
#df['index'] = df.index
#print (df)

#voici le doc contenant toutes les séqunces homologues backbone triées et filtrées sans les zéros
dfh1=pd.read_csv('chimere_backbone_liste_homologue1.txt',sep="\t")
dfh1['index'] = dfh1.index
print (dfh1)

#voici le doc contenant toutes les séqunces homologues backbone triées et filtrées sans les zéros
dfh2=pd.read_csv('chimere_backbone_liste_homologue2.txt',sep="\t")
dfh2['index'] = dfh2.index
print (dfh2)

#on considère la premiere colonne de la data frame comme une liste pour en faire la médiane

list_seq0_leftend1 = dfh1['seq0_leftend'].tolist()
list_seq0_leftend2 = dfh2['seq0_leftend'].tolist()



listenom=[]
listerapportd=[]
listerapportr=[]
taux_recombinaison=[]

j=0

for start in list_seq0_leftend1:
    j=j+1
    print(start)
    reqd_Index = dfh1[dfh1['seq0_leftend']==start].index.tolist()
    #recupéerer index de la sequence 
    print(reqd_Index)
    x = dfh1[dfh1['seq0_leftend']==start].seq1_leftend.tolist()
    #récuperer valeur de seq1leftendpour seq0=start
    p=statistics.median(x)
    #transforme la liste en variable 
    print(p)
    
    ale=randint(p, longueur2)
    
    seqdonneuse=seqRecord2[p:ale]
    longueur4=len(seqdonneuse)
    
    seqreceveuse=seqRecord1[0:start]+seqRecord1[ale:]
    longueur3=len(seqreceveuse)
    
    seqfinale=seqRecord1[0:start]+seqRecord2[p:ale]+seqRecord1[ale:]
    
    longueur5=len(seqfinale)

    rapportreceveuse=(longueur3/longueur5)*100
    listerapportr.append(rapportreceveuse)
    rapportdonneuse=(longueur4/longueur5)*100
    listerapportd.append(rapportdonneuse)
    taux_recombinaison.append(rapportdonneuse)
    #calcul le rapport de chaque parent

    

    
    sequence = seqfinale.seq
    gcrapport=GC(sequence)
    #donne le taux de gc de chaque sequence
    listenom.append(f'chimere{j}')
    titre=f">Chimere|{nom1}|#{j}\n"

    
    fichier1 = open(f'chimere{nom1}-{j}.fasta', "a")
    fichier1.write(f'{titre}')
    fichier1.write (f'{sequence}')
    fichier1.close()
    fichier2 = open(f'chimere{nom1}-{j}-infos.txt', "a")
    fichier2.write(f'{titre}')
    fichier2.write (f'Pourcentage de Receveuse:{rapportreceveuse}%\n')
    fichier2.write (f'Pourcentage de Donneuse:{rapportdonneuse}%\n')
    fichier2.write (f'Pourcentage de GC:{gcrapport}%\n')
    fichier2.write (f'Debut recombinaison: {p} \n')
    fichier2.write (f'Fin recombinaison : {ale} \n')
    fichier2.close()
    
    
data = {'Nom':listenom,'Pourcentage receveuse':listerapportr,'Pourcentage donneuse':listerapportd,'Taux de recombinaison':taux_recombinaison}
dfrap1 = pd.DataFrame(data)
dfrap1.to_csv(f'chimere{nom1}-{j}Resumé.csv')







listenom2=[]
listerapportd2=[]
listerapportr2=[]
taux_recombinaison2=[]
g=0
for start in list_seq0_leftend2:
    g=g+1
    print(start)
    reqd_Index = dfh2[dfh2['seq0_leftend']==start].index.tolist()
    #recupéerer index de la sequence 
    print(reqd_Index)
    x = dfh2[dfh2['seq0_leftend']==start].seq1_leftend.tolist()
    #récuperer valeur de seq1leftendpour seq0=start
    p=statistics.median(x)
    #transforme la liste en variable 
    print(p)
    ale=randint(p, longueur1)
    
    seqdonneuse=seqRecord1[p:ale]
    longueur4=len(seqdonneuse)
    
    seqreceveuse=seqRecord2[0:start]+seqRecord2[ale:]
    longueur3=len(seqreceveuse)
    
    seqfinale=seqRecord2[0:start]+seqRecord1[p:ale]+seqRecord2[ale:]
    
    longueur5=len(seqfinale)

    rapportreceveuse=(longueur3/longueur5)*100
    listerapportr2.append(rapportreceveuse)
    rapportdonneuse=(longueur4/longueur5)*100
    listerapportd2.append(rapportdonneuse)
    taux_recombinaison2.append(rapportdonneuse)
    #calcul le rapport de chaque parent
    
    

    sequence = seqfinale.seq
    gcrapport=GC(sequence)
    #donne le taux de gc de chaque sequence
    
    listenom2.append(f'chimere{g}')
    titre=f">Chimere|{nom2}|#{g}\n"
    fichier1 = open(f'chimere{nom2}-{g}.fasta', "a")
    fichier1.write(f'{titre}')
    fichier1.write (f'{sequence}')
    fichier1.close()
    
    fichier2 = open(f"chimere{nom2}-{g}-infos.txt", "a")
    fichier2.write(f'{titre}')
    fichier2.write (f'Pourcentage de Receveuse:{rapportreceveuse}%\n')
    fichier2.write (f'Pourcentage de Donneuse:{rapportdonneuse}%\n')
    fichier2.write (f'Pourcentage de GC: {gcrapport}% \n')
    fichier2.write (f'Debut recombinaison: {p} \n')
    fichier2.write (f'Fin recombinaison : {ale} \n')
    fichier2.close()

    
    
data2 = {'Nom':listenom,'Pourcentage receveuse':listerapportr,'Pourcentage donneuse':listerapportd,'Taux de recombinaison':taux_recombinaison}
dfrap2 = pd.DataFrame(data2)
dfrap2.to_csv(f'chimere{nom2}-Resumé.csv')


