############################################################################Importation des modules############################################################################
import pandas as pd
import statistics
import openpyxl
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import time
from Bio import SeqFeature
from dna_features_viewer import (
    GraphicFeature,
    GraphicRecord,
    CircularGraphicRecord,
)
import matplotlib.pyplot as plt
############################################################################Importation des modules############################################################################
seqRecord = SeqIO.read('1_S210_consensus_S1A1-3_.fasta', 'fasta')
seql=len(seqRecord)


df=pd.read_csv('SNPs',sep="\t")

long=input('Entrez la taille en nucléotide de la zone : ')
Longueur=int(long)



for i in df.index:
    if df.sequence_3_PosInContg[i]==0:
        df.drop(i, inplace=True)




df.sort_values(by=['sequence_3_GenWidePos3'], inplace=True)


liste= df['sequence_3_GenWidePos3'].tolist()
l=len(liste)


#####on a éliminé les valeures nulles, on s'interessent mtn suelement aux snps
####on va mtn calculé la densité moyenne des snps
den=(l/seql)
densite=(seql/l)
dens=int(densite)


Feature=[]


lstart=[]
for i in range(seql):
    if i%Longueur == 0:
        lstart.append(i)
lstart.append(seql)

len2=len(lstart)
len3=len2-1
lcompteur=[]

for t in range(len(lstart)):
    g=t+1
    if g<=len2:
        compteur=0
        for y in range(len(liste)):
            if lstart[t]<=liste[y]<=lstart[g]:
                compteur=compteur+1
    #if compteur>testden:
        #Feature.append(GraphicFeature(start=lstart[t], end=lstart[g],color="#ffcccc"))
    #else:
        #Feature.append(GraphicFeature(start=lstart[t], end=lstart[g],color="#cffccc"))
        
        lcompteur.append(compteur)



    
data = {'Start':lstart,'SNPs':lcompteur}
dfinal = pd.DataFrame(data)
dfinal['index'] = dfinal.index
dfinal.to_excel(f'total_SNPs.xlsx')



moyatti=den*Longueur
moyatt=int(moyatti)




dlarge=dfinal.nlargest(10, ['SNPs'])

newsizestart=[]
newsizeend=[]
for i in dlarge.index:
    h=0
    f=0
    h=Longueur+(dlarge.SNPs[i]*Longueur)
    f=h+Longueur
    newsizestart.append(h)
    newsizeend.append(f)    

data = {'Start':newsizestart,'End':newsizeend}
dnewlarge = pd.DataFrame(data)


fichier1 = open(f'Infos.txt', "a")
fichier1.write(f'Taille de la zone: {Longueur} pb\n')
fichier1.write(f'Densité moyenne de snps sur le génome total {den}\n')
fichier1.write(f'Cela correspont en moyenne à 1 SNP tous les {dens} nucléotides \n')
fichier1.write(f'Moyenne attendue : {moyatt} SNPs sur chaque zone\n')
fichier1.write(f'{dnewlarge} \n')
fichier1.close()


#record = GraphicRecord(sequence_length=seql, features=Feature)
#ax, _ = record.plot(figure_width=25)
#ax.figure.savefig("SNPs.png")
print("It's over!")




        





