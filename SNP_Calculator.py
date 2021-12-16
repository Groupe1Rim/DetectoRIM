############################################################################Importation des modules############################################################################
import pandas as pd
import statistics
import openpyxl
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import time
from Bio import SeqFeature



############################################################################################################################################################################################
###########################################################################################################################################################################################
############################################################################################################################################################################################

débutduprogramme = time.perf_counter()
#Permet de calculer le temps d'action du programme


#####################Version interactive###############################
#Seqinteret=input('Entrez le nom de votre séquence:  ')
#seqRecord3 = SeqIO.read('Seqinterest', 'fasta')
#seql=len(seqRecord3)
#tableau=input('Entrez le nom de votre document contenant les SNPs:  ')
#df=pd.read_csv('tableau',sep="\t")



#####################Version devellopeur###############################
seqRecord3 = SeqIO.read('Sequence4s11ououi.fasta', 'fasta')
seql=len(seqRecord3)


#print(seql)
#On importe la séquence d'interet et on calcule sa taille
df=pd.read_csv('SNPS',sep="\t")
#On récupère le Document contenant les snps
number=[1,2,3]
for i in number:
    indexNames = df[ df[f'sequence_{i}_PosInContg'] == 0 ].index
    df.drop(indexNames , inplace=True)

#cette partie permet de supprimer les lignes avec des valeurs à zéro (C'est normalement impossible, mais mieux vaut ne pas prendre de risque)

    

snp_pattern= df[f'SNP_pattern'].tolist()
Position=[]
A=[]
B=[]
R=[]
for label, row in df.iterrows():
    
    position=df.sequence_3_PosInContg[label]
    Position.append(position)
    snp=df.SNP_pattern[label]
    patternpA=snp[0]
    A.append(patternpA)
    patternpB=snp[1]
    B.append(patternpB)
    patternR=snp[2]
    R.append(patternR)

data = {'SNP':snp_pattern,'A':A,'B':B,'R':R,'Position':Position}
dfr = pd.DataFrame(data)
dfr['index'] = dfr.index


for i in range(len(dfr)):
    if dfr.loc[i,"A"]==dfr.loc[i,"B"]:
        dfr.drop(i,0,inplace=True)
        
dsort=dfr.sort_values(by=['Position'])
dsort['index'] = dsort.index


for i in dsort.index:
    if dsort.A[i]!= dsort.R[i]and dsort.A[i]!= dsort.B[i] and dsort.B[i]!= dsort.R[i]:
        dsort.drop(i,0,inplace=True)
simple=[]
infsnp=[]
Positions=[]


for i in dsort.index:
    Positions.append(dsort.Position[i])
    if dsort.A[i]==dsort.R[i]:
        simple.append('A')
    if dsort.B[i]==dsort.R[i]:
        simple.append('B')


        
tui=len(Positions)

tur=tui+1
Gapinf=[]
Gapsup=[]

for u in range(tui):
    
    ty=u+1
    tu=u-1
    if ty>=tui:
        ty=u
    else:
        ty=ty

    if tu<0:
        tu=u
    else:
        tu=tu
    
    tinf=Positions[u]-Positions[tu]
    Gapinf.append(tinf)
    
    tsup=Positions[ty]-Positions[u]
    Gapsup.append(tsup)
    
    
   
data = {'Origine':simple,'Position':Positions,'Gapinf':Gapinf,'Gapsup':Gapsup}
dfinal = pd.DataFrame(data)
dfinal.to_csv(f'SNPsfinal.csv')
dfinal.to_excel(f'SNPsfinal.xlsx')


dfy=pd.read_csv('SNPsfinal.csv',sep=",")
dfy['index'] = dfy.index

yu=len(dfy.index)
yut=yu-1

listesnp = dfy['Position'].tolist()
listenature=[]
listepos=[]
list_ = list(range(seql))
oip=list_[0::1000000]


yu=len(dfy.index)
yht=(yu/seql)*100
#print(f"Densité SNPs sur le génome entier: {yht}%")

listexcl=[]
lstart=[]

for r in oip:
    pou=0
    t=r+1000000
    jr=r/1000000
    jt=t/1000000
    
    lstart.append(f'{jr}-{jt}')
    
    for i in listesnp:
        if r<i<t:
          pou=pou+1
    if t == 12000000:
       tuy=(pou/(seql-11000000))*100
       listexcl.append(tuy)
    else:
        tuy=(pou/1000000)*100
        listexcl.append(tuy)


    
    #print(f"Densité SNPs entre la position{r} et {t}: {tuy}")
listexcl.append(yht)
lstart.append('densitémoyenne')
data = {'Fenêtre en MgB':lstart,'Densité en SNP(en %)':listexcl}
pio = pd.DataFrame(data)
pio.to_excel(f'SNPsdensity.xlsx')


