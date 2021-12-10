############################################################################Importation des modules############################################################################
import pandas as pd
import statistics
import openpyxl
from dna_features_viewer import (
    GraphicFeature,
    GraphicRecord,
    CircularGraphicRecord,
)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import time
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature



############################################################################################################################################################################################
###########################################################################################################################################################################################
############################################################################################################################################################################################

débutduprogramme = time.perf_counter()
#Permet de calculer le temps d'action du programme


#####################Version interactive###############################
Seqinteret=input('Enter the name of the recombinant sequence:  ')
seqRecord3 = SeqIO.read('Seqinterest', 'fasta')
seql=len(seqRecord3)
tableau=input('Enter the name of the document containing all the SNPs:  ')
df=pd.read_csv('tableau',sep="\t")



#####################Version devellopeur###############################
#seqRecord3 = SeqIO.read('Sequence4s11ououi.fasta', 'fasta')
#seql=len(seqRecord3)


print(seql)
#On importe la séquence d'interet et on calcule sa taille
#df=pd.read_csv('SNPS',sep="\t")
#On récupère le Document contenant les snps
ahhh=input('Which is the first parent sequence?: ')
beeeh=input('Which is the second parent sequence?: ')
stringe=input('stringence?:')
stringe=int(stringe)
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
dfinal.to_excel(f'total_SNPs.xlsx')


dfr=pd.read_csv('SNPsfinal.csv',sep=",")



snpsvrai=[]
position2=[]
gapinf=[]
gapsup=[]
listeinfo=[]
for k in dfr.index:
    if dfr.Gapinf[k]<stringe or dfr.Gapsup[k]<stringe:
        snpsvrai.append(dfr.Origine[k])
        position2.append(dfr.Position[k])
        gapinf.append(dfr.Gapinf[k])
        gapsup.append(dfr.Gapsup[k])
        
data={'Origine':snpsvrai,'Position':position2,'Gapinf':gapinf,'Gapsup':gapsup}
dd= pd.DataFrame(data)
dd['index'] = dd.index
#dd.to_csv(f'snpfin.csv')
#dd.to_excel(f'snpfin.xlsx')

dd=dd.drop(dd.index[[0]])

listeoriginechange=[]
listepositionchange=[]
listeinfochange=[]


for i in dd.index:    
    if dd.Gapinf[i]<stringe and dd.Gapsup[i]<stringe:
        listeinfo.append("continu")     
    else:
        listeinfo.append("change")
        listeinfochange.append("change")
        listepositionchange.append(dd.Position[i])
        listeoriginechange.append(dd.Origine[i])


data={'Origine':listeoriginechange,'Position':listepositionchange,'info':listeinfochange}
dd3=pd.DataFrame(data)
dd3['index'] = dd3.index
#dd3.to_csv(f'Infos.csv')
#dd3.to_excel(f'Infos.xlsx')

listestart=[]
listeend=[]
listecool=[]


listeposit= dd3[f'Position'].tolist()
for x in listeposit[0::2]:
    listestart.append(x)
for y in listeposit[1::2]:
    listeend.append(y)
for p in listestart:
    ty=dd3[dd3[f'Position']==p].Origine.tolist()
    listecool.append(ty)

lons=len(listestart)
lone=len(listeend)


if lons!=lone:
    if lons<lone:
        listestart.append(listestart[lons-1])    
    if lons>lone:
        listeend.append(listeend[lone-1])


data={'Origine':listecool,'start':listestart,'end':listeend}
dd4=pd.DataFrame(data)
dd4['index'] = dd4.index
#dd4.to_csv(f'Fin.csv')
#dd4.to_excel(f'Fin.xlsx')
Feature=[]
Feature2=[]

yu=len(dd4.index)

yut=yu-1


dd4.loc[0]=[ dd4.Origine[0], 0,0,0 ]
dd4.loc[yut]=[dd4.Origine[yut], seql,seql,yut]

yuj=yut-1
Profit=[]
for b in dd4.index:
    m=b+1
    if m>yut:
        m=b
    else:
        m=m
    q=dd4.end[m]-dd4.end[b]
    Profit.append(q)

    
    
    


dd4['Taille'] = Profit

dd4.to_csv(f'Fin2.csv')
dd4.to_excel(f'Sorted_SNPs.xlsx')


tailleA=0
tailleB=0




for k in dd4.index:
    if dd4.Origine[k]==['A']:
        tailleA=tailleA+dd4.Taille[k]
    if dd4.Origine[k]==['B']:
        tailleB=tailleB+dd4.Taille[k]

identityA=((tailleA/seql)*100)
print(f'{ahhh}: {identityA}%')
print(tailleA)
identityB=((tailleB/seql)*100)    
print(f'{beeeh}: {identityB}%')
print(tailleB)

fichier1 = open(f'Recombinant_sequence.txt', "a")
fichier1.write(f'Pourcentage d ADN provenant de {ahhh}: {identityA}%\n')
fichier1.write(f'Pourcentage d ADN provenant de {beeeh}: {identityB}%\n')
fichier1.close()

lengthA = 0
lengthB = 0
elemsA = []
elemsB = []
ro=len(dd4.index)
ro=ro-1

for j in range(ro):
    if dd4.Origine[j]==['A']:
       distance =  dd4["start"][j+1]-dd4["start"][j]
       lengthA+=distance
       elemsA.append(j)
       if dd4.Origine[j+1]==['B']:
            startpoint = dd4["start"][elemsA[0]]
            endpoint = dd4["start"][elemsA[0]]+lengthA
            Feature.append(GraphicFeature(start=startpoint, end=endpoint,color="#cffccc",label=f'{ahhh}'))
            elemsA = []
            lengthA = 0
    if dd4.Origine[j]==['B']:
        distance =  dd4["start"][j+1]-dd4["start"][j]
        lengthB+=distance
        elemsB.append(j)
        if dd4.Origine[j+1]==['A']:
           startpoint = dd4["start"][elemsB[0]]
           endpoint = dd4["start"][elemsB[0]]+lengthB
           Feature.append(GraphicFeature(start=startpoint, end=endpoint,color="#ffcccc",label=f'{beeeh}'))
           elemsB = []
           lengthB = 0

#Feature.append(GraphicFeature(start=11456842, end=11536381 ,color="#ffcccc",label=f'{beeeh}'))


record = GraphicRecord(sequence_length=seql, features=Feature)

ax, _ = record.plot(figure_width=25)
ax.figure.savefig("Recombinant_sequence.png")
print("It's over!")





finduprogamme = time.process_time()
temps=finduprogamme-débutduprogramme
print(f'Elapsed time : {temps}s')
temps2=temps/60
print(f'Elapsed time : {temps2}min')



