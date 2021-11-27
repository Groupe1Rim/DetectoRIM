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
#Seqinteret=input('Entrez le nom de votre séquence:  ')
#seqRecord3 = SeqIO.read('Seqinterest', 'fasta')
#seql=len(seqRecord3)
#tableau=input('Entrez le nom de votre document contenant les SNPs:  ')
#df=pd.read_csv('tableau',sep="\t")



#####################Version devellopeur###############################
seqRecord3 = SeqIO.read('Sequenceconsensusmix2versionest.fasta', 'fasta')
print(seqRecord3)
seql=len(seqRecord3.seq)
print(seql)
#On importe la séquence d'interet et on calcule sa taille
df=pd.read_csv('SNPSMIX4S_11V2',sep="\t")
#On récupère le Document contenant les snps
ahhh=input('Qui est la premiere séquence?: ')
beeeh=input('Qui est la seconde séquence?: ')


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
#dfinal.to_excel(f'SNPsfinal.xlsx')


dfr=pd.read_csv('SNPsfinal.csv',sep=",")



snpsvrai=[]
position2=[]
gapinf=[]
gapsup=[]
listeinfo=[]
for k in dfr.index:
    if dfr.Gapinf[k]<5 or dfr.Gapsup[k]<5:
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
    if dd.Gapinf[i]<5 and dd.Gapsup[i]<5:
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

#dd4.to_csv(f'Fin2.csv')
#dd4.to_excel(f'Fin2.xlsx')


tailleA=0
tailleB=0




for k in dd4.index:
    if dd4.Origine[k]==['A']:
        tailleA=tailleA+dd4.Taille[k]
    if dd4.Origine[k]==['B']:
        tailleB=tailleB+dd4.Taille[k]

identityA=((tailleA/seql)*100)
print(f'{ahhh}: {identityA}%')
identityB=((tailleB/seql)*100)    
print(f'{beeeh}: {identityB}%')

fichier1 = open(f'R4_S211v2infos.txt', "a")
fichier1.write(f'Pourcentage d ADN provenant de {ahhh}: {identityA}%\n')
fichier1.write(f'Pourcentage d ADN provenant de {beeeh}: {identityB}%\n')
fichier1.close()


for j in dd4.index:
    if dd4.Origine[j]==['A']:
        
       Feature.append(GraphicFeature(start=dd4["start"][j], end=dd4["end"][j],color="#cffccc",label='A'))
       
    if dd4.Origine[j]==['B']:
        Feature.append(GraphicFeature(start=dd4["start"][j], end=dd4["end"][j],color="#ffcccc",label='B'))






record = GraphicRecord(sequence_length=seql, features=Feature)
print("Je vais commencer à tracer le graphique, cela risque de prendre quelques minutes :)" )
ax, _ = record.plot(figure_width=250)
ax.figure.savefig("4_S211v2.png")
print("C'est fini!")





finduprogamme = time.process_time()
temps=finduprogamme-débutduprogramme
print(f'Temps écoulé : {temps}s')
temps2=temps/60
print(f'Temps écoulé : {temps2}min')



