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




débutduprogramme = time.perf_counter()
df=pd.read_csv('bacbkone4s_11modif.txt',sep="\t")

df['index'] = df.index



seqRecord = SeqIO.read('4_S211_consensus_S1A1-3_.fasta', 'fasta')
seqRecord2 = SeqIO.read('4_S211_consensus_S1D4-23_.fasta', 'fasta')
print(seqRecord)
seql1=len(seqRecord.seq)
print(seql1)

seql2=len(seqRecord2.seq)
print(seqRecord2)
print(seql2)


sequenceconsensus=seqRecord.seq
sequenceconsensus=list(seqRecord.seq.strip())

print(len(sequenceconsensus))







dj=pd.read_csv('bacbkone4s_11modif.txt',sep="\t")
dj['index'] = dj.index


for y in dj.index:
    A2=(dj.seq0_leftend[y])+(dj.seq0_rightend[y])
    u=y-1
    if A2==0:

        s = dj[dj[f'index']==y].seq1_leftend.tolist()
        start1=statistics.median(s)
        e = dj[dj[f'index']==y].seq1_rightend.tolist()
        end1=statistics.median(e)
        recevdebut=dj.seq0_rightend[u]
        
        insertion=seqRecord2.seq[start1:end1]
        sequenceconsensus.insert(recevdebut, insertion)
        
        


print(len(sequenceconsensus))
seq2 = ''.join(str(elem) for elem in sequenceconsensus)
fichier1 = open(f'Sequenceconsensusmix2v2.fasta', "a")
fichier1.write(f'>Sequenceconsensusmix2v2 Streptomyces\n')
fichier1.write (f'{seq2}')
fichier1.close()



seqRecordp = SeqIO.read('Sequenceconsensusmix2v2.fasta', 'fasta')
print(seqRecordp)
seql=len(seqRecordp.seq)
print(seql)








finduprogamme = time.process_time()
temps=finduprogamme-débutduprogramme
print(f'Temps écoulé : {temps}s')
temps2=temps/60
print(f'Temps écoulé : {temps2}min')




