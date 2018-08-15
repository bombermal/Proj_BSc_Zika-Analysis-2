import pandas as pd
from Bio import SeqIO

#insert PDBS - Nodes
def readPDBs(condition = False):
    """
        Read my PDBs files
    """
    src = "_nodes.txt"
    if condition:
        src = "_edges.txt"
        
    dataList = ["5gs6", "5IY3", "5k6k"]
    nsFiles = {}
    nsFiles["5JMT"] = pd.read_csv("read/NS3/5JMT"+src, sep="\t", low_memory=False)    
    nsFiles["5TMH"] = pd.read_csv("read/NS5/5TMH"+src, sep="\t", low_memory=False)
    for ii in dataList:
        nsFiles[ii.upper()] = pd.read_csv("read/NS1/"+ii+src, sep="\t", low_memory=False)
    
    for key in nsFiles:
        addColumnsAndFillNa(nsFiles[key], key)
    return nsFiles

def addColumnsAndFillNa(df, name):
    """
        Add two columns to dataframes for the Database use
    """
    df["pdbName"] = name
    df["dbIndex"] = ''
    df.fillna('', inplace = True)
    
    
def fileRead(path, ext):
    """
       Read samples files and return a list with all values
    """
    path += ext
    temp = []
    with open(path) as fasta_file:
        for seqRecord in SeqIO.parse(fasta_file, 'fasta'):
            temp.append(seqRecord)
    return temp 