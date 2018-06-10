import Node as nd
import pandas as pd
import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool
from Bio.SeqUtils import seq3, seq1 as seqX
from Bio import SeqIO
from tqdm import tqdm

def makeNodes(seq, pdb=True):
    """
        Make list of nodes from a string
    """
    temp = []
    for ii in range(len(seq)):
        if pdb:
            temp.append(nd.node(seqX(seq.loc[ii, "Residue"]), seq.loc[ii, "Degree"]))
        else:
            temp.append(nd.node(seq[ii], -10))

            
    return temp        
 
    
def workDivide(dfWork):
    """
        Divide the DataFrame ina a list of rows, to use in multiprocess
    """
    listWork = []
    for i in range(len(dfWork)):
        listWork.append(dfWork.iloc[i])
    return listWork

 
def joinRows(listFinal):
    """
        Join all rows in a single dataframe
    """
    temp = listFinal[0]
    for ii in range(1,len(listFinal)):
        temp = pd.concat([temp, listFinal[ii]])
    temp.reset_index(inplace=True)
    temp.drop("index", axis=1, inplace=True)
    return temp

def matchDegrees(seq1, seq2):
    for ii, jj in zip(seq1, seq2):
        if ii.getAmino() == jj.getAmino():
            jj.setDegree(ii.getDegree())
        elif ii.getAmino() == "-":
            jj.setDegree(-2)
        else:
            jj.setDegree(ii.getDegree())
            
    return seq1, seq2


def fileRead(path, ext):
    """
       Read files and return a list with all values
    """
    path += ext
    temp = []
    with open(path) as fasta_file:
        for seqRecord in SeqIO.parse(fasta_file, 'fasta'):
            temp.append(seqRecord)
    return temp

def readPDBs():
    """
        Read my PDBs files
    """
    dataList = ["5gs6", "5IY3", "5k6k"]
    nsFiles = {}
    nsFiles["5JMT"] = pd.read_csv("read/NS3/5JMT_nodes.txt", sep="\t", low_memory=False)
    nsFiles["5TMH"] = pd.read_csv("read/NS5/5TMH_nodes.txt", sep="\t", low_memory=False)
    for ii in dataList:
        nsFiles[ii.upper()] = pd.read_csv("read/NS1/"+ii+"_nodes.txt", sep="\t", low_memory=False)
    
    return nsFiles

def makeDF(columns, sourceList):
    """
        Make a dataFrame
        
        Parameters
        ----------
        columns : list
            List of columns names
        souceList: list
            Data to populate DF
        
        Return
        ------
        pandas.DataFrame
    """
    baseDf = pd.DataFrame(columns=columns)
    
    for itr in sourceList:
        temp = itr.description.split("|")
        baseDf.loc[len(baseDf)] = [temp[0].strip(), temp[1], temp[2], temp[3], str(itr.seq)]
    fixDateTime(baseDf)
    return baseDf

def fixDateTime(df):
    """
        Fix Column DateTime in DataFrame
    """
    for ii, atual in enumerate(df.Date):
        df.Date.loc[ii] = pd.to_datetime(atual, format='%Y/%m/%d', errors='ignore')

def fillPosition(df):
    temp = pd.DataFrame()
    for jj, kk in enumerate(df):
        temp.loc[0, str(jj)+"POS"] = kk
        
    return temp

def makeCountDF(df, proteinID, amino=True):
    selectDF = df[df.Protein == proteinID].Seq.tolist()
    selectDF = [splitCell(x, amino) for x in selectDF]

    temp = pd.DataFrame()
    poolSize = mp.cpu_count()
    pool = ThreadPool(poolSize)
    results = []
    
    for ii in tqdm(pool.imap_unordered(fillPosition, selectDF ), total=len(selectDF)):
        temp = pd.concat([temp, ii], axis=0, sort=False)
    pool.close()
    pool.join()
    
    temp = temp.reset_index().drop(["index"], axis=1)

    if amino:
        temp = temp.apply(pd.Series.value_counts)
    else:
        temp = temp.drop["-"]
        
    return temp
    

def splitCell(seq, amino):
    seq = seq.split("|")
    if amino:
        seq = [x.split(",")[0] for x in seq]
    else:
        seq = [x.split(",")[1] for x in seq]
    return seq

def countAndSavePoli(dataList, df, amino=True):
    for ii in dataList:
        temp = makeCountDF(df, ii, amino)
        if amino:
            path = "read/PolimorfDF/dfAmino"+ii+".csv"
        else:
            path = "read/PolimorfDF/dfDegree"+ii+".csv"
        temp.to_csv(path, sep="\t", index=True)
        print("PDB: "+ii+" saved!")
        
def readSavedPoli(dataList, amino=True):
    temp = []
    for ii in dataList:
        print("Reading PDB: "+ii)
        if amino:
            path = "read/PolimorfDF/dfAmino"+ii+".csv"
        else:
            path = "read/PolimorfDF/dfDegree"+ii+".csv"
        temp.append(pd.read_csv(path, sep="\t", index_col=0))
    return temp

def findPoli(df, size=4):
    temp = []
    names = df.columns.tolist()
    for ii in names:
        aux = len(df[ii].value_counts())
        if aux >= size:
            temp.append(ii)
    return df[temp].copy()

def cleanData(df, dataList):
    dataCut = {'5GS6' : 1452, '5IY3': 767, '5K6K': 1454, '5JMT': 1835, '5TMH': 3655}
    cleanedDF = pd.DataFrame()
    strangeDF = pd.DataFrame()
    for ii in dataList:
        temp = df[(df.Protein == ii) & (df.Len <= dataCut[ii])]
        temp2 = df[(df.Protein == ii) & (df.Len > dataCut[ii])]
        cleanedDF = pd.concat([cleanedDF, temp])
        strangeDF = pd.concat([strangeDF, temp2])
        
    return cleanedDF, strangeDF
    