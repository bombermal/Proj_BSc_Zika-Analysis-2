import matplotlib.pyplot as plt
import WorkFunctions as wFunc
import pandas as pd
import numpy as np

def plotGraph(df, title, color=True, ceiling=1000000, floor=0):
    color_dict = {'-' : 'Black', 'A' : '#29a329', 'C' : '#cca300', 'D' : '#3399ff', 'E' : '#3399ff', 'F' : '#29a329', 'G' : '#cca300', 'H' : '#004d99', 'I' : '#29a329', 'K' : '#004d99', 'L' : '#29a329', 'M' : '#29a329', 'N' : '#a64dff', 'P' : '#cca300', 'Q' : '#a64dff', 'R' : '#004d99', 'S' : '#a64dff', 'T' : '#a64dff', 'U': '#cca300', 'V' : '#29a329', 'W' : '#29a329', 'X' : '#FF0000', 'Y' : '#29a329'}
    df = df[(df >= floor) & (df <= ceiling)]
    if color:
        ax = df.T.plot(kind="bar",  stacked=True, figsize=(200,10), title=title, color=[color_dict.get(x, '#333333') for x in df.index])
    else:
        ax = df.T.plot(kind="bar",  stacked=True, figsize=(200,10), title=title, colormap="tab20")
    plt.legend(loc=2, prop={'size': 15})
    plt.show()
    
    
"""
Plot graph showing the polimorphic regions.
    X - axis: base position
    Y - axis: Aminoacid count
    
    dataList: List of protein titles
    poliLIstAmino: List of aminoacids counts for each sample grouped by proteins
    mode:
        if "all", plot all graphs for all proteins
        else: "ProteinTitle" will define wich protein will be plotted
            
    color: 
        True: plot graph separating the aminoacids by chemic properties
        False: each color will represent one tipe of aminoacid
            
    size: is a int number that determine how many aminoacids will be needed in a postion
        for this this position be rotulated as a polymorphic base
"""
def plotPolimorph(dataList, poliListAmino, mode, color=True, size=4):
    if mode == "all":
        for ii, jj in zip(dataList, poliListAmino):
            print("Overview PDB {}:".format(ii))
            plotGraph(countPoli(jj), ii, color)
            print("Polimorphism PDB {}:".format(ii))
            plotGraph(selectColumnsForPlot(countPoli(jj),ii, size), ii, color)
    else:
        idx = dataList.index(mode)
        print("Overview PDB {}:".format(mode))
        plotGraph(countPoli(poliListAmino[idx]), mode, color)
        print("Polimorphism PDB {}:".format(mode))
        plotGraph(selectColumnsForPlot(countPoli(poliListAmino[idx]),mode, size), mode, color)

def countPoli(df):
    letters = [ '-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V' , 'W', 'X', 'Y']
    columns = df.iloc[:,1:].columns
    temp = pd.DataFrame(index=letters, columns=columns)
    for ii in columns:
        tt = df[ii].value_counts()
        idx = tt.index[0]
        temp.loc[idx, ii] = tt[0]
        
    return temp
        
def plotCoverRange(df):
    df.Cover.plot(kind="hist", title="Overview identities", colormap="tab20")
    plt.show()        

def flatDF(df):
    temp = []
    for col in df:
        for jj in df[col].tolist():
            temp.append(jj)
    return pd.Series(temp).dropna()

def plotDegrees(dataList, nsFiles, poliListDegree, poliListAmino, mode, polimorphCut=4, showDF=False):
    if mode == "all":
        for idx, ii in enumerate(dataList):
            proteinID = (idx, ii)
            plotSingleDegree(nsFiles, poliListDegree[idx], poliListAmino, proteinID, polimorphCut, showDF)
    else:
        idx = dataList.index(mode)
        proteinID = (idx, mode)
        plotSingleDegree(nsFiles, poliListDegree[idx], proteinID, polimorphCut, showDF)
        
def plotSingleDegree(nsFiles, poliListDegree, poliListAmino, proteinID, polimorphCut, showDF):
    idx = proteinID[0]
    proteinID = proteinID[1]
    poliListDegreeSelect = flatDF(poliListDegree[wFunc.findPoli(poliListAmino[idx], polimorphCut).columns].copy())

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
    ax1.hist(nsFiles[proteinID].Degree.tolist())
    ax1.set_title("Degrees distribution on PDB "+proteinID)  
    ax1.set_xlabel("Degree")
    ax1.set_ylabel("Frequency")
    ax1.grid(axis="x", linestyle='dotted')
    #ax1.set_xlim(0,15)
    ax2.hist(poliListDegreeSelect)
    ax2.set_title("Degrees distribution on polymorphic bases for "+proteinID)  
    ax2.set_xlabel("Degree")
    ax2.set_ylabel("Frequency")
    ax2.grid(axis="x", linestyle='dotted')
    #plt.xticks(np.arange(0, 15, step=1))
    plt.tight_layout()
    plt.savefig("imgs/Degree_"+proteinID+".png", dpi=600)
    plt.show()
    if showDF:
        print(poliListDegree.value_counts())
        
def plotBtw(file, name):
    name = str.upper(name[:4])
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
    ax1.hist(file.clusteringCoef.tolist())
    ax1.set_title("Clustering Coefficient. on PDB "+name)  
    ax1.set_xlabel("Coefficient")
    ax1.set_ylabel("Frequency")
    ax1.grid(axis="x", linestyle='dotted')
    ax2.hist(file.betweennessWeighted.tolist())
    ax2.set_title("Betweenness Weighted on PDB "+name)  
    ax2.set_xlabel("Betweenness")
    ax2.set_ylabel("Frequency")
    ax2.grid(axis="x", linestyle='dotted')
    plt.tight_layout()
    plt.savefig("imgs/Btween_"+proteinID+".png", dpi=600)
    plt.show()
    
def allPDBs(nsFiles):
    temp = []
    aux =[]
    for keys in nsFiles:
        temp.append(nsFiles[keys].Degree.tolist())

    for ii in temp:
        for jj in ii:
            aux.append(jj)

    fig, ax1 = plt.subplots(figsize=(10,5))
    ax1.hist(aux, 50)
    ax1.set_title("Degrees distribution on all PDBs")  
    ax1.set_xlabel("Degree")
    ax1.set_ylabel("Frequency")
    ax1.grid(axis="x", linestyle='dotted')
    plt.xticks(np.arange(0, 50, step=1))
    plt.savefig('imgs/PDBs_All.png', dpi=600)
    plt.tight_layout()
    plt.show()
    
def selectColumnsForPlot(jj, ii, size):
    poli = wFunc.findPoli(jj, size=size)
    temp = jj.copy()
    temp[temp >0] = 0
    for ii in poli.columns:
        temp[ii] = poli[ii].values
    
    return temp