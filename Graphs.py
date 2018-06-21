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
    
def plotPolimorph(dataList, poliListAmino, mode, color=True, size=4):
    if mode == "all":
        for ii, jj in zip(dataList, poliListAmino):
            print("Overview PDB {}:".format(ii))
            plotGraph(jj, ii, color)
            print("Polimorphism PDB {}:".format(ii))
            plotGraph(selectColumnsForPlot(jj,ii, size), ii, color)
    else:
        idx = dataList.index(mode)
        print("Overview PDB {}:".format(mode))
        plotGraph(poliListAmino[idx], mode, color)
        print("Polimorphism PDB {}:".format(mode))
        plotGraph(selectColumnsForPlot(poliListAmino[idx],mode, size), mode, color)
    
def plotCoverRange(df):
    df.Cover.plot(kind="hist", title="Overview identities", colormap="tab20")
    plt.show()        

def flatDF(df):
    temp = []
    for col in df:
        for jj in df[col].tolist():
            temp.append(jj)
    return pd.Series(temp).dropna()

def plotDegrees(dataList, poliListDegree, poliListAmino, mode, polimorphCut=4, showDF=False):
    if mode == "all":
        for idx, ii in enumerate(dataList):
            proteinID = (idx, ii)
            plotSingleDegree(poliListDegree[idx], poliListAmino, proteinID, polimorphCut, showDF)
    else:
        idx = dataList.index(mode)
        proteinID = (idx, mode)
        plotSingleDegree(poliListDegree[idx], proteinID, polimorphCut, showDF)
        
def plotSingleDegree(poliListDegree, poliListAmino, proteinID, polimorphCut, showDF):
    idx = proteinID[0]
    proteinID = proteinID[1]
    poliListDegreeSelect = flatDF(poliListDegree[wFunc.findPoli(poliListAmino[idx], polimorphCut).columns].copy())
    poliListDegree = flatDF(poliListDegree)

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10,5), sharex=True, sharey=True)
    ax1.set_title("Ovewall degrees distribution in "+proteinID)
    ax1.hist(poliListDegree)
    ax1.set_xlabel("Degree")
    ax1.set_ylabel("Frequency")
    ax1.set_xlim(0,14)
    ax1.grid(axis="x", linestyle='dotted')
    ax2.hist(poliListDegreeSelect)
    ax2.set_title("Degrees distribution on polimorphic bases for "+proteinID)  
    ax2.set_xlabel("Degree")
    ax2.set_ylabel("Frequency")
    ax2.grid(axis="x", linestyle='dotted')
    plt.xticks(np.arange(0, 14, step=1))
    plt.show()
    if showDF:
        print(poliListDegree.value_counts())
        
def selectColumnsForPlot(jj, ii, size):
    poli = wFunc.findPoli(jj, size=size)
    temp = jj.copy()
    temp[temp >0] = 0
    for ii in poli.columns:
        temp[ii] = poli[ii].values
    
    return temp