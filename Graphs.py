import matplotlib.pyplot as plt
import WorkFunctions as wFunc

def plotGraph(df, title, ceiling=1000000, floor=0):
    df = df[(df >= floor) & (df <= ceiling)]
    ax = df.T.plot(kind="bar",  stacked=True, figsize=(200,10), title=title, colormap="tab20")
    plt.legend(loc=2, prop={'size': 15})
    plt.show()
    
def plotCoverRange(df):
    df.Cover.plot(kind="hist", title="Overview identities", colormap="tab20")
    plt.show()
    
def plotPolimorph(dataList, poliListAmino, size):
    for ii, jj in zip(dataList, poliListAmino):
        print("Overview PDB {}:".format(ii))
        plotGraph(jj, ii)
        print("Polimorphism PDB {}:".format(ii))
        plotGraph(selectColumnsForPlot(jj,ii, size), ii)
        
        
def plotDegrees(dataList, poliListDegree, poliListAmino):
    for aa, (ii, jj) in enumerate(zip(dataList, poliListDegree)):
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
        title = "Degrees distribution in "+ii
        jj.T.plot(kind="hist", title=title, colormap="tab20", alpha=0.5, legend=False, ax=axes[0])
        aux =jj[wFunc.findPoli(poliListAmino[aa], 4).columns].copy()
        aux.T.plot(kind="hist", title=title, colormap="tab20", alpha=0.5, legend=False, ax=axes[1])
        plt.show()
        
def selectColumnsForPlot(jj, ii, size):
    poli = wFunc.findPoli(jj, size=size)
    temp = jj.copy()
    temp[temp >0] = 0
    for ii in poli.columns:
        temp[ii] = poli[ii].values
    
    return temp