import SmithWaterman2 as sw
import FinalDF as fd
import WorkFunctions as wFunc
from tqdm import tqdm
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing as mp
import pandas as pd


class work(object):
    sequence  = None
    pdbs      = None
    showAlign = False
    seqNode1  = []
    seqNode2  = []
    finalDF = None
    
    def getNode(self):
        return self.seqNode2, self.seqNode1
    
    def workDivide(self, dfWork):
        """
            Divide the DataFrame in a list of rows, to use in multiprocess
        """
        listWork = []
        for i in range(len(dfWork)):
            listWork.append(dfWork.iloc[i])
        return listWork

    def prepareWork(self, pdbs, sequence, single, showAlign=False, saveAlign=False):
        """
            Converte os PDBS em dicionários, o nome da proteína é a key,
            a sequencia o valor.

            Além disso, converte a sequencia para a lista de nós usada no mapeamento
        """        
        for key in pdbs:
            pdbs[key] = wFunc.makeNodes(pdbs[key])
        
        for ii in single:
            del pdbs[ii]
        
        tempSeq = {}
        sequence = sequence[["ID", "FullSeq"]].copy()
        for ii in range(len(sequence)):
            tempSeq[sequence.loc[ii, "ID"]] = wFunc.makeNodes(sequence.loc[ii, "FullSeq"], False)
        
        self.showAlign = showAlign
        self.saveAlign = saveAlign
        self.sequence = tempSeq
        self.pdbs = pdbs
        #self.finalDf = fd.finalDF(list(pdbs.keys()), sequence, ["ID", "Protein", "Seq", "Len"])
        final, finalDf = self.makeAllWork()
        
        return final, finalDf

    def singleWork(self, df):
        seqSample = self.sequence[df]
        pdbSample = self.pdbs
        
        seq1 = seqSample
        
        columns = ["Protein", "Sample_ID", "Seq", "Cover"]
        tempDf = pd.DataFrame(columns=columns)
        finalDf = pd.DataFrame(columns=["ID", "Protein", "Seq", "Len"])
        
        for key in self.pdbs:
            seq2 = self.pdbs[key]

            obj = sw.smithWaterman()
            temp2, temp1, self.seqNode2, self.seqNode1 = obj.constructor(2, -1, -1, seq1, seq2, self.showAlign, self.saveAlign)
            tempDf.loc[len(tempDf)] = [key, df, temp1, obj.getCut()]
            
            #A treta tá aqui
            aux = "|".join(str(x.getAmino())+","+str(x.getDegree()) for x in self.seqNode1)
            #self.finalDf.setCell(df, key, aux, len(aux))
            finalDf.loc[len(finalDf)] = [df, key, aux, len(aux)]
        
        return tempDf, finalDf


    def makeAllWork(self):
        """
            Start the process pipeline

        """
        finalDf = pd.DataFrame(columns=["ID", "Protein", "Seq", "Len"])
        poolSize = mp.cpu_count()
        pool = ThreadPool(poolSize)
        results = []
        
        for ii, jj in tqdm(pool.imap_unordered(self.singleWork, self.sequence), total=len(self.sequence)):
            results.append(ii)
            finalDf = pd.concat([finalDf, jj])
        pool.close()
        pool.join()
        
        return wFunc.joinRows(results), finalDf
    
    def getPdbs(self):
        return self.pdbs
    
    def makeFasta(self, df, dataList):
        for ii in dataList:
            temp = df[df.Protein == ii]
            for jj in tqdm(temp.iterrows()):
                path = "read/Aligned/"+ii+".fa"
                f = open(path, "a")
                try:
                    f.write(">{}\n{}\n".format(jj[1].Sample_ID, jj[1].Seq ))
                finally:
                    f.close()
