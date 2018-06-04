import pandas as pd

class finalDF:
    finalDF = None
    
    def __init__(self, dataList, df, columns):
        self.finalDF = pd.DataFrame(columns=columns)
        for ii in dataList:
            for jj in range(len(df)):
                self.finalDF.loc[len(self.finalDF)] = [ df.loc[jj, "ID"], ii, 0, 0]
                
    def getDf(self):
        return self.finalDF
    
    def setCell(self, sampleId, proteinID, val, size):
        idx = self.finalDF[(self.finalDF.ID == sampleId) & (self.finalDF.Protein == proteinID)].index[0]
        self.finalDF.loc[idx, "Seq"] = val
        self.finalDF.loc[idx, "Len"] = size
        
    def concatCell(self, df):
        self.finalDF = pd.concat([self.finalDF, df])