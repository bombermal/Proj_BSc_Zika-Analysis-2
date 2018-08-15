import pandas as pd
import sqlFunctions as sqlF
import workFunctions as wFunc
import readAndSaveFiles as rSF

def populateAll():
    protein = "read/Protein - C_NS5"
    fastaProSequences = rSF.fileRead(protein, ".fa")
    columns=["ID", "Host", "Region", "Date", "FullSeq"]
    #Make DataFrames
    fastaProDF = wFunc.makeDF(columns, fastaProSequences)
    nsFiles = rSF.readPDBs()
    #insert PDB Names
    pdbNames = ["5JMT","5TMH", "5GS6", "5IY3", "5K6K" ]
    populatePdbNameList(pdbNames)

    for key in nsFiles:
        sqlF.insertDataframe(nsFiles[key], "pdbs_nodes")
        
    nsFiles = rSF.readPDBs(True)
    for key in nsFiles:
        sqlF.insertDataframe(nsFiles[key], "pdbs_edges")

    #insert samples on Db
    sqlF.insertDataframe(fastaProDF, "samples")
    
def populatePdbNameList(names):
    columns = sqlF.returnColumnsList("pdbs_names_list")
    for ii in names:
        sqlF.insertData("pdbs_names_list", columns, [(ii, '')])