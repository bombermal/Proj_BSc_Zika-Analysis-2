# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 17:37:26 2017

@author: Ivan Alisson Cavalcante Nunes de Lima
"""
class node:
    """
        Class tha stores each string in a node, with a character, a index before alignment and a index after the alignment
    """
    amino = None
    degree = None
    aminoId = None
    
    def __init__(self, amino, degree, aminoId):
        self.amino = amino
        self.degree = degree
        self.aminoId = aminoId

    def getAmino(self):
        return self.amino
    
    def setAmino(self, val):
        self.amino = val
        
    def getDegree(self):
        return self.degree
    
    def setDegree(self, val):
        self.degree = val
    
    def getAminoId(self):
        return self.aminoId
    
    def setAminoId(self, val):
        self.aminoId = val
        
    def getAll(self):
        return self.amino, self.degree, self.aminoId
    
    def setAll(self, amino, degree, aminoId):
        self.amino = amino
        self.degree = degree
        self.aminoId = aminoId
