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
    
    def __init__(self, amino, degree):
        self.amino = amino
        self.degree = degree

    def getAmino(self):
        return self.amino
    
    def setAmino(self, val):
        self.amino = val
        
    def getDegree(self):
        return self.degree
    
    def setDegree(self, val):
        self.degree = val
        
    def getAll(self):
        return self.amino, self.degree
    
    def setAll(self, amino, degree):
        self.amino = amino
        self.degree = degree
