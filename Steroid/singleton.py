# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 13:50:04 2022

@author: antoine
"""

class SingletonMeta(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(SingletonMeta, cls).__call__(*args, **kwargs)
        return cls._instances[cls]
    
    
class Singleton(metaclass = SingletonMeta):
    def __init__(self):
        self.dataDeleted = []
        self.dataOfStarPassages = []
        
    def getDataDeleted(self):
        return self.dataDeleted
    
    def getDataStarsPassages(self):
        return self.dataOfStarPassages
    
    def appDataDeleted(self, obj):
        self.dataDeleted.append(obj)
        
    def appDataStarsPassages(self, obj):
        self.dataOfStarPassages.append(obj)
        
    def clearDataDeleted(self):
        self.dataDeleted = []
    
    def clearDataStarsPassages(self):
        self.dataOfStarPassages = []
        
    def clear(self):
        self.clearDataDeleted()
        self.clearDataStarsPassages()
        
        