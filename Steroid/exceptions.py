# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:29:08 2022

@author: antoine
"""

class EmptyListError(Exception):
    def __init__(self, message = "The sequence is empty"):
        super().__init__(message)
        
class NoneDectector(Exception):
    def __init__(self, message = "Decetor = None in Potometry object"):
        super().__init__(message)
        
class ApertureOnTheEdge(Exception):
    def __init__(self, imageIdx):
        super().__init__("apperture reach the edge on frame " + str(imageIdx))