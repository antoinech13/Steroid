# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:29:08 2022

@author: antoine
"""
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry
import cv2
import numpy as np
import os
from astropy.io import fits
from astropy.time import Time
from astropy.table import QTable
from astropy.stats import sigma_clipped_stats

from scipy.optimize import curve_fit

from exceptions import *
from utils import *

class Star:
    def __init__(self, s, eps = 5):
        self.s = s
        self.eps = eps
        
    def __add__(self, other):
        if isinstance(other, Star):
            return Star(self.s + other.s)
        
        return Star(self.s + other)
    
    def __sub__(self, other):
        if isinstance(other, Star):
            return Star(self.s - other.s)
        return Star(self.s - other)
    
    def __truediv__(self, other):
        if isinstance(other, Star):
            return Star(self.s / other.s)
        return Star(self.s / other)
    
    def __mul__(self, other):
        if isinstance(other, Star):
            return Star(self.s * other.s)
        return Star(self.s * other)
    
    def __eq__(self, other):
        if isinstance(other, Star):
            return (np.abs(self.s - other.s) <= self.eps).all()
        return (np.abs(self.s - other) <= self.eps).all()
    
    def __str__(self):
        return str(self.s[0]) + " " + str(self.s[1])
    
    def __len__(self):
        return self.s.shape[0]
    
    def __pow__(self, p):
        return self.s**p
    
    def toNumpyt(self):
        return np.asarray(self.s)
    
    def distance(self, other):
        return np.sum((np.asarray(self.s) - np.asarray(other.s))**2)**0.5
    
class Triangle:
    def __init__(self, s1, s2, s3, eps = 2):
        self.s1 = s1
        self.s2 = s2
        self.s3 = s3
       
        self.eps = eps
        
    def __add__(self, other):
        return Triangle(self.s1 + other.s1, self.s2 + other.s2, self.s3 + other.s3)
    
    def __sub__(self, other):
        return Triangle(self.s1 - other.s1, self.s2 - other.s2, self.s3 - other.s3)
    
    def __truediv__(self, other):
        return Triangle(self.s1 / other.s1, self.s2 / other.s2, self.s3 / other.s3)
        
    def __mul__(self, other):
        return Triangle(self.s1 * other.s1, self.s2 * other.s2, self.s3 * other.s3)
    
    def __eq__(self, other):
        return (np.abs(self.d1() - other.d1()) <= self.eps) and (np.abs(self.d2() - other.d2()) <= self.eps) and (np.abs(self.d3() - other.d3()) <= self.eps)
    
    def __str__(self):
        return "star 1: " + str(self.s1[0]) + " " + str(self.s1[1]) + " star 2: " + str(self.s2[0]) + " " + str(self.s2[1]) + " star 3: " + str(self.s3[0]) + " " + str(self.s3[1]) + "\nD1: " + str(self.d1()) + " D2: " + str(self.d2()) + " D3: " + str(self.d3()) 
    
    def d1(self):
        return np.sum((self.s2 - self.s1)**2)**0.5
        
    def d2(self):
        return np.sum((self.s3 - self.s1)**2)**0.5
    
    def d3(self):
        return np.sum((self.s3 - self.s2)**2)**0.5
            
    def buildVect(self):
        v1 = self.s2 - self.s1
        v2 = self.s3 - self.s1
        v3 = self.s3 - self.s2
        return v1, v2, v3
        
    def getRotationAngle(self, other):
        'modification ligne sign = np.sign(np.mean([sign1, sign2, sign3]))'
        
        v1, v2, v3 = self.buildVect()
        v1b, v2b, v3b = other.buildVect()
        
        scal1 = np.sum(v1*v1b)
        scal2 = np.sum(v2*v2b)
        scal3 = np.sum(v3*v3b)
        
        sign1 = np.sign(np.cross(v1,v1b))
        sign2 = np.sign(np.cross(v2,v2b))
        sign3 = np.sign(np.cross(v3,v3b))
        
        # print("sign,", sign1, sign2, sign3)
        
        c1 = scal1 / (self.d1() * other.d1()) 
        c2 = scal2 / (self.d2() * other.d2())
        c3 = scal3 / (self.d3() * other.d3())
       
        theta1 = np.arccos(c1)
        theta2 = np.arccos(c2)
        theta3 = np.arccos(c3)
        
        # print("sign,", sign1, sign2, sign3,c1,c2,c3, theta1,theta2,theta3)
        
        if c1 > 1:
            theta1 = 0
        elif c1 < -1:
            theta1 = np.pi
        if c2 > 1:
            theta2 = 0
        elif c2 < -1:
            theta2 = np.pi
        if c3 > 1:
            theta3 = 0
        elif c3 < -1:
            theta3 = np.pi
        
        
        #----------------case where two vectors are at PI --------
        if sign1 == 0:
            sign1 = 1
        if sign2 == 0:
            sign2 = 1
        if sign3 == 0:
            sign3 = 1
        #---------------------------------------------------------
        
        sign = np.sign(np.mean([sign1, sign2, sign3]))
        
        angle1 = -sign*theta1
        angle2 = -sign*theta2
        angle3 = -sign*theta3
     
        # print('ang', (angle1 + angle3 + angle2)/3)
        return (angle1 + angle3 + angle2)/3
    
    def computeDistance(self, other):
        return np.mean([self.s1 - other.s1, self.s2 - other.s2, self.s3 - other.s3], axis=0)
    
    def correctRot(self, angle, center):

        self.s1 = self.s1 - center 
        self.s2 = self.s2 - center 
        self.s3 = self.s3 - center 
     
        rot = np.asarray([[np.cos(angle), -np.sin(angle)],[np.sin(angle), np.cos(angle)]])
      
        self.s1 =  np.dot(self.s1, rot) + center
        self.s2 =  np.dot(self.s2, rot) + center
        self.s3 =  np.dot(self.s3, rot) + center

class Patern:
    def __init__(self, t1, t2, t3, t4, t5):
        self.t1 = t1
        self.t2 = t2
        self.t3 = t3
        self.t4 = t4
        self.t5 = t5
        
    def __add__(self, other):
        return Patern(self.t1 + other.t1, self.t2 + other.t2, self.t3 + other.t3, self.t4 + other.t4)
    
    def __sub__(self, other):
        return Patern(self.t1 - other.t1, self.t2 - other.t2, self.t3 - other.t3, self.t4 - other.t4)
    
    def __truediv__(self, other):
        return Patern(self.t1 / other.t1, self.t2 / other.t2, self.t3 / other.t3, self.t4 / other.t4)
        
    def __mul__(self, other):
        return Patern(self.t1 * other.t1, self.t2 * other.t2, self.t3 * other.t3, self.t4 * other.t4)
    
    def __eq__(self, other):
        return self.t1 == other.t1 and self.t2 == other.t2 and self.t3 == other.t3 and self.t4 == other.t4 and self.t5 == other.t5
    
    def __str__(self):
        return 'Triangle 1 \n ' + str(self.t1) + ' \nTriange 2 \n ' + str(self.t2) + ' \nTriangle 3 \n ' + str(self.t3) + ' \nTriangle 4 \n ' + str(self.t4) + '\nTriangle 5 \n ' + str(self.t5) + '\n'
    
    def computeDistance(self, other):
        return np.mean([self.t1.computeDistance(other.t1), self.t2.computeDistance(other.t2), self.t3.computeDistance(other.t3), self.t4.computeDistance(other.t4)], axis = 0)
    
    def computeAngle(self, other):
        
        
        a1 = self.t1.getRotationAngle(other.t1)
        a2 = self.t2.getRotationAngle(other.t2)
        a3 = self.t3.getRotationAngle(other.t3)
        a4 = self.t4.getRotationAngle(other.t4)
        a5 = self.t5.getRotationAngle(other.t5)
        # print("patern", a1,a2,a3,a4,a5,  np.median([a1,a2,a3,a4,a5]))
        return np.median([a1,a2,a3,a4,a5])

    def correctRot(self, angle, center):
        self.t1.correctRot(angle, center)
        self.t2.correctRot(angle, center)
        self.t3.correctRot(angle, center)
        self.t4.correctRot(angle, center)
        self.t5.correctRot(angle, center)

class SeqManager:
    def __init__(self, seq):
        
        if seq != None and len(seq) == 0:
            raise EmptyListError("The sequence list is empty")
        
        self.seq = seq

        
        self.darkM = 0
        self.biasM = 0
        self.flatM = 1
        
        self.darkExp = None
        self.exposurKey = None

    def __len__(self):
        return len(self.seq)
    
    def __str__(self):
        return self.seq

    def getPath(self, idx):
        return self.seq[idx]
    
    def getFileName(self, idx):
        return os.path.basename(self.getPath(idx))

    def getImg(self, i = 0):
        return Fit(self.seq[i], self.darkM, self.flatM, self.biasM, self.darkExp, self.exposurKey)

    def getHDU(self, i = 0, HDU = 0):
        return self.getImg(i).getHDU(HDU)
    
    def getInfo(self, i = 0):
        return self.getImg(i).getInfo()

    def getHeader(self, i = 0, HDU = 0):
        return self.getImg(i).getHDU(HDU).header
    
    def getExpo(self, i, key, HDU = 0):
        return self.getImg(i).getExposure(key)
    
    def getData(self, i = 0, idx_HDU = 0):
        return np.asarray(self.getImg(i).getHDU(idx_HDU).data)
    
    def getCenter(self, idx_img = 0, idx_HDU = 0):
        return self.getImg(idx_img).getCenter(idx_HDU)
    
    def getImgShape(self, idx = 0, idx_HDU = 0):
        return self.getData(idx, idx_HDU).shape
    
    def getTime(self, key, forma, idx = 0, HDU = 0):
        return self.getImg(idx).getTime(key, forma, HDU)
    
    def pop(self, idx = -1):
        '''
        Parameters
        ----------
        idx : int, optional
            DESCRIPTION. delete element of the sequence.The default is -1.

        Returns
        -------
        None.

        '''
    
        self.seq.pop(idx)
        
    def setDark(self, dark, darkExp = None, exposurKey = None):
        self.darkM = dark
        self.darkExp = darkExp
        self.exposurKey = exposurKey
    
    def setBias(self, bias):
        self.biasM = bias
    
    def setFlat(self, flat):
        self.flatM = flat
        
    def setBiasDarkFlat(self, bias, dark, flat, darkExp = None, exposurKey = None):
        self.setBias(bias)
        self.setDark(dark, darkExp, exposurKey)
        self.setFlat(flat)
    

    def getAvg(self, default = 0, keyExp = None, biasM = 0):
        
        if self.seq == None:
            return default
        
        img = self.getImg()
        data = np.zeros(img.getData().shape)
        
        if keyExp != None:
            expo = float(img.getHeader()[keyExp])
        
        scaleFactor = 1
        for i in range(len(self.seq)):
            img = self.getImg(i)
            
            if keyExp != None:
                exp = float(img.getHeader()[keyExp])
                scaleFactor = expo / exp
                
            data += (img.getData() - biasM) * scaleFactor
        return (data / len(self.seq)) 
    
    def histogram(self, idx = 0, idx_HDU = 0):
        """
        Parameters
        ----------
        idx : int, optional
            DESCRIPTION. idx of image in the sequence. The default is 0.
        idx_HDU : int, optional
            DESCRIPTION. Index of the hdu in fits file.The default is 0.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.getImg(idx).histogram(idx_HDU)
        
        
        
        

class Appertures:
    def __init__(self, positions, idxOfStars = None, r = 3, ri = 6, re = 8):
        self.positions = positions
        self.idxOfStars = idxOfStars
        self.aperture = CircularAperture(positions, r)
        self.annulus_aperture = CircularAnnulus(positions, ri, re)
        self.annulus_mask = self.annulus_aperture.to_mask(method='center')
        
        self.apmask = self.aperture.to_mask('center')
        
    def __str__(self):
        return "appertur pos: " + self.positions + " idx: " + self.idxOfStars
    
    def photom(self, img, key, forma, center = False, exposure = None):
        '''

        Parameters
        ----------
        img : Fit
            Fit object (see data_structurs/Fit).
        key : string
            Keyword of the time in the header of fits file.
        forma : string
            Format of the time in the header. refer to astropy.Time formats (jd = julian days, isot = YYYY-MM-DDTHH:MM:SS, etc...) .

        Returns
        -------
        phot : astropy.table.table.QTable
            table with :
                
            id   xcenter   ycenter    aperture_sum  Time     annulus_median   aper_bkg  aper_sum_bkgsub 
                   pix      pix                                                                                               
           int32 float64   float64      float64     str18       float64        float64    float64  .

        '''
        
        data = img.getReducedData()

        bkgMedian = []

        for i,mask in enumerate(self.annulus_mask):
            annulus_data = mask.multiply(data)
            array1D = annulus_data[mask.data > 0]
            _, median_sigclip, _ = sigma_clipped_stats(array1D)
            bkgMedian.append(median_sigclip)

        bkgMedian = np.array(bkgMedian)
        phot = aperture_photometry(data, self.aperture)
        
        if center:
            phot['Time'] = img.getTime(key, forma).jd1 + img.getTime(key, forma).jd2 + Utils.stojd(exposure/2.0)
        else:
            phot['Time'] = img.getTime(key, forma).jd1 + img.getTime(key, forma).jd2
        
        phot['annulus_median'] = bkgMedian
        phot['aper_bkg'] = bkgMedian * self.aperture.area
        phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
        # for col in phot.colnames:
        #     phot[col].info.format = '%.8g'  # for consistent table output
        return phot
    
    # def isInCenter(p):
    #     if Utils.dist()
    
    

class Fit():
    def __init__(self, path, dark = 0, flat = 1, bias = 0, darkExp = None, exposurKey = None):
        self.path = path
        self.darkM = dark
        self.flatM = flat
        self.biasM = bias

        self.file = fits.open(self.path)
        self.darkExp = darkExp
        self.exposurKey = exposurKey
        
    def getHDU(self, i = 0):
        return self.file[i]
    
    def getInfo(self):
        return self.file.info()
    
    def getHeader(self, HDU = 0):
        return self.getHDU(HDU).header
    
    def getExposure(self, key, HDU = 0):
        return float(self.getHeader()[key])
    
    def getData(self, idx_HDU = 0):
        return np.asarray(self.getHDU(idx_HDU).data, dtype = np.float64)
    
    def getCenter(self, idx_HDU = 0):
        return (self.getShape(idx_HDU) / 2)
    
    def getShape(self, idx_HDU = 0):
        return np.asarray(self.getData(idx_HDU).shape)
    
    def getReducedData(self, HDU = 0):
        
        data = self.getData()
        
        scaleFactor = 1
        
        if self.darkExp != None and self.exposurKey != None:
            scaleFactor = float(self.getHeader()[self.exposurKey]) / self.darkExp
        
        return ((data - self.biasM) - (self.darkM)*scaleFactor) / (self.flatM / np.nanmean(self.flatM)) 

    
    def getTime(self, key, forma, HDU = 0):
        
        val = str(self.getHeader(HDU)[key])
        precision = len(val.split('.')[-1])
        
        return Time(val, format = forma, precision = precision)
        
    
    def getDataReScale(self, forma, Data = None):
        if Data == None:
            Data = self.getData()
        Data = np.asarray(Data)
        return ( np.iinfo(forma).max *  (Data - np.min(Data)) / np.max(Data) ).astype(forma)

    def setFile(self, path):
        self.path = path
        self.file = fits.open(self.path)
        
    def getTresh(self, reduced = False, display = False):
        
        if reduced:
            y,x = self.reducedHistogram()
        else:
            y,x = self.histogram()
        
        
        x = x[1:-100]
        y = np.log10(y)
        y[np.isneginf(y)] = 0
        maxx = np.max(y[5:])
        idxOfMax = np.where(y[5:] == maxx)[0][0]
        
        
        #-----------first step: compute the gaussian----------
        yGauss = np.concatenate((y[:idxOfMax + 1], y[:idxOfMax + 1][::-1])) # take the first part of the histgram gaussian and build the total gaussian by miror
        
        startingPoint = [maxx, x[idxOfMax], 1]
        
        # def func(x, a, x0, sigma):
        #     return a*np.exp(-(x-x0)**2/(2*sigma**2))
        def func(x, a, x0, sigma):
            r = np.log10(a*np.exp(-(x-x0)**2/(2*sigma**2)))
            r[r<0] = 0
            return r
        
        popt, pcov = curve_fit(func, x[:idxOfMax*2+2], yGauss, p0 = startingPoint)
        
        def linear(x, a, b):
            return a*x + b
        
        def combineModel(x, a, b):
            x1 = x[:idxOfMax]
            x2 = x[idxOfMax:]
            
            comp1 = np.log10(popt[0]*np.exp(-(x1-popt[1])**2/(2*popt[2]**2)))
            comp2 = np.log10(popt[0]*np.exp(-(x2-popt[1])**2/(2*popt[2]**2)))
            comp1[comp1 < 0] = 0
            comp2[comp2 < 0] = 0
            
            comp3 = linear(x2, a, b)
            
            comp4 = comp2
            flag = comp3 > comp2
            comp4[flag] = comp2[flag] + comp3[flag]
            res= np.concatenate((comp1, comp4))
             
            res[res<0] = 0
            
            return res, flag
        
        def func2(x, a, b):
            return combineModel(x, a, b)[0]
        
        popt2, pcov2 = curve_fit(func2, x, y[1:-99])
        yn, flag = combineModel(x, *popt2)
        
        residu = (y[1:-99] - yn)[:-10000]
        
        if display:
            print("gaussian fit a, x0, sigma:", popt )
            print("linear fit a, b:", popt2)  
            print("tresh:",  x[idxOfMax + np.argmax(flag)])
            plt.figure()
            plt.plot(x, y[1:-99], label = 'hist')
            plt.plot(x, yn, label = 'fit')
            plt.legend()
            plt.show()
        
        
        
        return x[idxOfMax + np.argmax(flag)]#x[np.where(np.abs(residu) == np.max(np.abs(residu[int(popt[1])-1::])))[0][0]]
    
        
        
    
    def findStars(self, tresh = None, onReduced = False):
        
        if onReduced:
            image = self.getReducedData().astype(np.uint16)
        else:
            image = self.getData().astype(np.uint16)
       
       
        image[image> 2**16 -1] = 0
       
        if tresh == None:
            tresh =  1.5*np.median(image)
        
        (thresh, im_bw) = cv2.threshold(image, tresh, 2**8-1, cv2.THRESH_BINARY)
        
        im_bw = cv2.medianBlur(im_bw, 5)
     
        contours, hierarchy = cv2.findContours(np.uint8(im_bw),cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
        center = []
        for c in contours:
            try:
               M = cv2.moments(c)
               cX = int(M["m10"] / M["m00"])
               cY = int(M["m01"] / M["m00"])
               center.append(np.array([cX, cY]))
        
            except:
                pass
        return np.asarray(center)
    
    
    def histogram(self, idx_HDU = 0):
        return np.histogram(self.getData(idx_HDU), bins=np.linspace(0,2**16,2**16+1))
    
    def reducedHistogram(self, idx_HDU = 0):
        return np.histogram(self.getReducedData(idx_HDU).astype(int), bins=np.linspace(0,2**16,2**16+1))
    
         
         
class TimeStruct:
    def __init__(self, jd1, jd2):
        self.jd1 = jd1
        self.jd2 = jd2
        
    def __lt__(self, other):
        return other.jd1 > self.jd1 or (other.jd1 == self.jd1 and other.jd2 > self.jd2)
    
    def __gt__(self, other):
        return other.jd1 < self.jd1 or (other.jd1 == self.jd1 and other.jd2 < self.jd2)
    
    def __eq__(self, other):
        return other.jd1 == self.jd1 and other.jd2 == self.jd2
    
    def __add__(self, other):
        if isinstance(other, TimeStruct):
            return TimeStruct(self.jd1 + other.jd1, self.jd2 + other.jd2)
        elif isinstance(other, int):
            return TimeStruct(self.jd1 + other, self.jd2)
        elif isinstance(other, float):
            return TimeStruct(self.jd1 + int(other), self.jd2 + (other - int(other)))
        
    def __iadd__(self, other):
        return self + other
    
    def __truediv__(self, other):
        return TimeStruct(self.jd1 / other, self.jd2 / other)

    def __repr__(self):
        return 'time: ' + str(self.jd1 + self.jd2)