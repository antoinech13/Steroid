# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:35:12 2022

@author: antoine
"""

import numpy as np
import pandas as pd

from corrector import *
from data_structurs import SeqManager, Appertur, TimeStruct

from astropy.time import Time


class Occult(Corrector):
    def __init__(self, seq = None, flatSeq = None, biasSeq = None, darkSeq = None):
        
        if seq != None:
            super().__init__(SeqManager(seq), flatSeq, biasSeq, darkSeq)
            print(self.getHeader())
            self.key = str(input("ENTER THE KEY OF A TIME REF WITH | SEPARATOR IF NEEDED (exemple JD | jd or DATE-OBS | isot. Refere to Time.FORMATS from astropy.time): "))
            self.format = self.key.split("|")[-1].replace(' ', '')
            self.key =  self.key.split("|")[0].replace(' ', '')
        
        
        self.apperturs = []
        self.refStars = []
        self.objectOfInterst = []
        self.results = []
    
    
    def start(self, r, ri, re, center = False, correct = False):
        
        print("SELECT STAR OF INTERREST")
        self.selectAp(r, self.objectOfInterst)
        print("SELECT REFERENCES STARS")
        self.selectAp(r, self.refStars)
        self.objectOfInterst = Utils.centred(self.objectOfInterst, self.getImg().getData(), r)
        fig,ax = plt.subplots(1)
        
        for i in range(len(self)):
            
            img = self.getImg(i)
            dat = img.getReducedData()
            
            ax.clear()
            
            ax.imshow(dat, cmap='Greys', vmin = np.median(dat[dat>0])*0.5, vmax = np.median(dat[dat>0]) / 0.5)
            
            nexRefStar = Utils.centred(self.refStars, dat, r)
            dr = self.driftByCenter(nexRefStar)
            self.refStars = nexRefStar
            
            if center:
                self.objectOfInterst = Utils.centred(self.objectOfInterst, dat, r)
            else:
                self.objectOfInterst = self.objectOfInterst - dr
           
            for j in self.refStars:
                c = Circle(j, radius = r, fill=False)
                ax.add_patch(c)
                c = Circle(j, radius = ri, fill=False)
                ax.add_patch(c)
                c = Circle(j, radius = re, fill=False)
                ax.add_patch(c)
                
            for j in self.objectOfInterst:
                c = Circle(j, radius = r, fill=False)
                ax.add_patch(c)
                c = Circle(j, radius = ri, fill=False)
                ax.add_patch(c)
                c = Circle(j, radius = re, fill=False)
                ax.add_patch(c)
                
            plt.pause(0.01)
            
            
            stars = np.concatenate((self.refStars, self.objectOfInterst))
            self.apperturs.append(Appertur(stars, r = r, ri = ri, re = re))
            self.results.append(self.apperturs[-1].Photom(self.getImg(i), self.key, self.format))
            
            print("image: ", i, '\n', self.apperturs[-1].Photom(img, self.key, self.format))
    
    def driftByCenter(self, newPos):
        newPos = np.array(newPos)
        return np.mean(np.array(self.refStars) - newPos, axis = 0)
        
    def selectAp(self, r, s):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        isClose = [False]
        img = self.getData()
        
        ax.imshow(img, cmap='Greys', vmin = np.nanmedian(img)*0.5, vmax = np.nanmedian(img)/0.5)
            
        plt.connect('button_press_event', lambda event: Utils.getPosFromImage(event, s, isClose))
        plt.show(block=False)
        
        while not isClose[0]:
            plt.pause(0.01)
            
        
    def phot(self):
        stars = []
        for res in self.results:
            star = []
            for i in range(len(res)):
                star.append(res[i][7])
               
            stars.append(star)
        return np.asarray(stars)
    
    def starsPhot(self):
        
        stars = []
        for res in self.results:
            star = []
            for i in range(len(self.refStars)):
                star.append(res[len(self.objectOfInterst) + i][7])
            stars.append(star)
        return np.asarray(stars)
                
        
    def astPhot(self):
        
        asts = []
        for res in self.results:
            ast = []
            for i in range(len(self.objectOfInterst)):
                ast.append(res[i][7])
               
            asts.append(ast)
        return np.asarray(asts)
    

    def extractTime(self):
        
        # x = np.zeros((len(self)))
        # for i in range(len(self)):
        #     x[i] = self.getImg(i).getTime(self.key)
        # return x
        
        x = []
        for i in range(len(self.results)):
            t = Time(self.results[i][0][4], format = 'jd', precision = 9)
            x.append(TimeStruct(t.jd1, t.jd2))
        return np.asarray(x)
    
    # def timeFormat(self, tArr, forma):
    #     newArr = []
    #     for i in tArr.to_value(forma):
    #         if forma == 'hms':
    #             newArr = tArr.to_value('isot').split('T')[-1]
    #         else:
    #             newArr = tArr.to_value(forma)
    #     return newArr
    
    def timeFormat(self, tArr, forma):
        newArr = np.zeros(tArr.shape, dtype= object)
        for j in range(tArr.shape[0]):
            if forma == 'hms':
                newArr[j] = Time(tArr[j].jd1, tArr[j].jd2, format = 'jd', precision = 9).to_value('isot').split('T')[-1]
            else:
                newArr[j] = Time(tArr[j].jd1, tArr[j].jd2, format = 'jd', precision = 9).to_value(forma)
                
        return newArr

    def plot(self, yRange = None, binning = 1, astsSelection = None, starsSelection = None, inMag = False, forma = 'jd', xtick = None):
      
        if starsSelection is not None and (isinstance(starsSelection, int) or isinstance(starsSelection, float)):
            starsSelection = [starsSelection]
            
        if astsSelection is not None and (isinstance(astsSelection, int) or isinstance(astsSelection, float)):
            astsSelection = [astsSelection]
        
        plt.figure()
        if inMag:
            stars = -2.5*np.log10(Utils.binn(binning, self.phot()))
        else:
            stars = Utils.binn(binning, self.phot())
        print(stars.shape)
        
        x = Utils.binn(binning, self.extractTime())
        x = self.timeFormat(x, forma)
        
        for i in range(stars.shape[1]):
            if i < len(self.refStars):
                if starsSelection is not None:
                    if i in starsSelection:
                        plt.scatter(x, stars[:, i], label = 'star ' + str(i))
                else:
                    plt.scatter(x, stars[:, i], label = 'star ' + str(i))
            else:
                if astsSelection is not None:
                    if i in  astsSelection:
                        plt.scatter(x, stars[:, i], label = 'occulted ' + str(i))
                else:
                    plt.scatter(x, stars[:, i], label = 'occulted ' + str(i))
        
        if inMag:
            plt.gca().invert_yaxis()
   
        if xtick != None:
           plt.xticks(np.linspace(plt.xlim()[0], plt.xlim()[-1], xtick))
          
        if yRange != None:
            plt.ylim(yRange)
            
        plt.legend()
        plt.show()
        
    def plotDif(self, refS = 0, ast = -1, yRange = None, binning = 1, resc = True, forma = 'jd', xtick = None, inMag = True):  
        if inMag:
            print("result in Mag")
            asts = -2.5*np.log10(Utils.binn(binning, self.astPhot()))
            stars = -2.5*np.log10(Utils.binn(binning, self.starsPhot()))
        else:
            print("result in ADU")
            asts = Utils.binn(binning, self.astPhot())
            stars = Utils.binn(binning, self.starsPhot())
       
        
        x = Utils.binn(binning, self.extractTime())
        x = self.timeFormat(x, forma)
        
        if ast == -1:
            plt.figure()
            for i in range(asts.shape[1]):
                plt.scatter(x, asts[:,i] - stars[:, refS], label='ast - s' + str(refS+1), linewidths=1, marker='x')
        else:
            plt.figure()
            plt.scatter(x, asts[:,ast] - stars[:, refS], label='ast - s' + str(refS+1), linewidths=1, marker='x')
        
        of = np.zeros(stars.shape[1])
        if resc:
            of = Utils.rescalMulti((asts.T  - stars[:, refS]).T, (stars.T - stars[:, refS]).T)
        
        for j in range(stars.shape[1]):
            if j != refS:
                plt.scatter(x, (stars[:,j] - stars[:, refS]) + of[j], label='s' + str(j+1) +' - s' + str(refS+1), linewidths=1, marker='x')
       
        plt.legend()        
        if yRange != None:
            plt.ylim(yRange)

        if xtick != None:
            plt.xticks(np.linspace(plt.xlim()[0], plt.xlim()[-1], xtick))
        
        if inMag:
            plt.gca().invert_yaxis()
            
        plt.show()
    
    def toAOTACsv(self, pName, forma = 'jd'):
        
        if len(self.results) == 0:
            print('empty result list')
            return 
        
        x = Utils.binn(1, self.extractTime())
        # x = Time(x, format = 'jd', precision = 9)
        x = self.timeFormat(x, forma)
        
        header = ["FrameNo", "Time (UT)"]
        
        for i in range(len(self.results[0])):
            header.append("Signal ("+ str(i+1) + ")")
            header.append("Background ("+ str(i+1) +")")
        
        resultats = []
        for i, res in enumerate(self.results):
            resultat = [i,x[i]]
            
            for re in res:
                resultat.append(re[3])
                resultat.append(re[-2])
            
            resultats.append(resultat)
        
        df = pd.DataFrame(resultats, columns= header)
        df.to_csv(pName, index=False)
        
        
    def toCsv(self, path):
    
        header = self.buildCsvHeader()
        
        csv = []
        for table in self.results:
            csvRow = np.zeros(( len(table) * len(table[0]) ), dtype = object)
            for i, row in enumerate(table):
                if i < len(self.objectOfInterst):
                    csvRow[len(table)-1-i] = 1
                else:
                    csvRow[len(table)-1-i] = 0
                
                for j in range(len(row)-1):
                    if j == 0 or j == 1:
                        csvRow[(j+1)*len(table) + i] = row[1+j].value
                    # elif j == 3:
                    #     csvRow[(j+1)*len(table) + i] = row[1+j]
                    else:
                        csvRow[(j+1)*len(table) + i] = row[1+j]
            csv.append(csvRow)
            
        df = pd.DataFrame(csv, columns= header)
        df.to_csv(path, index=False)
            
        return df
    
    def buildCsvHeader(self):
        
        NOfAp = len(self.results[0])
        NbOfRowElements = len(self.results[0][0])
        header = [0]*NOfAp*NbOfRowElements
        
        for i in range(NOfAp):
            header[i] = 'isOcullt (' + str(i+1) + ')'
    
        for i in range(NbOfRowElements-1):
            for j in range(NOfAp):
                if i == 0:
                    header[NOfAp + i*NOfAp + j] = "xcenter (" + str(j+1) + ')'
                elif i == 1:
                    header[NOfAp + i*NOfAp + j] = "ycenter (" + str(j+1) + ')'
                elif i == 2:
                    header[NOfAp + i*NOfAp + j] = "aperture_sum (" + str(j+1) + ')'
                elif i == 3:
                    header[NOfAp + i*NOfAp + j] = "Time (" + str(j+1) + ')'
                elif i == 4:
                    header[NOfAp + i*NOfAp + j] = "annulus_median (" + str(j+1) + ')'
                elif i == 5:
                    header[NOfAp + i*NOfAp + j] = "aper_bkg (" + str(j+1) + ')'
                elif i == 6:
                    header[NOfAp + i*NOfAp + j] = "aper_sum_bkgsub (" + str(j+1) + ')'
            
            
        return header
   
    def readCsv(self, path):
       data = pd.read_csv(path)
       
       header = data.columns.values
       qTableHeader = ['id','xcenter','ycenter','aperture_sum','Time','annulus_median','aper_bkg','aper_sum_bkgsub']


       NofAp = len([i for i in header if 'isOcullt' in i])
       Nast = len([i for i in data.iloc[0,:NofAp] if i == 1])
       
       self.results = []
       
       for i in range(data.shape[0]):
           qTArr = np.zeros((NofAp, len(qTableHeader)), dtype = object)           
           
           for row in range(qTArr.shape[0]):
               for col in range(qTArr.shape[1]):
                   
                 
                   if col == 0:
                       qTArr[row, col] = row + 1
                   else:
                       qTArr[row, col] = data.iloc[i, col*NofAp+row] 
                       
                       
                   # elif col == 1:
                   #     qTArr[col, row] = data.iloc[i, NofAp+row] 
                   # elif col == 2:
                   #     qTArr[col, row] = data.iloc[i, 2*NofAp+row]
                   # elif col == 3:
                   #     qTArr[col, row] = data.iloc[i, 3*NofAp+row]
           self.results.append(QTable(qTArr, names = qTableHeader))
        
        
if __name__ == '__main__':

    import glob
    # ----------------ex occult-----------------
    path = 'C:\\Users\\antoine\\OneDrive\\Documents\\PHD\\lightcurve\\entrainement\\Adorea\\reduced\\target-c1\\'
    path = 'C:\\Users\\antoi\\OneDrive\\Documents\\PHD\\lightcurve\\entrainement\\22-04-12PST3c/'
    res = 'C:\\Users\\antoi\\OneDrive\\Documents\\PHD\\lightcurve\\entrainement\\PST3c\\res\\'
    seq = glob.glob(path + "*1999*.fit*")
    
    # bias = glob.glob(path + "*bias*.fit*")
    # dark = glob.glob(path + "*dark*.fit*")

    
    occult = Occult(seq)
    r = 3
    occult.start(r, r*1.5, r*2, False)
    