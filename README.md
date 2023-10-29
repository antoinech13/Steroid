# Steroid

Steroid is a python package dedicated to help users to developpe their own automatic or semi-automatic procedure in way to correct frame and/or to performe photometry.

  ## Summary

  

  ## Code Structur


  ## Classes Description

  in this section, we will details all classes of the code and we will details some methods that could be usefull for users

  <details>

  <summary> 
    
  ### Corrector 
  
  </summary>
  

  **Description:**

  Astronomical images from the same sequence are rarely aligned with each other. It is common to observe, at best, a drift in both the x and y directions between each image, and at worst, a field rotation. This misalignment can have several origins, but as the  main cause, we can note the type of telescope mount (equatorial or azimuthal), as well as the quality of the mechanics, the presence or absence of guiding, the alignment, a meridian flip, etc. The "Corrector" class is therefore aimed at estimating the drift and the rotation angle between each image in the same sequence. It also provides several functionalities to the user, allowing either simple position correction or direct image correction. In the case of photometric studies, it is preferable not to correct the images directly. Indeed, due to the discrete nature of an image, a rotation of it will introduce undesirable artifacts in the image. In the case of an amateur who simply wants to do astrophotography, it is possible to directly correct the image and apply interpolation to obtain an image without visible artifacts.

  **Constructor:**
  
  ***Corrector(seqManager, flatSeq = None, biasSeq = None, darkSeq = None, exposurKey = None):***
   
  The constructor of the class *Corrector* take, as input: 
  
  -  (mandotory) a sequence of images (see the data structure SeqManager)
  -  (optional) a list of path (list of string) for the flat sequence, a list of path (list of string) for the bias sequence, a list of path (list of string) for the dark sequence and a string which correspond to the fit header key of the exposure (usually exposure is store in fits header under the key EXPOSURE or EXPTIME)


  **Methods:**
  
  
   ***getImgShape(idx = 0, idx_HDU = 0):*** 
 
  -  description: return the shape of an image of the sequence
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (tuple)

  ***getImgCenter(idx_img = 0, idx_HDU = 0):***

  -  description: return the coordinate of the center of an image of the sequence
  -  input: (INT) idx_img of image in the sequence, (INT) idx_HDU in the image
  -  return: (tuple)

  ***getImg(idx = 0):***

  -  description: return an object of type Fit (see data structure Fit)
  -  input: (INT) idx of image in the sequence
  -  return: (Fit)

  ***getData(idx = 0, HDU = 0)***

  -  desciption: return the data of the raw image.
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array)

  ***getReducedData(idx = 0, HDU = 0)***

  -  desciption: return the data of the reduced image.
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array)
  
  ***getHeader(idx = 0, HDU = 0):***
  
  -  desciption: return the header of the image at idx in the sequence.
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (STRING)

  ***histogram(idx = 0 , idx_HDU = 0):***
  
  -  desciption: return the histogram of the image at idx in the sequence and at the HDU idx.
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array) histogram, (numpy.array) bin edgesarray (see numpy.histogram)

  ***getStarsListOfImg(idx):***

  -  desciption: return the position of stars detected at the iamge idx.
  -  return: (numpy.array) stars position

  ***computeImagesCorrection(offsetTreshStarsDetection = 0, treshOnReduced = False)***

  -  description: compute the drift and the angle of rotation for each images of the sequence and store them in two list of lengh = to the sequence lengh
  -  input: (FLOAT) a offset that can be add to adjust treshold value. (BOOLEAN) if the treshold should be estimated on reduced images or raw
    
  ***medDrif(idx):***

  - description: drift is estimated between all stars detected. this function return the median value of the image at idx
  - input: (INT) idx of image in the sequence
  - oupt: (array) 2d arry of drift in both axis

  ***avgDrif(idx):***

  - description: drift is estimated between all stars detected. this function return the average value of the image at idx
  - input: (INT) idx of image in the sequence
  - oupt: (array) 2d arry of drift in both axis

  
   ***medAng(idx):***

  - description: angle is estimated between all stars detected. this function return the median value of the image at idx
  - input: (INT) idx of image in the sequence
  - oupt: (FLOAT) angle of rotation

  ***avgAng(idx):***

  - description: angle is estimated between all stars detected. this function return the average value of the image at idx
  - input: (INT) idx of image in the sequence
  - oupt: (FLOAT) angle of rotation

  ***correctStarsFromRot(arrayToCorrect, idx, coefMultAngle = -1)***

  -   description: according to a given array of positions, this function will correct each position according to the drif and angle of the image idx. The coefMultAngle take 1 or -1 and only give the direction of rotation. (different value frome 1 or -1 will influence the angle of rotation)
  -   input: (2d array) array of position to correct, (INT) idx of image for which to correct, (INT) coeficient multiply to the angle
  -   return: (2d array) new position of objects located at positions arrayToCorrect, according to the image idx drift and angle.

  ***correctedImg(idx = 0, HDU_idx = 0):***

  -  desciption: return the image corrected
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array) the corrected image

  ***getSuperImg(idx_ims = None, HDU_idx = 0):***

  -  desciption: return the average combination of all images of the sequence (after correction)
  -  input: (INT) idx of image in the sequence, (INT) idx_HDU in the image
  -  return: (numpy.array) the combined image

  ***rejectBadData():***

  - description: reject all data where drift and/or angle was not found

 ***imshowstar(idx = 0):***

 -  description: method to display image at idx and show objects detected
 -  input: (INT) idx of the image in the sequence
 
***checkPatterns(idxOfImage = 0, patidx = None):***

-  description: method to display image idx and, if patterns idx (patidx) set to None, will show all patterns. If patterns idx set to a value, will only show the pattern selected
-  input: (INT) idxOfImage image idex in the sequence. (INT) patidx index of the patterns of all patterns of the image. If set to None, will show all patterns

 </details id="detector">

  <details>

  <summary> 
    
  ### Detector 
  
  </summary>
 


**Description:**

This class is dedicated to detect moving object. it's internaly stor a list of moving objects position and and other list of their speed along x and y axis. with inital poistions and speed, it's easy to determine the position of moving objects on each frames.


**Constructor:**

***Detector(imageSeq, flatSeq = None, biasSeq = None, darkSeq = None):***

 The constructor of the class *Detector* take, as input: 
  
  -  (mandotory) a list of path (list of string) of the main image sequence
  -  (optional) a list of path (list of string) for the flat sequence, a list of path (list of string) for the bias sequence and a list of path (list of string) for the dark sequence.

**Methods:**

***computeImagesDrift(offsetTreshStarsDetection = 0, treshOnReduced = False)***

  -  description: call the function computeImagesCorrection from *Corrector* than reject all bad data (without drift or/and angle values detected)
  -  input: (FLOAT) a offset that can be add to adjust treshold value. (BOOLEAN) if the treshold should be estimated on reduced images or raw

***findAsteroid(offsetTreshAstsDetection = 0, treshOnReduced = False, eps = 2):***

  -  description: will find slow moving object based on method comparing present object or not from one of the first and one of the last frame of the sequence. To help to not overdetect to much, this algorythm is helped by a convolutional neural network based on AlexNet. This method will feed a list of moving object position on the initial frame and a list a object speed along x and y axis.
  -  input: (FLOAT) a offset that can be add to adjust treshold value. (BOOLEAN) if the treshold should be estimated on reduced images or raw. (INT) epsilon which correspond to the tolerence

***fasterAst():***

  - description: return the index of the faster moving object
  - return: (INT) the index of the faster asteroid in the list

***slowestAst():***

  - description: return the index of the slowest moving object
  - return: (INT) the index of the slowest asteroid in the list

***getAstPositionAtImg(idx):***

  - description: return the position of moving objects on the image at the idx
  - input: (INT) idex of the image where to get positions
  - return: (numpy.array) array of moving objects position

***nofa():***

  - description: return the number of moving object detected
  - return: (INT) number of moving object detected

***astSpeed(idx = 0):***

  - description: return the speed of moving objects on the image at the idx
  - input: (INT) idex of the image where to get positions
  - return: (numpy.array) array of speed on x and y axis of moving objects


</details>

## Data structure description:

this section is dedicated to talk about some classes store in the file data_structurs. 

<details>

  <summary> 
    
  ### Triangle
  
  </summary>



**Description:**

This class store 3 stars and represent a triangle. this class overload the addition, substraction, division, multiplication, comparaison et also \_\_str\_\_

**Constructor:**

***Triangle(s1, s2, s3, eps = 2):***

s1, s2 and s2 are (numpy.array). eps is a tolerence used in the \_\_eq\_\_ to estimate if two triangles are equal or not


**Methods:**

***d1():***

-  description: return the eucledian distance between s1 and s2
-  return: (FLOAT) distance between s1 and s2

***d2():***

-  description: return the eucledian distance between s1 and s3
-  return: (FLOAT) distance between s1 and s3

  
***d3():***

-  description: return the eucledian distance between s2 and s3
-  return: (FLOAT) distance between s2 and s3

***buildVect():***

- desciption: build tree vector v1, v2 and v3 between (s1, s2) , (s1, s3) and (s2, s3)
- return (numpy.array, numpy.array, numpy.array) three vector v1, v2 and v3

***getRotationAngle(other):***

-  description: compute the angle between the triangle and an other one. Cauntion!!! this method do not check if both triangles are the same
-  input: (Triangle) an other triangle
-  return: (FLOAT) the angle of rotation between both tirangles


***computeDistance(other):***

-  description: return the mean distance between the triangle and an other
-  input: (Triangle) and othee triangle
-  return: (numpy.array) mean distance in x and y of both triangles

***correctRot(angle, center):***

- description: rotate the position of s1, s2 and s3 of an angle according to a center of rotation
- input: (FLOAT) angle of rotation, (numpy.array) position of the center of rotation

</details>

<details>

  <summary> 
    
  ### Pattern:
  
  </summary>



**Desciption:**

this class store Triangles as a pattern. the addition, substraction, multiplication, division, comparaison and \_\_str\_\_ are overloaded

**Constructor**

***Pattern(t1, t2, t3, t4, t5):***

t1, t2, t3, t4 and t5 are Triangle object (see the datastructure class *Triangle*)

**Methods:**

***computeDistance(other):***

-  description: compute the mean distance between two Pattern
-  input: (Pattern) an other pattern to compute distance
-  return: (numpy.array) mean distance in x and y between the two pattern

***computeAngle(other):***

- description: compute the angle of rotation between two patterns
- input: (Pattern) and other pattern
- return: (FLOAT) angle of rotation between the two pattern


***correctRot(angle, center):***

-  description: rotate t1, t2, t3, t4 and t5 of an angle according to a center of rotation
-  input: (FLOAT) angle of rotation, (numpy.array) position of the center of rotation

</details>

<details>

  <summary> 
    
  ### SeqManager:
  
  </summary>



**Description:**

this class store list of images path from a same sequence

**Constructor:**

***SeqManager(seq):"""

seq is just a list of path of raw images (STRING)

**Methods:**

***getPath(idx):***

-  description: return the path of the image at idx
-  input: (INT) idx of the image of interest in the sequence
-  return: (STRING) return the path of the image

***getFileName(idx):***

-  description: return the name of the image at idx
-  input: (INT) the idx of the image of interest in the sequence
-  return: (STRING) the name of the image

***getImg(idx = 0):***

-  description: return on object *Fit* of the image at idx
-  input: (INT) idx of the image of interest
-  return: (Fit) a data structure of type Fit

***getHDU(idx = 0, HDU = 0):***

-  description: return the HDU of the image idx
-  input: (INT) idx of the image. (INT) HDU index
-  return (astropy.io.fits.hdu.image.PrimaryHDU) HDU of the image idx

***getInfo(idx = 0):***

- description: display info of the image at idx
- input: (INT) idx of the image of interest

***getHeader(idx = 0, HDU = 0):***

-  description: retunr the header at the HDU and at the image idx
-  input: (INT) idx of image, (INT) index of the HDU of the image at idx
-  return: (astropy.io.fits.header.Header) header of the image idx at the hdu

***getExpo(idx, key, HDU = 0):***

-  description: return the exposure from the header of the image at idx and hdu. the exposure is determine according to the key
-  input: (INT) idx of the image of interest, (STRING) key in the header corresponding to the exposure, (INT) hdu index
-  return: (FLOAT) exposure

***getData(idx = 0, idx_HDU = 0):***

-  description: return the image idx and idx of hdu as array
-  input: (INT) index of the image of interest, (INT) idex of HDU
-  return: (numpy.array) the image

***getCenter(idx_img = 0, idx_HDU = 0):***

-  description: return the center of an image at idx and of HDU
-  input: (INT) image idex, (INT) image HDU
-  return: (numpy.array) coordinate of the center of the image at idx and at hdu

***getImgShape(idx = 0, idx_HDU = 0):***

-  description: return the shape of the image at idx and at hdu
-  intput: (INT) index of the image of interest. (INT) HDU index
-  return: (TUPLE) image shape

***getTime(key, forma, idx = 0, HDU = 0):***

-  description: get the time of the image at idx and hdu from the header using the key and forma. if in the header, the time is store as julian day, (exemple: JD=2458780) then key = JD and forma=JD. For more format, refere to Time.FORMATS from astropy.time
-  input: (STRING) key of the time in the header, (STRING) format of the time in the header (refere to Time.FORMATS from astropy.time), (INT) idx of the image, (INT) idx of HDU
-  return: (astropy.time.core.Time) time of the image

***pop(idx = -1):***

-  description: delete an image at the position idx from the sequence. by default idx is set to -1 so the last image is delete
-  input: (INT) idx of image to delete

***histogram(idx = 0, idx_HDU = 0):***

-  description: return the histogram of the image idx and hdu at idx_hdu
-  input: (INT) index of the image of interest, (INT) idx of the HDU
-  return: (numpy.array) histogram, (numpy.array) bin edgesarray (see numpy.histogram) 

</details>

## Photometry:


### Description:

In the context of our studies on asteroids from the main belt, we have used the various functions of Steroid to design our own photometric software capable of automatic or semi-automatic processing. With this new tool, we have significantly increased our ability to produce light curves. This software is included in Steroid in the photometry.py file. In this section, we will show how to use it.


<details>

  <summary> 
    
  ### Methods:
  
  </summary>


**Constructor:**

***Photometry(detector = None)***

Photometry take only one optional paramters of type (Detector). Why optional? because *Photometry* also include some functions to save photometry but also some function to load. Indeed if users want to rework on some lightcurves already processed they don't need to redo all the work. *Photometry* can reload previous lightcurves. In this kind of situation, user don't need any *Detector* as the photometry was already done. He just need to build an empty *Photometry* object and use the method ***readCsv(path)***.

**Methods:**

***start(nbOfStars, center = True, maxVal = 30000, starPassageOfs = 15000)***

- description: launch the photometry according to some input parameters
- input: (INT) number of reference stars (only in case of automatic procedure), (BOOLEAN) center or not appertures of center of brightness, (FLOAT) maximum value that automatic reference stars automaticly selected should not overstep. (FLOAT) the threshold to detect stars in the context of stars passages

***plotDif(refS = 0, ast = -1, yRange = None, binning = 1, resc = True, forma = 'jd', xtick = None, inMag = True, rmExtremPoint = False, cStd = 2, deg = 4, displayRmFit = False, starPassage = False)*** 

- description: Perform plot of differential photometry
- input: - (INT) refS is the index of the star selected as reference \
  &emsp;&emsp;&ensp; - (INT) ast is, in the case of multiple asteroids, the index of the asteroid that we want to plot. \
  &emsp;&emsp;&emsp;&ensp;  If set to -1 all asteroids will be plot  \
  &emsp;&emsp;&ensp; - (list) yRange range of y axis \
  &emsp;&emsp;&ensp; - (INT) binning. use to bin lightcurve. automaticly choose if set to -1 \
  &emsp;&emsp;&ensp; - (BOOLEAN) resc. rescale stars lightcurves close to asteroid's lightcurves \
  &emsp;&emsp;&ensp; - (STRING) forma. format of the time.  refere to Time.FORMATS from astropy.time \
  &emsp;&emsp;&ensp; - (array) xticks. new x ticks \
  &emsp;&emsp;&ensp; - (BOOLEAN) inMag. if True, y axis display in magnitude. if False. Y axis display in instrumental flux. \
  &emsp;&emsp;&ensp; - (BOOLEAN) rmExtremPoint. if True, will remove extreme points. To remove extrem points, the \
  &emsp;&emsp;&emsp;&ensp; algorythm will fit a polynome then normalise asteroid's lightcurves with the polynome. \
  &emsp;&emsp;&emsp;&ensp;  each points out of [median - C x Std, median + C x Std] are removed \
  &emsp;&emsp;&ensp; - (FLOAT) cStd. this correspond de C. \
  &emsp;&emsp;&ensp; - (INT) deg. Degree of the polynome. \
  &emsp;&emsp;&ensp; - (BOOLEAN) displayRmFit. if True, display more plot to monitore rmExtremPoint. \
  &emsp;&emsp;&ensp; - (BOOLEAN) starPassage. if True, will remove star's passages \

***toDat(path, filename, binning = 1, forma = 'mjd', refS = -1, deg = 4, cStd = 2, displayRmFit = False)***

-  description: write files with extention .1, .2, .3 and .4. for each of them, the first column is the time. for others columns .1 correspond to data in instrumatal flux, .2 correspond to data in magnitude, .3 correspond to differential photometry and .4 correspond to differential photometry with averaged references star.
- input: (STRING) path correspond to the path where to save those files. (STRING) correspond to the name to give to files. Other parameters are the same than ***plotDif***.

***log(path, name = "log.txt")***
  
- description: will write a log file with information on: - data rejected, - star passages data, fwhm detected on each frames...
- input: (STRING) path to save the log file. (STRING) name give to the log file. do not forget the extention


***toGif(path)***

- description: toGif write a .gif image of all frames with appertures.
- input: (STRING) path + file name with extention (ex: r"C:/.../myGif.gif")

***toCsv(path)***

- description: write a csv file which summery the photometry. It's can be use as a backup with the method ***readCsv*** (see below).
- input: (STRING) path + file name with extention (ex: r"C:/.../myCsv.csv")
  
***readCsv(path)***   

- description: load a csv file produced with the method ***toCsv***. can be use to rework plots
- input: (STRING) path + file name with extention (ex: r"C:/.../myCsv.csv")

</details>


<details>

  <summary> 
    
  ### How to use ?
  
  </summary>




**First step:**

The first thing to use photometry is to import it:

    from photometry import Photometry

*Photometry* object constructor take, as an optional parameter, an object *Detector*. 

If the photometry was not yet done, then the users need to provide an object *Detector*. It is therefore important to include it:

    from detector import Detector


**Second step:** 

The next step is to build an object *Detector* (constructor is describ in the [Detector](#detector) section). To do this, we will use glob

     import glob

an exemple of a piece of code that can be use to build an object *Detector*:

  ~~~

#-----------set up all list of path for the raw data and bias, dark and flat data--------------

    path = r"C:\...\directory_of_your_data/"
   
    seq = glob.glob(path + "*target_repetable_name_pattern*.f*t*") #return a list of path of all file which contain target_repetable_name_pattern in their name
 
    dark = glob.glob(path + "*dark_repetable_name_pattern*.f*t*") #return a list of path of all file which contain dark_repetable_name_pattern in their name
    flat = glob.glob(path + "*flat_repetable_name_pattern*.f*t*") #return a list of path of all file which contain flat_repetable_name_pattern in their name
    bias = glob.glob(path + "*bias_repetable_name_pattern*.f*t*") #return a list of path of all file which contain bias_repetable_name_pattern in their name

#-----------------------check if bias, dark and flat data was found-----------------------------

    if len(dark) == 0:        # Check if dark data was found. If not, set up dark variable to the optional default value of Detector attribute darkSeq 
        dark = None
        print('DARK EMPTY')
    if len(bias) == 0:        # Check if bias data was found. If not, set up dark variable to the optional default value of Detector attribute biasSeq 
        bias = None
        print('BIAS EMPTY')
    if len(flat) == 0:        # Check if flat data was found. If not, set up dark variable to the optional default value of Detector attribute flatSeq 
        flat = None
        print('FLAT EMPTY')

#-------------------construct Detector object----------------------------------------------------

    d = Detector(seq, flatSeq = flat, biasSeq = bias , darkSeq = dark)
  ~~~

**Third step:**

The next step is to launch *Detector* methods to correct images and to detect asteroids:

    d.computeImagesDrift(offsetTreshStarsDetection = 0, treshOnReduced = False)
    d.findAsteroid(offsetTreshAstsDetection = 0, treshOnReduced = False, eps = 2)

***computeImagesDrift*** could me internally called in ***findAsteroid*** but we choose to let it like this to give more flexibility to users. Indeed, In this process, ***computeImagesDrift*** will take more time has it has to detect stars from all images of the sequence. Also it's the proccesse the less sensitive to the detection treshold. Indeed, it's only need 5 common stars on each frames to be able to correct images. According to this, users can earn time of execution setting up a high value of **offsetTreshStarsDetection**. For the same reason, **treshOnReduced** can be set to False.  

On the other hand, ***findAsteroid*** is really sensityve to the data quality and to the treshold. More close to the optimal value the treshold will be and better the algorythm will perform. Therefor **offsetTreshStarsDetection** should be small for small adjustement and **treshOnReduced** should be set to True.

**Fourth step**

The final step is to build an object *Photometry* and to launch the photometric processe.

    phot = Photometry(d)
    phot.start(nbOfStars, center = True, maxVal = 30000, starPassageOfs = 15000)

**nbOfStars** is mandatory for now but it's only used in the case that you choose automatic reference stars selection. It's the number of reference stars that the code will search. 

**center** is a boolean and is a paramter to allow the code to center appertures on the "center of intensity" (in reference to the center of masse equation where the masse was changed by the intensity of pixels) or not. To be clear, apperture position will all time be set according to the initial positions, to the image drift and angle of rotation and for moving object, to the speed. But, if center is set to true, after placing all appertures, the code will simply perform a centring. In case of starpassages, the apperture will probably stay fixe on the stars the time that the astroid will pass in front but after, it will come back centred on the asteroid when the star will leave the apperture feild. Moreover, with the algorythm set up to delete stars passages, all the time where the apperture will stay focuse on the star will not be present on the final lightcurve.

**maxVal** correspond to the maximum value that reference stars automaticly select should not overstep.

**starPassageOfs** have the same function as **treshOnReduced** for objects and methodes dedicated to detect objects. This one is dedicated to stars Passages detection. It's, by default, set to 15000, which detect only bright stars but can be set much lower to detect fainter stars.



**Conclusion**

The final code should look like this:

 ~~~

from photometry import Photometry
from detector import Detector

import glob

#-----------Set up all list of path for the raw data and bias, dark and flat data--------------

path = r"C:\...\directory_of_your_data/"
   
seq = glob.glob(path + "*target_repetable_name_pattern*.f*t*") 

dark = glob.glob(path + "*dark_repetable_name_pattern*.f*t*") 
flat = glob.glob(path + "*flat_repetable_name_pattern*.f*t*")
bias = glob.glob(path + "*bias_repetable_name_pattern*.f*t*") 


if len(dark) == 0:      
   dark = None
   print('DARK EMPTY')
if len(bias) == 0:        
   bias = None
   print('BIAS EMPTY')
if len(flat) == 0:        
   flat = None
   print('FLAT EMPTY')

#-------------------Construct Detector object----------------------------------------------------

d = Detector(seq, flatSeq = flat, biasSeq = bias , darkSeq = dark)

#----------------Estimate images correction and perform asteroids detection----------------------

d.computeImagesDrift(offsetTreshStarsDetection = 0, treshOnReduced = False)
d.findAsteroid(offsetTreshAstsDetection = 0, treshOnReduced = False, eps = 2)

#---------------Construct Photometry object and launch the photometry----------------------------

phot = Photometry(d)
phot.start(nbOfStars, center = True, maxVal = 30000, starPassageOfs = 15000)

 ~~~

**Caution:** in case where users use semi-automatic procedures, the selection on images are done with the left click and when selection is finish press the right click.


</details>
