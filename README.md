# Steroid

Steroid is a python package dedicated to help users to developpe their own automatic or semi-automatic procedure in way to correct frame and/or to performe photometry.

  ## Summary

  

  ## Code Structur

  ## Classes Description

  in this section, we will details all classes of the code and we will details some methods that could be usefull for users

  ### Corrector

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


### Detector

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


## Data structure description:

this section is dedicated to talk about some classes store in the file data_structurs. 

### Triangle

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


### Pattern:


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

### SeqManager:

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

