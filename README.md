Simplified Gaia scanning law calculation 

The gaiarot modules implement the approximate Gaia scanning law and allows to determine the number of visits of a given location on the sky

Usage: 
> print (gaiarot.coveragemap(nside=2,nyears=0.2))                                 
Out[3]: 
(array([ 0,  1,  0, 13,  4,  4,  5,  3,  0,  3,  4,  2,  6,  9,  3,  3,  2,
         0,  1,  2,  2,  4,  0, 13,  1,  0,  2,  1, 12,  5,  0,  4,  3,  3,
         4,  0,  2,  8,  3,  7, 14,  1,  0,  0,  5,  4,  5,  5]),
 [array([], dtype=float64),
  array([0.5053]),
  array([], dtype=float64),
  array([57.13745, 57.38675, 57.63615, 66.11575, 66.36505, 66.6145 ,
          5.25805, 57.31295, 57.56235, 57.81185, 
          66.2913 , 66.5407 ,
the coveragemap function returns the number of times the hpx was observed and the list of arrays of observing times for each hpx

