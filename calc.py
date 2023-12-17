from math import atan2,sin,cos,asin

alpha = (17+13/60+28.5/3600)*15
delta = (3+32/60+40/3600)
obsphi = 52.2297
obslam = 21.0122

juliandate = 2460053.5 # days
gmst = 206.8412555028 # deg
gmst+=obslam
gmst+=alpha
gmst%=360

print(alpha)