__author__ = "Felipe Zapata"

#==========> Standard libraries and third-party <===============
from collections import namedtuple
from pymonad     import curry
from pyparsing   import *
import numpy     as np

#==================> Internal modules <====================
from qmworks.parsers.parser import floatNumber,minusOrplus,natural,point,skipAnyChar

#====================> <=========================
