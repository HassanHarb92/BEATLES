#!/usr/bin/python

from __future__ import division
import sys
import math
import cmath
import numpy as np
from numpy import genfromtxt
import csv
from decimal import Decimal
import os

# This file contains quotes from the BEATLES songs and/or songs written by individual members of the BEATLES
# 
# The function is called at the end of the python script and it will randomly pick one quote from the list below
# 

def LyricsLibrary(n):

    print "\nDaily Dose of Beatles!\n"
    if (n == 1):
       print "Imagine there's no heaven\nIt's easy if you try\nNo hell below us\nAbove us only sky\nImagine all the people\nliving for today\n"

    elif (n==2):
       print "You may say I'm a dreamer,\nbut I'm not the only one.\nI hope someday you'll join us.\nAnd the world will live as one.\n"

    elif (n == 3):
       print "Imagine there's no countries\nIt isn't hard to do\nNothing to kill or die for\nAnd no religion too\nImagine all the people living life in peace\n"

    elif (n == 4):
       print "Imagine no possessions\nI wonder if you can\nNo need for greed or hunger\nA brotherhood of man\n"

    elif (n == 5):
       print "Here come old flat top\nHe come groovin' up slowly\nHe got joo joo eyeballs\nHe one holy roller\nHe got hair down to his knee\nGot to be a joker\nHe just do what he please\n"
 
    elif (n == 6):
       print "He bad production\nHe got walrus gumboot\nHe got Ono sideboard\nHe one spinal cracker\nHe got feet down below his knee\n"

    elif (n == 7):
       print "He wear no shoeshine\nHe got toe jam football\nHe got monkey finger\nHe shoot Coca-Cola\nHe say I know you, you know me\n"

    elif (n == 8):
       print "One thing I can tell you is\nYou got to be free\n"

    elif (n == 9):
       print "It's been a hard day's night, and I've been working like a dog\nIt's been a hard day's night, I should be sleeping like a log\n"

    elif (n == 10):
      print "When I'm home everything seems to be right\nWhen I'm home feeling you holding me tight, tight, yeah\n"

    elif (n == 11):
      print "Hey Jude, don't make it bad\nTake a sad song and make it better\nRemember to let her into your heart\nThen you can start to make it better\n"

    elif (n == 12):
      print "And anytime you feel the pain\nHey Jude, refrain\nDon't carry the world upon your shoulders\n"      

    elif (n == 13):
      print "Hey Jude, don't let me down\nYou have found her, now go and get her\nRemember to let her into your heart\nThen you can start to make it better\n"

    else:
       print "Imagine there's no countries..."




# Songs covered:
#
#  1. Imagine
#  2. Come together
#  3. A Hard Day's Night 


