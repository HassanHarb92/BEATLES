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

    elif (n == 14):
      print "Let me take you down\b'Cause I'm going to Strawberry Fields\bNothing is real\nAnd nothing to get hung about\nStrawberry Fields forever\n"
    
    elif (n == 15):
       print "Living is easy with eyes closed\nMisunderstanding all you see\nIt's getting hard to be someone\nBut it all works out\nIt doesn't matter much to me\n"

    elif (n == 16):
       print "No one I think is in my tree\nI mean it must be high or low\nThat is you can't, you know, tune in\nBut it's all right\nThat is, I think, it's not too bad\n"

    elif (n == 17):
       print "In Penny Lane, there is a barber showing photographs\nOf every head he's had the pleasure to know\nAnd all the people that come and go\nStop and say, Hello\n"

    elif (n == 18):
       print "On the corner is a banker with a motorcar\nAnd little children laugh at him behind his back\nAnd the banker never wears a mac\nIn the pouring rain, very strange\n"

    elif (n == 19):
       print "Penny Lane is in my ears and in my eyes\nThere beneath the blue suburban skies\n"

    elif (n == 20):
       print "In Penny Lane there is a fireman with an hourglass\bAnd in his pocket is a portrait of the Queen\bHe likes to keep his fire engine clean\nIt's a clean machine\n"  

    elif (n == 21):
       print "In Penny Lane, the barber shaves another customer\nWe see the banker sitting waiting for a trim\nAnd then the fireman rushes in\nFrom the pouring rain, very strange\n"
 
    elif (n == 22):
       print "I read the news today, oh boy\nAbout a lucky man who made the grade\nAnd though the news was rather sad\nWell, I just had to laugh\nI saw the photograph\n"

    elif (n == 23):
       print "He blew his mind out in a car\nHe didn't notice that the lights had changed\nA crowd of people stood and stared\nThey'd seen his face before\nNobody was really sure if he was from the House of Lords\n"

    elif (n == 24):
       print "I saw a film today, oh boy\nThe English Army had just won the war\nA crowd of people turned away\nBut I just had to look\nHaving read the book\nI'd love to turn you on\n"

    elif (n == 25):
       print "You think you've lost your love\nWell, I saw her yesterday\nIt's you she's thinking of\nAnd she told me what to say\n"
  
    elif (n == 26):
       print "She says she loves you\nAnd you know that can't be bad\nYes, she loves you\nAnd you know you should be glad\n"

    elif (n == 27):
       print "She said you hurt her so\nShe almost lost her mind\nBut now she said she knows\nYou're not the hurting kind\n"

    elif (n == 28):
       print "You don't need me to show the way, love\nWhy do I always have to say, love\nCome on, come on, come on, come on\nCome on, come on, come on, come on\nPlease, please me, whoa yeah, like I please you\n"

    elif (n == 29):
       print "I don't want to sound complaining\nBut you know there's always rain in my heart\nI do all the pleasing with you\nIt's so hard to reason with you\nwhoa yeah, why do you make me blue?\n"

    elif (n == 30):
       print "Dear Prudence, won't you come out to play?\nDear Prudence, greet the brand new day\n"
  
    elif (n == 31):
       print "Dear Prudence, open up your eyes\nDear Prudence, see the sunny skies\nThe wind is low, the birds will sing\nThat you are part of everything\nDear Prudence, won't you open up your eyes?\n"

    elif (n == 32):
       print "Dear Prudence, let me see you smile\nDear Prudence, like a little child\nThe clouds will be a daisy chain\nSo let me see you smile again\nDear Prudence, won't you let me see you smile?\n"

    else:
       print "Imagine there's no countries..."



# Songs covered:
#
#  Following up on this list:
#  
#  https://www.vulture.com/2017/06/all-213-beatles-songs-ranked-from-worst-to-best.html
#
#  1. Imagine
#  2. Come together
#  3. A Hard Day's Night 
#  4. Hey Jude
#  5. Strawberry Fields Forever
#  6. Penny Lane
#  7. A Day in the Life
#  8. She Loves You
#  9. Please Please Me
#  10. Dear Prudence
#  
