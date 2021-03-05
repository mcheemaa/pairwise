
from load import *
import numpy as np
import argparse


mainX = []
mainY = []
gapPen = -1

def alignment (firstSeq, secondSeq, checki):

    matrixC = np.zeros(((len(firstSeq) + 1), (len(secondSeq) + 1)))
    position = np.zeros(((len(firstSeq) + 1), (len(secondSeq) + 1)))


    position[0][0] = None

    for i in range (1, (len(firstSeq) + 1)):
        matrixC[i][0] = i * gapPen
        position[i][0] = None
    for j in range (1, (len(secondSeq)+ 1)):
        matrixC[0][j] = j * gapPen
        position[0][j] = None

    for i in range (1, (len(firstSeq) + 1)):
        for j in range (1, (len(secondSeq)+ 1)):
            diag = matrixC[i - 1][j - 1] + int(matchScore(firstSeq[i - 1], secondSeq[j - 1], checki))
            left = matrixC[i - 1][j] + gapPen
            top = matrixC[i][j - 1] + gapPen

            matrixC [i][j] = max(diag, top, left)

            # diag = 1 match , left = -1, top = -2

            if (matrixC[i][j] == diag):
                position[i][j] = 1

            elif (matrixC[i][j] == left):
                position[i][j] = -1

            elif (matrixC[i][j] == top):
                position[i][j] = -2

    fseq = ""
    sseq = ""

    i = len(firstSeq) 
    j = len(secondSeq) 

    while i > 0 and j > 0:

        if (position[i][j] == 1):
            fseq = fseq + firstSeq[i - 1]
            sseq = sseq + secondSeq[j - 1]
            i = i - 1
            j = j - 1

        elif (position[i][j] == -1):
            fseq = fseq + firstSeq[i - 1]
            sseq = sseq + '_'
            i = i - 1

        elif (position[i][j] == -2):
            sseq = sseq + secondSeq[j - 1]
            fseq = fseq + '_'
            j = j - 1

        '''cur = matrixC[i][j]
        curleft = matrixC[i - 1][j]
        curtop = matrixC[i][j - 1] 
        curdiag = matrixC[i - 1][j - 1]


        if cur == curdiag + int(matchScore(firstSeq[i - 1], secondSeq[j - 1], checki)):
            fseq = fseq + firstSeq[i - 1]
            sseq = sseq + secondSeq[j - 1]
            i = i - 1
            j = j - 1

        elif cur == curleft + gapPen:
            fseq = fseq + firstSeq[i - 1]
            sseq = sseq + '_'
            i = i - 1

        elif cur == curtop + gapPen:
            sseq = sseq + secondSeq[j - 1]
            fseq = fseq + '_'
            j = j - 1'''

    while i > 0:
        fseq = fseq + firstSeq[i - 1]
        sseq = sseq + '_'
        i = i - 1

    while j > 0:
        sseq = sseq + secondSeq[j - 1]
        fseq = fseq + '_'
        j = j - 1

   

    return [fseq[::-1], sseq[::-1], matrixC[-1][-1]]

def outfunc(align1, align2, best):

    score = 0
    identity = 0

    for i in range (len(align1)):

       
        if (align1[i] == align2[i] and align1[i] != "_" and align2[i] != "_"):

            identity = identity + 1

        
    print(best)
    print((identity/(len(align1))) * 100)



def matchScore(first, second, checki):

    matrixScorew(checki)


    firstcoord = 0
    secondcoord = 0

    for i in range (len(mainX)):
        if (first == mainX[i]):
            firstcoord = i
    
    for j in range (len(mainX)):
        if (second == mainX[j]):
            secondcoord = j
    
    return(mainY[firstcoord][secondcoord + 1])
    



def matrixScorew (check):

    if (check == "F"):
        matFile = "dnaMatrix2"
    elif (check == "T"):
        matFile = "BLOSUM45"


    counter = 0
    with open(matFile) as f:
        lineList = f.readlines()
            
        for lines in lineList:
            tempList = lines.split()
               
            if (counter == 2):
                global mainX
                mainX = tempList
            elif (counter > 2 and counter < (len(lineList) - 1)):
                global mainY
                mainY.append(tempList)
            elif (counter == len(lineList) - 1):
                global gapPen
                #print(tempList[0])
                gapPen = int(tempList[0])

            counter = counter + 1


def semig (firstSeq, secondSeq, checki):


    matrixC = np.zeros(((len(firstSeq) + 1), (len(secondSeq) + 1)))
    position = np.zeros(((len(firstSeq) + 1), (len(secondSeq) + 1)))


    position[0][0] = None

    for i in range (1, (len(firstSeq) + 1)):
        matrixC[i][0] = 0
        position[i][0] = None
    for j in range (1, (len(secondSeq)+ 1)):
        matrixC[0][j] = 0
        position[0][j] = None




    for i in range (1, (len(firstSeq) + 1)):
        for j in range (1, (len(secondSeq)+ 1)):
            diag = matrixC[i - 1][j - 1] + int(matchScore(firstSeq[i - 1], secondSeq[j - 1], checki))
            
            if (i == len(firstSeq) or j == len(secondSeq)):
                left = matrixC[i - 1][j] 
                top = matrixC[i][j - 1]
            else:
                left = matrixC[i - 1][j] + gapPen
                top = matrixC[i][j - 1] + gapPen

            matrixC [i][j] = max(diag, top, left)

            # diag = 1 match , left = -1, top = -2, 


            if (matrixC[i][j] == diag):
                position[i][j] = 1

            elif (matrixC[i][j] == left):
                position[i][j] = -1

            elif (matrixC[i][j] == top):
                position[i][j] = -2


    fseq = ""
    sseq = ""

    i = len(firstSeq) 
    j = len(secondSeq) 

    while i > 0 and j > 0:

        if (position[i][j] == 1):
            fseq = fseq + firstSeq[i - 1]
            sseq = sseq + secondSeq[j - 1]
            i = i - 1
            j = j - 1

        elif (position[i][j] == -1):
            fseq = fseq + firstSeq[i - 1]
            sseq = sseq + '_'
            i = i - 1

        elif (position[i][j] == -2):
            sseq = sseq + secondSeq[j - 1]
            fseq = fseq + '_'
            j = j - 1

        

    while i > 0:
        fseq = fseq + firstSeq[i - 1]
        sseq = sseq + '_'
        i = i - 1

    while j > 0:
        sseq = sseq + secondSeq[j - 1]
        fseq = fseq + '_'
        j = j - 1

    #print(position)
    #print(matrixC)
    #outfunc(fseq[::-1], sseq[::-1], matrixC[-1][-1])
    #print(fseq[::-1])
    #print(sseq[::-1])

    return [fseq[::-1], sseq[::-1], matrixC[-1][-1]]


def local (firstSeq, secondSeq, checki):

    #print(len(firstSeq))

    matrixC = np.zeros(((len(firstSeq) + 1), (len(secondSeq) + 1)))
    position = np.zeros(((len(firstSeq) + 1), (len(secondSeq) + 1)))


    #print(gapPen)
    #print(matrixC)

    #print(len(firstSeq))

    position[0][0] = None

    for i in range (1, (len(firstSeq) + 1)):
        matrixC[i][0] = 0
        position[i][0] = None
    for j in range (1, (len(secondSeq)+ 1)):
        matrixC[0][j] = 0
        position[0][j] = None

    maxScore = 0
    coordI = 0
    coordJ = 0


    for i in range (1, (len(firstSeq) + 1)):
        for j in range (1, (len(secondSeq)+ 1)):
            diag = matrixC[i - 1][j - 1] + int(matchScore(firstSeq[i - 1], secondSeq[j - 1], checki))
            left = matrixC[i - 1][j] + gapPen
            top = matrixC[i][j - 1] + gapPen
            neg = matrixC[0][0]
            matrixC [i][j] = max(diag, top, left, neg)




            if (matrixC[i][j] > maxScore):

                maxScore = matrixC[i][j]
                coordI = i
                coordJ = j



            # diag = 1 match , left = -1, top = -2, neg = -3
            #print(matrixC[i][j])

            if (matrixC[i][j] == diag):
                position[i][j] = 1

            elif ((max(diag, top, left) < 0)):
                matrixC [i][j] = neg
                position [i][j] = -3

            elif (matrixC[i][j] == left):
                position[i][j] = -1

            elif (matrixC[i][j] == top):
                position[i][j] = -2

           


    fseq = ""
    sseq = ""

    i = coordI
    j = coordJ 

    while i > 0 and j > 0:

        if (position[i][j] == 1):
            fseq = fseq + firstSeq[i - 1]
            sseq = sseq + secondSeq[j - 1]
            i = i - 1
            j = j - 1

        elif (position[i][j] == -1):
            fseq = fseq + firstSeq[i - 1]
            sseq = sseq + '_'
            i = i - 1

        elif (position[i][j] == -2):
            sseq = sseq + secondSeq[j - 1]
            fseq = fseq + '_'
            j = j - 1

        elif (position[i][j] == -3):
            fseq = fseq + firstSeq[i - 1]
            sseq = sseq + secondSeq[j - 1]
            break

    #print(position)
    #print(matrixC)
    #outfunc(fseq[::-1], sseq[::-1], matrixC[-1][-1])
    #print(fseq[::-1])
    #print(sseq[::-1])
    return [fseq[::-1], sseq[::-1], matrixC[-1][-1]]



def main():

    p = argparse.ArgumentParser()

    p.add_argument('-i', '--first', type = str)
    p.add_argument('-j', '--second', type = str)
    p.add_argument('-p', '--protein', type = str)
    p.add_argument('-atype', '--gs', type = str)
    p.add_argument('-o', '--output', type = str)


    check = p.parse_args()

    seq1 = loadSeq(check.first)
    seq2 = loadSeq(check.second)

    dnacheck = check.protein
    globalcheck = check.gs

    if (globalcheck == "G"):

        if (len(seq1) >= len(seq2)):
            lis = alignment(seq1, seq2, dnacheck)

        else:
            lis = alignment(seq2, seq1, dnacheck)

    if (globalcheck == "S"):

        if (len(seq1) >= len(seq2)):
            lis = semig(seq1, seq2, dnacheck)

        else:
            lis = semig(seq2, seq1, dnacheck)

    

    output = open(check.output, 'w')

    output.write("Sequence 1: " + seq1 + "\n")
    output.write("Sequence 2: " + seq2 + "\n")
    output.write("Alignment 1: " + lis[0] + "\n")
    output.write("Alignment 2: " + lis[1] + "\n")
    output.write("Score: " + str(lis[2]) + "\n")

    score = 0
    identity = 0

    for i in range (len(lis[0])):

       
        if (lis[0][i] == lis[1][i] and lis[0][i] != "_" and lis[1][i] != "_"):

            identity = identity + 1

    sameS = ((identity/len(lis[0]))) * 100

    output.write("Identities: " + str(identity) + " / " + str(len(lis[0])) + " " + str(sameS))

    
       
if __name__ == "__main__":
    main()
