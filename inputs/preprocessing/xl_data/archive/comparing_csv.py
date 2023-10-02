import csv
import os, sys
#for these dataset keyword used for searching the database doesnot affect the rows in the file; i.e., same xlink will be there in same order
with open('micos60.csv', 'r') as f1, open(sys.argv[1], 'r') as f2:
    fileone = f1.readlines()
    filetwo = f2.readlines()

with open('micos60.csv', 'a') as outFile:
    for line in filetwo: #if the row is absent in micos60.csv, then it will be added to it
        if line not in fileone:
            outFile.write(line)
