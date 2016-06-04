#!/bin/python

import sys
import subprocess

def generate_matrix():
    rts = ""
    exe = sys.argv[2]
    args = " --generate " + sys.argv[2]
    matrix_size = int(sys.argv[2])
    proc = subprocess.Popen([exe + args], stdout=subprocess.PIPE, shell=True)
    result, err = proc.communicate()
    file = open("matrix.dat", "w")
    file.write(result.decode("ascii"))
    file.close()
    result = result.decode("utf-8")
    row = 2
    rts+= "1 "
    for c in result:
        if c != '\n':
            rts += c + " "
        else:
            if row <= matrix_size:
                rts += c + str(row) + " "
                row+=1
    return rts


content = "set V:="
size_matrix = int(sys.argv[2])
for i in range(1, size_matrix+1):
    content+=str(i) + " "
content+= ";\n"
content+= "param G:  "
for i in range(1, size_matrix+1):
    content+=str(i) + " "
content+= ":=\n"
content+= generate_matrix()

print(content)

