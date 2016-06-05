#!/bin/python

import sys
import subprocess

def generate_matrix():
    rts = ""
    file = open(sys.argv[1], "r")
    result = file.read()
    file.close()
    row = 2
    rts+= "1 "
    matrix_size = 0
    current_row = 0
    for c in result:
        if c != '\n':
            rts += c + " "
            if current_row == 0:
                matrix_size+=1
        else:
            if row <= matrix_size:
                rts += c + str(row) + " "
                row+=1
                current_row+=1
    return [rts, matrix_size]


matrix_data, size_matrix = generate_matrix()

content = "set V:="
for i in range(1, size_matrix+1):
    content+=str(i) + " "
content+= ";\n"
content+= "param G:  "
for i in range(1, size_matrix+1):
    content+=str(i) + " "
content+= ":=\n"
content+= matrix_data

print(content)

