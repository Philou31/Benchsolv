#!/bin/bash

import os
import sys
import math

def count_non_zeros(f):
	count=0
	header=1
	row=[]
	col=[]
	arrow=[]
	fstr=open(f, "r")
	while 1:
		line=fstr.readline()
		if line=="":
			break
		if line[0]=="%":
			continue
		array=line.split()
		if header:
			m=int(array[0])
			n=int(array[1])
			nz=int(array[2])
			for i in range(0, m):
				row.append(0)
			for i in range(0, n):
				col.append(0)
			for i in range(0, min(m, n)):
				arrow.append(0)
			header=0
		else:
			if not count%2000000:
				print "count: " + str(count)
			i=int(array[0])-1
			j=int(array[1])-1
			val=float(array[2])
			row[i]+=1
			col[j]+=1
			arrow[min(i, j)]+=1
			count+=1
	fstr.close()
	return m,n,nz,row,col,arrow

#COMPUTE NUMBER OF BLOCKS AND NZ
def compute_nrows(m, n, nz_chunk, row, col, arrow, current):
	# At least one block
	current_nz=nz_chunk-row[current]
	for j in range(current+1, m):
		if row[j] > current_nz:
			break
		else:
			current_nz-=row[j]
	return current_nz, j

def compute_ncols(m, n, nz_chunk, row, col, arrow, current):
	return compute_nrows(n, m, nz_chunk, col, row, arrow, current)

def compute_narrows(m, n, nz_chunk, row, col, arrow, current):
	# At least one block
	current_nz=nz_chunk-arrow[current]
	for j in range(current+1, min(m, n)):
		nz_arrow=arrow[j]
		if nz_arrow > current_nz:
			break
		else:
			current_nz-=nz_arrow
	if current_nz < 0:
		print "CURRENT NZ NEGATIVE"
	return current_nz, j

compute_nblocks = [
	compute_nrows,
	compute_ncols,
	compute_narrows
]

# OUTPUT REST IN LAST CHUNK
def last_rows_chunk(m, n):
	return m;

def last_cols_chunk(m, n):
	return n;

def last_arrows_chunk(m, n):
	return min(m, n);

output_last_blocks_chunk = [
	last_rows_chunk,
	last_cols_chunk,
	last_arrows_chunk
]

def output_chunks_per_block(m, n, nz, nz_chunk, row, col, arrow, option, output):
	total=0
	total_blocks=0
	current=0
	# DIRECTORY CREATION
	directory = "/".join(output.split("/")[:-1])
	if directory != "" and not os.path.exists(directory):
		print "Creating directory: " + directory
		os.makedirs(directory)
	# DISTRIBUTION IN CHUNKS
	f=open(output, "w")
	for i in range(0, nb_chunk-1):
		# COMPUTE NUMBER OF BLOCKS AND NZ
		current_nz, last=compute_nblocks[option](m, n, nz_chunk, row, col, arrow, current)
		# OUTPUT CURRENT CHUNK
		f.write(str(i)+"\t"+str(current+1)+"\t"+str(last)+"\t"+str(int(nz_chunk-current_nz))+"\n");

		total+=nz_chunk-current_nz
		nblocks=last-current
		total_blocks+=nblocks

		#PREPARE NEXT CHUNK
		current=last
		nz-=nz_chunk-current_nz
#		# LOG
#		print "current nz: " + str(nz_chunk-current_nz)
#		print "current nblocks: " + str(nblocks)
#		print ""
	# OUTPUT REST IN LAST CHUNK
	last = output_last_blocks_chunk[option](m, n)
	f.write(str(nb_chunk-1)+"\t"+str(current+1)+"\t"+str(last)+"\t"+str(int(nz))+"\n");

	total+=nz
	nblocks=last-current
	total_blocks+=nblocks

#	# LOG
#	print "current nz: " + str(nz)
#	print "current nblocks: " + str(nblocks)
#	print ""
	# TOTAL LOG
	print "total nz outputted: " + str(total)
	print "total blocks outputted: " + str(total_blocks)

	f.close()


#MAIN
f=sys.argv[1]
nb_chunk=int(sys.argv[2])
output=sys.argv[3]
#header=(sys.argv[4]=="true")

print "Getting matrix values..."
m,n,nz,row,col,arrow = count_non_zeros(f)
print "m: " + str(m)
print "n: " + str(n)
print "nz: " + str(nz)
print "length row: " + str(len(row))
print "length col: " + str(len(col))
#print row
#print col
print "nb_chunk: " + str(nb_chunk)
print "nz_chunk: " + str(math.floor(nz/nb_chunk))
print ""
print "Computing chunks..."
for i in range(0, 3):
	print "Option " + str(i)
	output_chunks_per_block(m,n,nz,math.floor(nz/nb_chunk),row,col,arrow,i,output+"_"+str(i))
