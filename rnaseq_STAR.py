#bowtie2_cutnrun.py
from datetime import datetime
import os
import glob
import shutil
import sys
from subprocess import check_output

startTime = datetime.now()

# Change directory within container
# Change directory within container
nav1 = os.path.isdir('/data/trim_out')
nav2 = os.path.isdir('/data/fastq')
if nav1 == True:
	os.chdir('/data/trim_out')
	print 'Current working directory is:' +  os.getcwd()
	print '\n'
elif nav2 == True:
	os.chdir('/data/fastq')
	print 'Current working directory is:' +  os.getcwd()
	print '\n'
else:
	os.chdir('/data')
	print 'Current working directory is:' +  os.getcwd()
	print '\n'

#prepare input files
input_files = sorted(glob.glob('*val*'))
input_files_R1 = []
input_files_R2 = []

if len(input_files) == 0:
	print 'No quality-controlled input files from trim_galore, checking input folder for fastq output...'
	print '\n'
	input_files = sorted(glob.glob('*.fastq'))
	if len(input_files) == 0:
		print 'No uncompressed fastq files detected, looking for fastq.gz...'
		print '\n'
		input_files = sorted(glob.glob('*.fastq.gz'))
		if len(input_files) == 0:
			print 'No valid input files detected. Exiting...'	
			sys.exit()
		else:
			print 'Decompressing .gz files'
			#create fastq filename from .fastq.gz
			temp_list = [f.replace('.gz', '') for f in input_files]
			#run system command to decompress file with zcat
			for i in range(len(input_files)):
				temp_str= 'zcat ' + input_files[i] + ' > ' + temp_list[i]
				check_output(temp_str, shell=True)
			input_files = temp_list

for item in input_files:
	if '_R1_' in item and '_R2_' in item:
		print "Input file: " + item + " contains both strings 'R1' and 'R2'. Not including..."
		sys.exit()
	elif '_R1_' in item and '_R2_' not in item:
		input_files_R1.append(item)
	elif '_R1_' not in item and '_R2_' in item:
		input_files_R2.append(item)
	else:
		print "Input file: " + item + "does not contain string 'R1' or 'R2'. Not including..."

if len(input_files_R1) != len(input_files_R2):
	print 'Unequal numbers of files assigned as R1 and R2. Check naming convention. Exiting...'
	sys.exit()

if (len(input_files_R1) + len(input_files_R2)) != len(input_files):
	print 'Not all of input files assigned as R1 or R2. Exiting...'
	sys.exit()



print 'Files assigned as R1:'
print '\n'.join(input_files_R1)
print '\n'

print 'Files assigned as R2:'
print '\n'.join(input_files_R2)
print '\n'

#make sam output names
STAR_prefixes = []

for i in range(len(input_files_R1)):
	STAR_prefix = input_files_R1[i].split('_R1')[0] + '_'
	STAR_prefixes.append(STAR_prefix)

print 'Output STAR prefixes:'
print '\n'.join(STAR_prefixes)
print '\n'


sam_names = []

for i in range(len(input_files_R1)):
	sam_name = input_files_R1[i].split('_R1')[0] + '_Aligned.out.sam'
	sam_names.append(sam_name)

STAR_logs = []

for i in range(len(input_files_R1)):
	STAR_log = input_files_R1[i].split('_R1')[0] + '_Log.out'
	STAR_logs.append(STAR_log)

for i in range(len(input_files_R1)):
	STAR_log = input_files_R1[i].split('_R1')[0] + '_Log.final.out'
	STAR_logs.append(STAR_log)

for i in range(len(input_files_R1)):
	STAR_log = input_files_R1[i].split('_R1')[0] + '_Log.progress.out'
	STAR_logs.append(STAR_log)

gene_counts = []

for i in range(len(input_files_R1)):
	gene_count = input_files_R1[i].split('_R1')[0] + '_ReadsPerGene.out.tab'
	gene_counts.append(gene_count)


print 'Output sam filenames:'
print '\n'.join(sam_names)
print '\n'

#define bowtie2 index for hg19
print 'Data will be aligned to hg19'

hg19_index = '/genomes/STAR/'
print 'STAR index for data files found at:'
print hg19_index
print '\n'

#run bowtie for hg19
print 'Running STAR alignment for hg19 followed by read counts'
print '\n'
for i in range(len(input_files_R1)):
	print 'count = ' +str(i)
	print '\n'
	#create string for system command
	temp_str = 'STAR --runThreadN 16 --genomeDir ' + hg19_index + ' --readFilesIn ' + input_files_R1[i] + ' ' + input_files_R2[i] + ' --outFileNamePrefix ' + STAR_prefixes[i] + ' --quantMode GeneCounts --twopassMode Basic --outSAMunmapped Within'

	print temp_str

	check_output(temp_str, shell=True)
	print '\n'


#make new directory for output
os.mkdir('/data/sams')
os.mkdir('/data/sams/logs')
os.mkdir('/data/sams/counts')
print 'Current working directory is:' +  os.getcwd()
print '\n'

#copy files to output folder
output_dir = '/data/sams/'
print 'Moving sam files to output folder'
print '\n'
for i in range(len(sam_names)):
	shutil.move(sam_names[i], output_dir)

output_dir = '/data/sams/logs'
print 'Moving sam files to output folder'
print '\n'
for i in range(len(STAR_logs)):
	shutil.move(STAR_logs[i], output_dir)

output_dir = '/data/sams/counts'
print 'Moving sam files to output folder'
print '\n'
for i in range(len(gene_counts)):
	shutil.move(gene_counts[i], output_dir)


print 'Alignment Runtime (hh:mm:ss): ' + str(datetime.now() - startTime)
print '\n'

###SAM conversion to bam, bedgraph, and BigWig

os.mkdir('/data/bams')
os.chdir('/data/sams')
print 'Current working directory is:' +  os.getcwd()
print '\n'

import pybedtools
from pybedtools import BedTool
import pandas as pd
from pybedtools.helpers import chromsizes
from pybedtools.contrib.bigwig import bedgraph_to_bigwig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

startTime = datetime.now()

#######################################################
## spike in or average normalization with file generation#########
#######################################################

print 'Converting sams to bams, bedgraph, and BigWig without spike-in normalization.'
print '\n'
datafiles = sorted(glob.glob('*.sam'))

print '\n'
print 'Data files loaded:'
print '\n'.join(datafiles)
print '\n'


##convert to bam format
print 'Converting to bam format'
print '\n'

bam_names = []

for i in range(len(sam_names)):
	bam_name = sam_names[i].split('_Aligned')[0] + '.bam'
	bam_names.append(bam_name)

print '\n'
print '\n'.join(bam_names)
print '\n'
print '\n'
print 'SAM to BAM'
print '\n'

bam_string = []

for i in range(len(bam_names)):
        bam_string.append('samtools view -b -S ' + datafiles[i] + ' > /data/bams/' + bam_names[i])

for item in bam_string:
        check_output(item, shell = True)

datafiles = bam_names

##Sort and index
print '\n'
print 'Sorting bams'
print '\n'

sorted_bam_names = []
for i in range(len(bam_names)):
	sorted_bam_name = 'sorted.' + bam_names[i]
	sorted_bam_names.append(sorted_bam_name)

bam_names_sorted = [f.replace('.bam', '') for f in sorted_bam_names]

sort_string = []

for i in range(len(bam_names)):
		sort_string.append('samtools sort /data/bams/' + bam_names[i] + ' /data/bams/' + bam_names_sorted[i])

for item in sort_string:
		check_output(item, shell = True)

print '\n'
print 'Bam files to index:'
print '\n'.join(bam_names_sorted)
print '\n'


print 'Indexing bams'
print '\n'             
os.chdir('/data/bams')
print 'Current working directory is:' +  os.getcwd()

index_string = []

for i in range(len(bam_names_sorted)):
        index_string.append('samtools index ' + sorted_bam_names[i])

for item in index_string:
        check_output(item, shell = True)

sorted_files = sorted(glob.glob('sorted*'))
os.mkdir('/data/bams/sorted_bams')
output_dir = '/data/bams/sorted_bams'
for i in range(len(sorted_files)):
	shutil.move(sorted_files[i], output_dir)

	
##generate bigwig files from bam files
print 'Generating bigwig files'
print '\n'

os.chdir('/data/bams/sorted_bams')

print 'Current working directory is:' +  os.getcwd()
print '\n'

#generate bigwig file names
bw_names = [f.replace('bam', 'bw') for f in datafiles]


bigwig_string = []

for i in range(len(sorted_bam_names)):
		bigwig_string.append('bamCoverage -b '+ sorted_bam_names[i] + ' -o /data/bams/sorted_bams/' + bw_names[i])

for item in bigwig_string:
		check_output(item, shell = True)


print 'Finished generating bigwig files:'
print '\n'
print 'whole insert bigwig files:' + '\n' + '\n'.join(bw_names)
print '\n'

os.mkdir('/data/bigwigs')

print 'Moving bigwigs to output folder'
print '\n'
output_dir0 = '/data/bigwigs'
for i in range(len(bw_names)):
	shutil.move(bw_names[i], output_dir0)


print 'Finished'
print '\n'
print 'Runtime (hh:mm:ss): ' + str(datetime.now() - startTime)
	
