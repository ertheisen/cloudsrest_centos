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
nav2 = os.path.isdir('/data/trim_out')
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
	print 'No quality-controlled input files from trim_galore, checking input folder for fastqc output...'
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
sam_names = []

for i in range(len(input_files_R1)):
	sam_name = input_files_R1[i].split('_L0')[0] + '.sam'
	sam_names.append(sam_name)

print 'Output sam filenames for data species:'
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
	temp_str = 'STAR --runThreadN 16 --genomeDir ' + hg19_index + '--readFilesIn ' + input_files_R1[i] + ' ' + input_files_R2[i] + ' --outFileNamePrefix ' + sam_names[i] + ' --quantMode GeneCounts --twopassMode Basic --outSAMunmapped Within'

	print temp_str

	check_output(temp_str, shell=True)
	print '\n'


#make new directory for output
os.mkdir('/data/sams')
print 'Current working directory is:' +  os.getcwd()
print '\n'

#copy files to output folder
output_dir = '/data/sams/'
print 'Moving sam files to output folder'
print '\n'
for i in range(len(sam_names)):
	shutil.move(sam_names[i], output_dir)

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

bam_names = [f.replace('sam', 'bam') for f in datafiles]

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

	
##generate bed files from bam files
print 'Generating bed files representing whole insert from paired end reads in the data files'
print '\n'

print 'Current working directory is:' +  os.getcwd()
print '\n'

#generate bed file names
bed_names = [f.replace('bam', 'bed') for f in datafiles]

#generate file names for length analysis
lengths_names = [f.replace('bam', 'lengths') for f in datafiles]

#generate bed files with bam_to_bed tool (makes bed12 format)
for i in range(len(datafiles)):
	temp_bed = BedTool(datafiles[i]).bam_to_bed(bedpe=True).to_dataframe()

	#need to strip out start and end position of whole insert (bed12 is both reads)
	#column names actually represent <chrom>, <start of insert>, <end of insert>
	temp_bed_stripped = temp_bed.iloc[:,[0,1,5]].sort_values(by = ['chrom', 'start', 'strand'])

	#calculate insert size as column 4 and save file with bed_name
	temp_bed_stripped['length'] = temp_bed_stripped['strand'] - temp_bed_stripped['start']

	temp_bed_stripped.to_csv(bed_names[i], sep="\t", header = False, index = False)

	#analyze lengths of inserts
	temp_lengths = temp_bed_stripped.groupby(by=['length'])['length'].count()

	temp_lengths.to_csv(lengths_names[i], sep="\t", header = [bed_names[i]], index = True, index_label='length')


print 'Finished generating bed files:'
print '\n'
print 'whole insert bed files:' + '\n' + '\n'.join(bed_names)
print '\n'

#generate normalized bedgraphs
print 'Current working directory is:' +  os.getcwd()
print '\n'

#generate bedgraph names
bg_names = [f.replace('bed', 'bg') for f in bed_names]

print 'Generating average normalized bedgraphs for all the bed files'
print '\n'
#count total number of reads in each bed file (before size selection)
read_count = []
for item in bed_names:
	read_count.append(BedTool(item).count())

print read_count

#calculate genome size
genome_file = chromsizes('hg19')
DF = pd.DataFrame.from_dict(genome_file, orient='index')
genome_size = DF[1].sum()
	
scaling_factor = []
for item in read_count:
	scaling_factor.append(float(genome_size) / read_count[i])

for i in range(len(bed_names)):
	count = str(read_count[i])	
	print '\n'
	print 'for ' + bed_names[i] + ' read count:'
	print count
	print '\n'
	print 'scaling factor for ' + bed_names[i] + ' is:'
	print scaling_factor[i]

#run bedtools genomecov to generate bedgraph files
for i in range(len(bg_names)):
	BedTool(bed_names[i]).genome_coverage(bg = True, genome = 'hg19', scale = scaling_factor[i]).moveto(bg_names[i])

print 'Finished generating bedgraph files:'
print '\n'
print 'whole insert bedgraph files:' + '\n' + '\n'.join(bg_names)
print '\n'

##make bigwig files
print 'Current working directory is:' +  os.getcwd()
print '\n'

print 'Generating big_wig files for all the bedgraphs'
print '\n'

#generate bigwig names
bw_names = [f.replace('bg', 'bw') for f in bg_names]

#run bedgraph_to_bigwig tool
for i in range(len(bg_names)):
	bedgraph_to_bigwig(BedTool(bg_names[i]), 'hg19', bw_names[i])

print 'Finished generating bigwig files:'
print '\n'
print 'whole insert bigwig files:' + '\n' + '\n'.join(bw_names)
print '\n'

os.mkdir('/data/beds')
os.mkdir('/data/beds/lengths')
os.mkdir('/data/bedgraphs')
os.mkdir('/data/bigwigs')

output_dir0 = '/data/beds'
output_dir1 = '/data/beds/lengths'
for i in range(len(bed_names)):
	shutil.move(bed_names[i], output_dir0)
for i in range(len(lengths_names)):
	shutil.move(lengths_names[i], output_dir1)

print 'Moving bedgraphs to output folder'
print '\n'
output_dir4 = '/data/bedgraphs'
for i in range(len(bg_names)):
	shutil.move(bg_names[i], output_dir4)

print 'Moving bigwigs to output folder'
print '\n'
output_dir7 = '/data/bigwigs'
for i in range(len(bw_names)):
	shutil.move(bw_names[i], output_dir7)


print 'Finished'
print '\n'
print 'Runtime (hh:mm:ss): ' + str(datetime.now() - startTime)
	
