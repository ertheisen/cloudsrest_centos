#bowtie2_cutnrun.py
from datetime import datetime
import os
import glob
import shutil
import sys
from subprocess import check_output
import pybedtools
from pybedtools import BedTool
import pandas as pd
from pybedtools.helpers import chromsizes
from pybedtools.contrib.bigwig import bedgraph_to_bigwig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


startTime = datetime.now()
#######################################
# Change directory within container
os.chdir('/data')
print 'Current working directory is:' +  os.getcwd()
#######################################
bam_files = sorted(glob.glob('*bam'))
if bam_files == 0:
	# Change directory within container
	os.chdir('/data/bams')
	print 'Current working directory is:' +  os.getcwd()
	bam_files = sorted(glob.glob('*bam'))
	if bam_files == 0:
		print 'No bam files detected. Exiting...'
		sys.exit()
	else:
		work_dir = 2
else:
	work_dir = 1

bai_files = sorted(glob.glob('*bai'))
if bai_files == 0:
	print 'Bam files will be indexed. Sorting and indexing...'
	for i in range(len(bam_files)):
		print '\n'
		#create string for system command to sort
		temp_str = 'samtools sort ' + bam_files[i] + ' sorted.' + bam_files[i].split('.')[0]
		print temp_str

		check_output(temp_str, shell=True)
		print '\n'

	sorted_files = sorted(glob.glob('sorted*'))
	os.mkdir('/data/sorted_bams')
	print 'Moving sorted bams to new directory'
	output_dir = '/data/sorted_bams'
	for i in range(len(sorted_files)):
		shutil.move(sorted_files[i], output_dir)
	os.chdir('/data/sorted_bams')
	for i in range(len(sorted_files)):
		print '\n'	
		#create string for system command to index
		temp_str = 'samtools index /data/sorted_bams' + sorted_files[i]
		print temp_str

		check_output(temp_str, shell=True)
		print '\n'

##############################################

###BAM conversion to BED, bedgraph, and BigWig

print 'Converting bams to beds, bedgraph, and BigWig without spike-in normalization.'
print '\n'

datafiles = bam_files
##generate bed files from bam files
print 'Generating bed files representing whole insert from paired end reads in the data files'
print '\n'

if work_dir == 1:
	os.chdir('/data')

if work_dir == 2:
	os.chdir('/data/bams')

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
output_dir2 = '/data/bedgraphs'
for i in range(len(bg_names)):
	shutil.move(bg_names[i], output_dir2)

print 'Moving bigwigs to output folder'
print '\n'
output_dir3 = '/data/bigwigs'
for i in range(len(bw_names)):
	shutil.move(bw_names[i], output_dir3)

print 'Finished'
print '\n'
print 'Runtime (hh:mm:ss): ' + str(datetime.now() - startTime)

#run plotFingerprint

startTime = datetime.now()

os.chdir('/data/sorted_bams')
print 'Current working directory is:' +  os.getcwd()

print 'Running plotFingerprint for all bam files in folder'
print '\n'
temp_str = 'plotFingerprint -b *.bam --smartLabels --minFragmentLength 20 -T "Mapping Fingerprints" --plotFile mapqc_fingerprint.png --outQualityMetrics mapqc.tab'

print temp_str

check_output(temp_str, shell=True)
print '\n'

print 'Finished'
print '\n'
print 'Runtime mapping QC (hh:mm:ss): ' + str(datetime.now() - startTime)
	
