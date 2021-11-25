#!/usr/bin/env python3

# 11/2021
# alexandrina.bodrug@chu-angers.fr

##############
#dependencies#
##############
#the script assumes you can call from command line: bcftools seekin ls file
#beware of bcftools version, --missing-to-ref is an option not present in older versions
#the script assumes the --vcf-dir contains gzipped vcfs and their tbi
#GraphViz libraries should be installed both for python, with pip or conda, AND on your system (dnf, apt-get, etc)

#modules
import sys
import os
if sys.version_info<(3,6,0): exit("\n[SYS ERROR] You need python 3.6 or later to run this script\n\n")
import argparse as ap
import os
import subprocess
import time
from graphviz import Digraph
import vcf
import gzip
import re
from PIL import Image
from fpdf import FPDF

#######################
#timestamp and version#
#######################
start_time = time.time()
VERSION = "Version 2.0"

#helper

if len(sys.argv) < 2:
	sys.argv.append("-h")

if sys.argv[1] in ["-h", "-help", "--help", "getopt", "usage"]:

	sys.exit('''
        {0}
        This script computes the kinship relationship between several vcfs with the programs
        seekin and KING. It outputs a picture of a graph representing the detected relationships.

        The script assumes your files are bgzipped and indexed with tabix. If it is not the case, run
        something like >for file in *vcf; do bgzip $file ; tabix -p vcf file.gz; done

        Example of command: python3 kinship.py -f vcfdir1 vcfdir2 -o outdir -p kinpic -t 1 -k 0.01 -v


        -f    --vcf-dir           directory(ies) containing vcf files

        -o    --output-dir        directory that will contain the output files and picture

        -p    --output-picture    prefix of picture (don't specify extension)

        -t    --threads           threads to use when running seekin or KING

        -k    --kinship-threshold threshold under which seekin won't consider there is a relationship

        -v    --verbose           run with more verbosity

        -s    --sex               compute the sex

        -h    --help              display this message

    '''.format(VERSION))

###############
#option parser#
###############
p = ap.ArgumentParser()
p.add_argument("-f", "--vcf-dir", type=str, required=True, action='append', nargs='+') #dir containing vcfs to test, multiple dirs can be passed
p.add_argument("-o", "--output-dir", type=str, required=False, default="outdir") #dir that will contain outputs
p.add_argument("-p", "--output-picture", type=str, required=False, default="kinpic") #name of the picture produced
p.add_argument("-t", "--threads", type=int, required=False, default=1) #name of the picture produced
p.add_argument("-k", "--kinship-threshold", type=float, required=False, default=0.01) #kinship thresholds
p.add_argument("-v", "--verbose", required=False, action='store_true', default=False) #verbosity
p.add_argument("-s", "--sex", required=False, action='store_true', default=False) #if one should compute the sex
p.add_argument("-r", "--rerun", required=False, action='store_true', default=False) #if one wants to rerun to make the graphs but not the calculus

args = p.parse_args()

###########
#functions#
###########
def verify_format(file_list):
    """
        Description: checks that the list of files provided are indeed vcfs
        Input / Output: list of files / return none
    """
    for file in file_list:
        cmd = "bcftools view -h " + file
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.communicate() # this line with make it so .returncode actually returns the return code, and not 'None'
        if args.verbose:
            print(file,'is indeed a vcf file.', process.returncode)
        if process.returncode != 0:
            print("ERROR: file is not vcf ", file, ". Try launching this command: >", cmd); exit()
    return 0
#
def run_kinship(file_list):
    """
        Description: merges all vcf files to prepare input for kinship
        Input / Output: list of vcf files (gzipped) / no return, just create kinship's output in the output_dir
    """
    emplacementScript = "/".join(os.path.realpath(__file__).split("/")[:-1])
    # if it's not a rerun, create output_dir and merged vcf, then run seekin
    if args.rerun == False:
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)

        cmd = "bcftools merge --missing-to-ref --threads "+ str(args.threads)+ " " + " ".join(file_list) +" -o "+ args.output_dir +"/kinship_merge.vcf"
        if args.verbose :
            sys.stderr.write("Launched bcftools merge: "+ cmd +"\n\n")
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        if process.returncode!=0:
            print("ERROR: in bcftools merge\n Try launching this to find out why: >" + cmd) ; exit()

        #run seekin
        cmd =emplacementScript+"/seekin kinship -i "+ args.output_dir +"/kinship_merge.vcf -r 0.3 -m 0.05 -d GT -p hom -t "+ str(args.threads) +" -w 1 -o "+ args.output_dir +"/kinship"
        if args.verbose :
            sys.stderr.write("Launched seekin kinship: "+ cmd +"\n\n")
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        process.wait()
        if "Finished!" not in str(process.communicate()):
            print("ERROR: in seekin\nTry launching this to find out why: >" + cmd) ; exit()
    return 0
#
def compute_sex(file_list):
    samplesex_hash = {}
    for file in file_list:
        reader = vcf.Reader(filename=file)
        sample = reader.samples
        samplename = sample[0]
        if args.sex:
            lAF_autosomes = []
            lAF_X = []
            #lAF_Y = []
            for record in reader:
               if re.match("chr[0-9]", record.CHROM) :
                   lAF_autosomes.extend(record.INFO['AF'])  #this match hits chr1 to chr22
               if record.CHROM == "chr"+"X" :
                   lAF_X.extend(record.INFO['AF'])
               #if record.CHROM == "chr"+"Y" :
                   #lAF_Y.extend(record.INFO['AF'])

            avg_AF_autosomes = sum(lAF_autosomes) / len(lAF_autosomes)
            avg_AF_X = sum(lAF_X) / len(lAF_X)

            interval = 0.10 # fake standard deviation
            if avg_AF_X >= (avg_AF_autosomes - interval) and avg_AF_X <= (avg_AF_autosomes + interval):
                sex='female'
                samplesex_hash[samplename] = sex
            else:
                sex='male'
                samplesex_hash[samplename] = sex
            if args.verbose :
                print(file,"is",sex)
        else:
            sex='unknown' #unknown sex
            samplesex_hash[samplename] = sex
    return samplesex_hash
#
def make_graph(kinfile, samplesex_hash):
    """
        Description: construct kinship graph based on seekin output. All edges included
        Input / Output: seekin output / graph
    """
    #extract nodes and edges
    file1 = open(kinfile, 'r')
    Lines = file1.readlines()
    nodes = []
    edges = []
    for line in Lines:
        line = line.strip()
        if ("ind1	ind2	nsnp	kinship" not in line):
            sample1, sample2, nsnp, kinscore = line.split()
            nodes.append(sample1)
            nodes.append(sample2)
            if (float(kinscore) > args.kinship_threshold):
                edges.append({'sample1':sample1, 'sample2':sample2, 'kinscore':kinscore})
    #test zone
    #print(samplesex_hash)
    #print(edges)
    #print(nodes)

    #make graĥviz graph
    gsk = Digraph(name='none', filename='kinship',comment='none', engine='circo', format='png')
    for node in nodes:
        if samplesex_hash[node] == 'male':
            gsk.node(node, shape='polygon', fillcolor='lightblue1', color='black', style='filled')
        elif samplesex_hash[node] == 'female':
            gsk.node(node, shape='ellipse', fillcolor='lightsalmon', color='black', style='filled')
        else:
            gsk.node(node, shape='diamond', fillcolor='ivory2', color='black', style='filled')
    for edge in edges:
        sample1 = edge['sample1']
        sample2 = edge['sample2']
        kinscore = edge['kinscore']
        gsk.edge(sample1, sample2, dir='none', color='black', penwidth='2', label=kinscore, fontname='Helvetica', fontsize='12', fontcolor='limegreen')

    #render picture
    ugsk = gsk.unflatten(stagger=3)
    ugsk.render(args.output_dir+'/'+args.output_picture+'_seekin', view=False)
    return 0


######
#main#
######
if __name__ == "__main__":
    # the command run by the user is printed to STDERR
    command = " ".join(sys.argv)
    if args.verbose :
        sys.stderr.write("Launched kinship2 command: {0}\n\n".format(command))

    # harvest file list
    file_list = []
    for dir in args.vcf_dir[0]:
        for file in os.listdir(dir):
            if "vcf" in file and "tbi" not in file:
                file_list.append(dir+"/"+file)
    # check if vcfs
    verify_format(file_list)
    #compute sex
    samplesex_hash = compute_sex(file_list)
    # compute kinship
    run_kinship(file_list)
    # compute graph
    make_graph(args.output_dir+"/kinship.kin", samplesex_hash)
