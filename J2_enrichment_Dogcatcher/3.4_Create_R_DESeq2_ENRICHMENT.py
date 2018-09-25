#!/usr/bin/env python
import argparse
import csv
import numpy as np
import os
import shlex
import shutil
import subprocess
import sys
import pandas as pd
import string
import glob

# Running in parallel
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Part II: Calculating run in/on for sense and antisense strands.')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--input_R_template_file', action= 'store', metavar='input_R_template_file', help='input_R_template_file. Example "R_template/R_template_new.R" ')

    parser.add_argument('--output_prefix', action= 'store', metavar='output_prefix', default = "Runion_out/", help='Output prefix. KEEP THIS THE SAME FROM PART II: Default: Runion_out/')
    parser.add_argument('--input_prefix', action= 'store', metavar='input_prefix', default = "Runion_out/", help='input prefix. MAKE THIS PART II OUTPUT_PREFIX : Default: Runion_out/')
    parser.add_argument('--cpus', action= 'store', dest='cpus', metavar='cpus', default= 1, type=int, help='Enter available cpus per node.  The more cpus the faster Runion performs. Default: "1"')


    args = parser.parse_args()


    input_R_template_file = args.input_R_template_file
    output_prefix = args.output_prefix
    input_prefix = args.input_prefix
    cpus = args.cpus


    #########################################
    #Create initial R subread and DSeq2 script

    def create_R_script():
        with open(input_R_template_file, "r") as infile, open(output_prefix + "/" + "DESeq2_initial.R", "w") as outfile:
            for line in infile:
                line = line.replace("outputprefix", output_prefix)
                line = line.replace("COL_DATA", "col_data.txt")
                line = line.replace("cpus", str(cpus))
                outfile.write(line)
    create_R_script()
