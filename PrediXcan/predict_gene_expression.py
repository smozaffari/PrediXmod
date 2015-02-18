#!/usr/bin/env python
from collections import defaultdict
import datetime
import gzip
import numpy as np
import os
import sqlite3
import sys

DATA_DIR = 'data'
GENE_LIST = os.path.join(DATA_DIR, 'genelist.txt')
DOSAGE_DIR = os.path.join(DATA_DIR, 'dosages')
BETA_FILE = os.path.join(DATA_DIR, 'betas.db')
OUTPUT_FILE = 'output_file.txt'


def get_all_dosages():
    print "ONLY PROCESSING CHR1 FOR TESTING PURPOSES"
    for chrfile in list(sorted(os.listdir(DOSAGE_DIR)))[:1]:
        print datetime.datetime.now(), "Processing %s"%chrfile
        for line in gzip.open(os.path.join(DOSAGE_DIR, chrfile)):
            arr = line.strip().split()
            rsid = arr[1]
            refallele = arr[3]
            dosage_row = np.array(map(float, arr[6:]))
            yield rsid, refallele, dosage_row

class GetApplicationsOf:
    def __init__(self):
        self.conn = sqlite3.connect(BETA_FILE)

    def query(self, sql, args=None):
        c = self.conn.cursor()
        if args:
            for ret in c.execute(sql, args):
                yield ret
        else:
            for ret in c.execute(sql):
                yield ret

    def __call__(self, rsid):
        for tup in self.query("SELECT gene, weight, eff_allele FROM weights WHERE raid=?", (rsid,)):
            yield tup
get_applications_of = GetApplicationsOf()

class TranscriptionMatrix:
    def __init__(self):
        self.D = None

    def update(self, gene, weight, ref_allele, allele, dosage_row):
        if self.D is None:
            self.gene_list = list(sorted([line.strip().split()[-1] for line in open(GENE_LIST)]))
            self.gene_index = { gene:k for (k, gene) in enumerate(self.gene_list) }
            self.D = np.zeros((len(self.gene_list), len(dosage_row))) # Genes x Cases
        if gene in self.gene_index:            
            multiplier = 1 if ref_allele == allele else -1
            self.D[self.gene_index[gene],] += dosage_row * weight * multiplier # Update all cases for that gene

    def save(self):
        with open(OUTPUT_FILE, 'w+') as outfile:
            outfile.write('\t'.join(self.gene_list) + '\n') # Nb. this lists the names of rows, not of columns
            outfile.write(str(self.D))

transcription_matrix = TranscriptionMatrix()
for rsid, allele, dosage_row in get_all_dosages():
    for gene, weight, ref_allele in get_applications_of(rsid):
        transcription_matrix.update(gene, weight, ref_allele, allele, dosage_row)
transcription_matrix.save()        
