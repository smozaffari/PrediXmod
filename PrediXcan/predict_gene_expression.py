#!/usr/bin/env python

from collections import defaultdict
import datetime
import gzip
import numpy as np
import os
import sqlite3
import sys
from cStringIO import *

DATA_DIR = 'data'
GENE_LIST = os.path.join(DATA_DIR, 'genelist.txt')
DOSAGE_DIR = os.path.join(DATA_DIR, 'dosages')
BETA_FILE = os.path.join(DATA_DIR, 'betas.db')
OUTPUT_FILE = 'output_file.txt'

# We start from scratch
if os.path.exists(OUTPUT_FILE):
    os.unlink(OUTPUT_FILE)

class DBClass:
    def __init__(self):
        self.conn = sqlite3.connect(BETA_FILE)
    
    def __call__(self, sql, args=None):
        c = self.conn.cursor()
        if args:
            for ret in c.execute(sql, args):
                yield ret
        else:
            for ret in c.execute(sql):
                yield ret

                    
betas = DBClass()

def all_genes():
    genes = defaultdict(list) # chromosome -> [genes]
    for line in open(GENE_LIST):
        chromosome, gene = line.strip().split()
        genes[chromosome].append(gene)
    for chromosome in sorted(genes):
        print datetime.datetime.now(), "Processing %s"%chromosome
        for gene in sorted(genes[chromosome]):
           yield gene, chromosome


def get_gene_data(gene):
    raids, weights, eff_alleles = [], [], []
    for (raid, weight, eff_allele) in betas('SELECT raid, weight, eff_allele FROM weights WHERE gene=? ORDER BY raid', (gene,)):
        raids.append(raid)
        weights.append(weight)
        eff_alleles.append(eff_allele)
    return raids, weights, eff_alleles        
        
    

class Chromosome:
    def __init__(self, label):
        self.refalele = {}
        self.dosagerow = {}
        for line in gzip.open(os.path.join(DOSAGE_DIR, '%s.dosage.gz'%label)):
            arr = line.strip().split()
            _rsid = arr[1]
            _refalele = arr[3]
            _dosagerow = np.array(map(float, arr[6:]))
            self.refalele[_rsid] = _refalele
            self.dosagerow[_rsid] = _dosagerow
        
class ChromosomeManager:
    def __init__(self):
        self.current_chromosome_label = None
        self.current_chromosome_object = None

    def __call__(self, chromosome_label):
        if not self.current_chromosome_label == chromosome_label:
            self.current_chromosome_label = chromosome_label
            self.current_chromosome_object = Chromosome(chromosome_label)
        return self.current_chromosome_object            
get_chromosome = ChromosomeManager()

def get_dosage_data(rsids, eff_alleles, chromosome):
    D = None
    chromosome = get_chromosome(chromosome)
    for idx, rsid in enumerate(rsids):
        if not rsid in chromosome.refalele:
            continue
        refalele = chromosome.refalele[rsid]
        multiplier = 1 if refalele == eff_alleles[idx] else -1
        dosagerow = multiplier * chromosome.dosagerow[rsid]
        if D is None: # We didn't know the number of samples until now
            D = np.zeros((len(rsids), len(dosagerow)))
        D[idx,] = dosagerow
    return D

for gene, chromosome in all_genes():
    rsids, weights, eff_alleles = get_gene_data(gene)
    D = get_dosage_data(rsids, eff_alleles, chromosome)
    l = len(weights)
    weights = np.array(weights)
    weights.shape = (1,len(weights))
    Tg = np.dot(weights, D)

    with open(OUTPUT_FILE, 'a+') as outfile:
        outfile.write("%s %s\n"%(gene, ' '.join(map(str, Tg.tolist()[0]))))


