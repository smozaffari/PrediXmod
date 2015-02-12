#!/usr/bin/env python

SOURCE_DIR = 'GTEx-model-results'
TARGET_DIR = 'generated_dbs'

import os
import sys
import sqlite3


def source_files(source_dir=SOURCE_DIR):
    "List all relevant source files"
    for x in sorted(os.listdir(source_dir)):
        if x.endswith('.allBetas.txt'):
            yield os.path.join(source_dir, x)


def data_rows_in(source_file):
    "Iterate over data rows in the source file, labeling fields and converting formats as required."
    def upconvert(x):
        for f in (int, float):
            try:
                return f(x)
            except ValueError:
                pass
        return x

    for k, line in enumerate(open(source_file)):
        if k == 0:
            header = line.strip().split()
        else:
            yield dict(zip(header, map(upconvert, line.strip().split())))

class DB:
    "This encapsulates a single SQLite DB (for a given source file and alpha)."
    def __init__(self, source_file, alpha, target_dir=TARGET_DIR):
        tissue_name = os.path.basename(source_file).split('.')[0]
        db_filename = os.path.join(target_dir, '%s_%s.db'%(tissue_name, alpha))
        if not os.path.exists(target_dir):
            os.mkdir(target_dir)

        if os.path.exists(db_filename):
            os.unlink(db_filename)

        self.connection = sqlite3.connect(db_filename)
        
        self("CREATE TABLE weights (raid TEXT, gene TEXT, weight DOUBLE, eff_allele CHARACTER, pval DOUBLE, N INTEGER, cis INTEGER)")
        self("CREATE INDEX weights_raid ON weights (raid)")
        self("CREATE INDEX weights_raid_gene ON weights (raid, gene)")


    def __call__(self, sql, args=None):
        c = self.connection.cursor()
        if args:
            c.execute(sql, args)
        else:
            c.execute(sql)

    def close(self):
        self.connection.commit()            

    def insert_row(self, row):
        self("INSERT INTO weights VALUES(?, ?, ?, ?, NULL, NULL, NULL)", (row['rsid'], row['gene'], row['beta'], row['ref']))
        

class MetaDB:
    "This handles all the DBs for each source file (tissue type)"
    def __init__(self, source_file):
        self.source_file = source_file
        self.dbs = {} # alpha -> DB object

    def insert_row(self, row):
        alpha = row['alpha']
        if alpha not in self.dbs:
            self.dbs[alpha] = DB(self.source_file, alpha)
        self.dbs[alpha].insert_row(row)

    def close(self):
        for db in self.dbs.values():
            db.close()
        


if __name__ == '__main__':
    for source_file in source_files():
        print "Processing %s..."%source_file
        meta_db = MetaDB(source_file=source_file)
        for row in data_rows_in(source_file):
            meta_db.insert_row(row)
        meta_db.close()
