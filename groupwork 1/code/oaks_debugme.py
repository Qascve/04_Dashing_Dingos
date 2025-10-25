#!/usr/bin/env python3

"""
Filters oak species from CSV data based on genus name matching.

This script reads a CSV file containing taxa data and extracts only the oak
species (genus Quercus) using fuzzy string matching to handle minor typos.

Author: Anaga Ambady (aa6725@ic.ac.uk)
Version: 1.0.0
Date: Oct 2025
"""

__author__ = 'Anaga Ambady (aa6725@ic.ac.uk)'
__version__ = '1.0.0'

import csv
import sys
import doctest
import difflib


def is_an_oak(name):
    """Returns True if name starts with 'quercus' or is close to it (fuzzy matching)
   
    >>> is_an_oak('Fagus sylvatica')
    False
    >>> is_an_oak('Quercus')
    True
    >>> is_an_oak('Pinus')
    False
    >>> is_an_oak('quercus ')
    True
    >>> is_an_oak('quercuss')
    True
    >>> is_an_oak('Quercuss')
    True
    >>> is_an_oak('Quercus robur')
    True
    >>> is_an_oak('quercus ilex')
    True
    """ 
    # Extract genus name (first word) and normalize
    genus = name.strip().split()[0].lower()
    
    # Use fuzzy matching to allow for minor typos
    # A ratio > 0.8 means high similarity to 'quercus'
    similarity = difflib.SequenceMatcher(None, genus, 'quercus').ratio()
    return similarity > 0.8


def main(argv): 
    with open('../data/TestOaksData.csv', 'r') as f, \
         open('../data/JustOaksData.csv', 'w', newline='') as g:
        
        taxa = csv.reader(f) 
        header = next(taxa)  # Skip header row
        
        csvwrite = csv.writer(g)
        csvwrite.writerow(header)  # Write header to output file
        
        oaks_count = 0
        
        for row in taxa:
            if not row:  # Skip empty rows
                continue
                
            genus = row[0].strip()
            print(f"The genus is: {genus}")
            
            if is_an_oak(genus):
                print('FOUND AN OAK!\n')
                csvwrite.writerow(row)
                oaks_count += 1
        
        print(f"\nTotal oaks found: {oaks_count}")
    
    return 0

    
if __name__ == "__main__":
    status = main(sys.argv)
    print("\nRunning doctests...\n")
    doctest.testmod(verbose=True)