# Sketch about Workflow of the Software based on the Flowchart on https://app.diagrams.net/#G1pZdeIWqxWHPq4vcbU_8pv2OpZDH51izH

# potential IMPORTS
import yaml
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import itertools

# 1. Check imports 
# 2. Compute absolute values of lineage number
# 3. Compute absolute values of aa mutation from its lineage respectively 
# 4. Compute frequency from 2. and 3. 
# 5. Handle output based on option args
#   5.1. Save as Matrix 
#       5.1.1 Save as Heatmap
#   5.2 Save frequency "on the fly to perform users tasks"
#       5.1.2 Create Consensus Sequence from user defined sequences


# ---------------------------- IMPORT ----------------------------
# input: yaml file (als Argument)
# output: dictionary d
# Dictonary als Datenstrukur
# enthält Information über markierte Indices und die zu betrachtenen Indices


# -------------------------- INIT ----------------------------
# Für die Datenstruktur
# 1. Um bis unterste Ebene zu iterieren: recursive_items()
# input: dictionary
# output: (key,) value

# ------ GET INDICES
# input: Dictionary d
# output:
#    forward_reverse: liste mit allen indices als (forward, reverse) tupel
#    selected: liste aus dictionaries aus den indices und den jeweilgen Samplenamen
#    index_for_matrix: liste, die sowohl dictonaries als auch tupel enthält
#   diejenige liste, die später für die Matrix eingelesen wird


# -------------------------- OUTPUT ----------------------------
# Editdistance (Levinthstein Distanz - minimale Editdistanz)
# input: sequenz a, b
# 1. Fall: Distanz zwischen zwei Strings
# 2. Fall: Distanz zwischen zwei Tupel
# 3. Fall: Distanz zwischen einem Dictonary und einem Tupel
# output: Distanz zwischen beiden Sequenzen (integer)


# -------------------------- OPTIONS 
# nachdem die matrix erstellt wurde
# input: indexsequenzen als liste
# output: array matrix in n*n, wobei n die anzahl der indices in der liste entsprechen

# -------------------------- Greating Consensus
# Maximale Distanz paarweiser indexsequenzen

