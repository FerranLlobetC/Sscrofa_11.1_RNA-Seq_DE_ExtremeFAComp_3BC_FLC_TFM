#!/home/fllobet/miniconda3/bin/python3

# Module import
import xlrd
import sys
import re
# Loading of the file
file = xlrd.open_workbook("/home/fllobet/Documents/FGJ_Fluidigm_FLC/Primers/DEFINITIVE/IDT_PRIMERS_ORDER.xlsx")

# Reading of the indicated sheet
file_sheet = file.sheet_by_index(0)

# The elements of the first columns are the primer names
primers = file_sheet.col(0)[1:]  # The first (0) element is the column header which is not any of the primers

""" Removal of the 'text' before each element + FW and RV tails """
for i in range(0, len(primers)):
    primers[i] = primers[i].value # Remvola of the 'text'
    primers[i] = re.sub("_[A-Z]+", "", str(primers[i]))  # Removal of the '_FW' and '_RV' tails

""" Removal of duplicated elements (FW and RV pairs) """
primers = set(primers)   # Conversion into a set: which only stores repeated elements once
primers = list(primers)  # Re-conversion into a list
len(primers) # CHECK
