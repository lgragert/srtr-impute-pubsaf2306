
import os
import pandas as pd
import hlagenie


# Initialize `GENIE` object with a version of the IMGT/HLA database
genie = hlagenie.init("3510")

# %%

HLA_B_ALLELES = pd.read_excel('B.xlsx', dtype=str)
print(f"Read {len(HLA_B_ALLELES)} HLA B alleles from the Excel file.")
print(HLA_B_ALLELES.head())

# %%

AA_POSITIONS = [103]

# Create a new column for each AA position
for pos in AA_POSITIONS:

    # Create a new column with the AA at the given position, and skip if error
    HLA_B_ALLELES[f'B_{pos}'] = None  # Initialize the column with None

    def safe_getAA(allele):
        try:
            if pd.notna(allele) and isinstance(allele, str) and '*' in allele:
                return genie.getAA(allele, pos)
            else:
                return None
        except:
            return None
    
    HLA_B_ALLELES[f'B_{pos}'] = HLA_B_ALLELES['Haplotype'].apply(safe_getAA)

print(f"Assigned AA positions {AA_POSITIONS} to the DataFrame.")
print(HLA_B_ALLELES.head())

# %%

output_filename = "Alleles_B_with_AA_assignment.xlsx"
HLA_B_ALLELES.to_excel(output_filename, index=False)
print(f"Saved the updated DataFrame with AA assignments to {output_filename}.")
