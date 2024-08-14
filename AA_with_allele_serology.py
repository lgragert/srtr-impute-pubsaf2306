import pandas as pd
from openpyxl.styles import PatternFill
from openpyxl import load_workbook
import hlagenie

# Go from serologic group to allele to amino acid pairings

# Get amino acid position with hlagenie: DRB1(16, 28, 30, 32, 37, 47, 57, 60, 67, 70, 71, 74, 85, 86)
#                                        DQB1 (9, 13, 23, 26, 30, 37, 38, 57, 66, 67, 70, 86, and 87)
def aminoacid(ag_allele, which_allele):
    if which_allele == "DRB1":
        AA_list = [16, 28, 30, 32, 37, 47, 57, 60, 67, 70, 71, 74, 85, 86]
    if which_allele == "DQB1":
        AA_list = [9, 13, 23, 26, 30, 37, 38, 57, 66, 67, 70, 86, 87]

    for AA in AA_list:
        D_B1_AA = []
        for row in range(len(ag_allele)):
            allele = ag_allele['HLA_Allele'][row]
            D_B1_AA.append(genie.getAA(allele, AA))
        ag_allele[which_allele + "_" + str(AA)] = D_B1_AA
    return ag_allele


# Transpose dataset to get the data to be what you want
def get_allele_format(optn_ag, ag_alleles):

    # Only keep the first ten DR alleles for now
    ag_alleles = ag_alleles['IMGT/HLA alleles'].str.split(pat='/', expand=True)
    ag_alleles = ag_alleles.iloc[:, :10]
    ag_alleles = pd.concat([optn_ag, ag_alleles], axis=1).reset_index(drop=True)

    # Want the rows to be indexed by allele
    transpose_ag_allele = ag_alleles.transpose()
    transpose_ag_allele.columns = transpose_ag_allele.iloc[0]
    transpose_ag_allele = transpose_ag_allele.iloc[1:]
    melt_ag_allele = pd.melt(transpose_ag_allele)
    melt_ag_allele = melt_ag_allele.dropna(ignore_index=True)
    melt_ag_allele = melt_ag_allele.rename(columns={'value': 'HLA_Allele'})
    return melt_ag_allele


# Two blank lines between each serological group
def blank_lines_insert(allele_aa, optn_ag_list):
    index_list = []
    for serology in optn_ag_list:
        index_list = allele_aa[(allele_aa["OPTN_Antigen"] == serology)].index.values
        # last_index_str = ''.join(map(str, index_list[-1:]))
        # last_index = int(last_index_str)
        last_index = index_list[-1]
        blank_line = pd.DataFrame(index=[0, 1], columns=allele_aa.columns)
        allele_aa = pd.concat([allele_aa.loc[:last_index], blank_line, allele_aa.loc[last_index + 1:]]).reset_index(drop=True)
    allele_aa = allele_aa[:-2]
    return allele_aa


# First start with list of allele serology with antigen type
ag_to_alleles = pd.read_csv('OPTN_antigens_to_alleles_CPRA.txt', sep='\t')

DRB1_ag_alleles = ag_to_alleles[(ag_to_alleles['IMGT/HLA alleles'].str.contains('DRB1'))]
DR_optn_ag = DRB1_ag_alleles['OPTN_Antigen']
DR_optn_ag_list = DR_optn_ag.values.tolist()

DRB1_alleles = get_allele_format(DR_optn_ag, DRB1_ag_alleles)

DQB1_ag_alleles = ag_to_alleles[(ag_to_alleles['IMGT/HLA alleles'].str.contains('DQB1'))]
DQ_optn_ag = DQB1_ag_alleles['OPTN_Antigen']
DQ_optn_ag_list = DQ_optn_ag.values.tolist()

DQB1_alleles = get_allele_format(DQ_optn_ag, DQB1_ag_alleles)


# Use HLAGenie to get specific amino acids at certain positions
genie = hlagenie.init("3510", imputed=True)

allele_aa_DRB1 = aminoacid(DRB1_alleles, "DRB1")
allele_aa_DQB1 = aminoacid(DQB1_alleles, "DQB1")
print(allele_aa_DRB1)
print(allele_aa_DQB1)


# Putting blank lines inbetween each serologic group
allele_aa_DRB1 = blank_lines_insert(allele_aa_DRB1, DR_optn_ag_list)
allele_aa_DQB1 = blank_lines_insert(allele_aa_DQB1, DQ_optn_ag_list)
color = "spring"

with pd.ExcelWriter("AA_polymorphism_with_allele_" + color + ".xlsx") as writer:
    allele_aa_DRB1.to_excel(writer, index=False, sheet_name='DRB1 Alleles & AA')
    allele_aa_DQB1.to_excel(writer, index=False, sheet_name='DQB1 Alleles & AA')
# color coding for the 20 amino acids by characteristics
# Polar = yellow [S, T, N, Q, C],           index colors: [26, 34, 43, 47, 51]
# + Charged = green [K, R, H],              index colors: [50, 42, 57]
# - Charged = purple [D, E],                index colors: [46, 24]
# Aromatic = pink [F, Y, W],                index colors: [45, 29, 33]
# Aliphatic = blue [G, A, V, M, I, L, P],   index colors: [27, 44, 40, 35, 41, 15, 49]


# Use to get amino acid color for each worksheet
def amino_acid_color(worksheet):
    max_column = worksheet.max_column
    for col in range(2, max_column + 1):
        column_letter = worksheet.cell(row=2, column=col).column_letter
        for cell in worksheet[column_letter]:
            if cell.value in AA_colors.keys():
                color_AA = PatternFill(patternType='solid', fgColor=AA_colors.get(cell.value))
                cell.fill = color_AA
            else:
                continue

    worksheet.column_dimensions['B'].width = 13
    return worksheet


wb = load_workbook("AA_polymorphism_with_allele_" + color + ".xlsx")
ws_DR = wb['DRB1 Alleles & AA']
ws_DQ = wb['DQB1 Alleles & AA']


# # Each AA color to fill in - with clustalx color scheme, these are the aRGB codes (different from colors in lines 70-75)
# AA_colors = {"S": "FF19CC19", "T": "FF19CC19", "N": "FF19CC19", "Q": "FF19CC19", "C": "FFE57F7F", "K": "FFE53319",
#              "R": "FFE53319", "H": "FF19B2B2", "D": "FFCC4CCC", "E": "FFCC4CCC", "F": "FF197FE5", "Y": "FF19B2B2",
#              "W": "FF197FE5", "G": "FFE5994C", "A": "FF197FE5", "V": "FF197FE5", "M": "FF197FE5", "I": "FF197FE5",
#              "L": "FF197FE5", "P": "FFCCCC00"}

# # # Spring color scheme
AA_colors = {"A": "ffababab", "C": "ff0893fe", "D": "fffc5977", "E": "ffff9c83", "F": "ff02f8b7", "G": "fff66c17",
             "H": "ffd2b907", "I": "ff47beff", "K": "ffffc56d", "L": "ff358b8f", "M": "ff4b9e93", "N": "ffeb684b",
             "P": "fffd75cf", "Q": "ffdeb16f", "R": "ffb38b2b", "S": "ffbb8187", "T": "ffd3add1", "V": "ff519db9",
             "W": "ff00d619", "Y": "ff38a549"}


ws_DR = amino_acid_color(ws_DR)
ws_DQ = amino_acid_color(ws_DQ)

wb.save(filename="AA_polymorphism_with_allele_" + color + ".xlsx")

