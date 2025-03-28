import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from openpyxl import load_workbook

excel_file = 'FIBERS_SRTR_AAMM_Bin_1_OVERALL_summary_table_v3.xlsx'
sheet_name = 'FIBERS_PARAMS'  
plot_filename = 'Pareto_Front_Plot_v3.png'

df = pd.read_excel(excel_file, sheet_name=sheet_name)

fig, ax = plt.subplots(figsize=(10, 6))
scatter = ax.scatter(df['Group Ratio'], df['Log-Rank Score'], c='blue')

# Label axes
ax.set_xlabel('Group Ratio')
ax.set_ylabel('Log-Rank Score')
ax.set_title('Pareto Front Plot')

# Add bin labels
for i, txt in enumerate(df['FIBERS Bin']):
    ax.annotate(txt, (df['Group Ratio'][i], df['Log-Rank Score'][i]), 
                xytext=(5, 5), textcoords='offset points')

# Add Adjusted HR labels
for i, txt in enumerate(df['Adjusted HR']):
    ax.annotate(f'{txt:.2f}', (df['Group Ratio'][i], df['Log-Rank Score'][i]), 
                xytext=(5, -10), textcoords='offset points')

# Function to identify Pareto front points
def is_pareto_efficient(costs):
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
        if is_efficient[i]:
            is_efficient[is_efficient] = np.any(costs[is_efficient]<c, axis=1)  # Keep any point with a lower cost
            is_efficient[i] = True  # And keep self
    return is_efficient

# Identify Pareto front points
costs = np.column_stack((-df['Group Ratio'], -df['Log-Rank Score']))  # Negative because we want to maximize
pareto_front = is_pareto_efficient(costs)

# Draw dotted line between Pareto front points
pareto_points = df[pareto_front].sort_values('Group Ratio')
ax.plot(pareto_points['Group Ratio'], pareto_points['Log-Rank Score'], 'r--', linewidth=2)

plt.tight_layout()
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
#plt.show()
plt.close()