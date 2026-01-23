
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')

replicate = "1"
filename = f"SRTR_AA_MM_matrix_epmm_fibers_grffail_{replicate}.txt.gz"

df = pd.read_csv(filename, sep='\t', compression='gzip')
print(f"Dataset loaded: {len(df):,} transplant pairs")
print(df.columns.tolist())

risk_cols = []
if 'Wiebe_Risk' in df.columns:
    risk_cols.append('Wiebe_Risk')
    print(f"Wiebe_Risk: {list(df['Wiebe_Risk'].unique())}")
if 'Bin_1_FIBERS_Risk' in df.columns:
    risk_cols.append('Bin_1_FIBERS_Risk')
    print(f"Bin_1_FIBERS_Risk: {list(df['Bin_1_FIBERS_Risk'].unique())}")

binary_cols = ['Wiebe_Low', 'Wiebe_Medium', 'Wiebe_High']
existing_binary = [col for col in binary_cols if col in df.columns]
if existing_binary:
    print(f"Other Wiebe columns: {existing_binary}")

outcome_cols = ['graftyrs', 'grf_fail', 'patyrs', 'died_dt']
existing_outcomes = [col for col in outcome_cols if col in df.columns]

if 'Wiebe_Risk' in df.columns:
    wiebe_counts = df['Wiebe_Risk'].value_counts()
    total = len(df)
    print("EpMM Categories and counts:")
    for category, count in wiebe_counts.items():
        pct = count/total*100
        print(f"     • {category}: {count:,} patients ({pct:.1f}%)")

if 'Bin_1_FIBERS_Risk' in df.columns:
    fibers_counts = df['Bin_1_FIBERS_Risk'].value_counts()
    total = len(df)
    print("FIBERS Categories and counts:")
    for category, count in fibers_counts.items():
        pct = count/total*100
        print(f"     • {category}: {count:,} patients ({pct:.1f}%)")

print(f"\nGraft failure:")
if 'grf_fail' in df.columns:
    grf_counts = df['grf_fail'].value_counts()
    total = len(df)
    failed = grf_counts.get(1, 0)
    print(f"   • Graft failures: {failed:,} out of {total:,} ({failed/total*100:.1f}%)")
    print(f"   • Functioning grafts: {total-failed:,} ({(total-failed)/total*100:.1f}%)")

# %%

if 'Wiebe_Risk' in df.columns and 'grf_fail' in df.columns:
    print(f"\nEplet GF:")
    wiebe_rates = df.groupby('Wiebe_Risk')['grf_fail'].agg(['count', 'sum', 'mean'])
    for risk_level in wiebe_rates.index:
        total = wiebe_rates.loc[risk_level, 'count']
        events = wiebe_rates.loc[risk_level, 'sum']
        rate = wiebe_rates.loc[risk_level, 'mean'] * 100
        print(f"     • {risk_level} risk: {events}/{total} ({rate:.1f}% graft failure)")

if 'Bin_1_FIBERS_Risk' in df.columns and 'grf_fail' in df.columns:
    print(f"\nFIBERS GF:")
    fibers_rates = df.groupby('Bin_1_FIBERS_Risk')['grf_fail'].agg(['count', 'sum', 'mean'])
    for risk_level in fibers_rates.index:
        total = fibers_rates.loc[risk_level, 'count']
        events = fibers_rates.loc[risk_level, 'sum']
        rate = fibers_rates.loc[risk_level, 'mean'] * 100
        print(f"     • {risk_level} risk: {events}/{total} ({rate:.1f}% graft failure)")

# %%

if 'Wiebe_Risk' in df.columns and 'Bin_1_FIBERS_Risk' in df.columns:
    crosstab = pd.crosstab(df['Wiebe_Risk'], df['Bin_1_FIBERS_Risk'], 
                            rownames=['Wiebe'], colnames=['FIBERS'])
    print(crosstab)
    
    exact_agreement = (df['Wiebe_Risk'] == df['Bin_1_FIBERS_Risk']).sum()
    total = len(df[(df['Wiebe_Risk'].notna()) & (df['Bin_1_FIBERS_Risk'].notna())])
    agreement_pct = exact_agreement / total * 100 if total > 0 else 0
    
    print(f"\nAgreement: {agreement_pct:.1f}%")
    print(f"     • Both methods assign same risk level: {exact_agreement:,} out of {total:,} patients")
    
    wiebe_unique = df['Wiebe_Risk'].unique()
    fibers_unique = df['Bin_1_FIBERS_Risk'].unique()

    if len(wiebe_unique) >= 2 and len(fibers_unique) >= 2:
        wiebe_high = (df['Wiebe_Risk'] == 'high').astype(int)
        fibers_high = (df['Bin_1_FIBERS_Risk'] == 'high').astype(int)
        
        both_high = ((wiebe_high == 1) & (fibers_high == 1)).sum()
        both_nonhigh = ((wiebe_high == 0) & (fibers_high == 0)).sum()
        binary_agreement = (both_high + both_nonhigh) / total * 100 if total > 0 else 0
        
        print(f"\nHigh vs Non-High agreement: {binary_agreement:.1f}%")
        print(f"     • Both classify as 'high': {both_high:,} patients")
        print(f"     • Both classify as 'non-high': {both_nonhigh:,} patients")
        
        wiebe_high_fibers_low = ((wiebe_high == 1) & (fibers_high == 0)).sum()
        wiebe_low_fibers_high = ((wiebe_high == 0) & (fibers_high == 1)).sum()
        
        print(f"\nDisagreement:")
        print(f"     • Wiebe 'high' but FIBERS 'low': {wiebe_high_fibers_low:,} patients")
        print(f"     • Wiebe 'non-high' but FIBERS 'high': {wiebe_low_fibers_high:,} patients")
        
# %%

def _plot_survival_curves(risk_var, title):
    
    plt.figure(figsize=(10, 6))
    
    risk_groups = sorted(df[risk_var].dropna().unique())
    
    print(f"\n{title.upper()} SURVIVAL CURVES:")
    
    for group in risk_groups:
        group_data = df[df[risk_var] == group]
        valid_data = group_data[(group_data['graftyrs'] > 0) & (group_data['graftyrs'].notna())]
        
        if len(valid_data) == 0:
            continue
        
        kmf = KaplanMeierFitter()
        kmf.fit(durations=valid_data['graftyrs'], 
                event_observed=valid_data['grf_fail'],
                label=f'{group} risk (n={len(valid_data)})')
        
        kmf.plot_survival_function()
        
        median_survival = kmf.median_survival_time_
        events = valid_data['grf_fail'].sum()
        print(f"     • {group} risk: {events} failures, median survival = {median_survival:.1f} years")
    
    plt.title(f'Graft Survival by {title}')
    plt.xlabel('Years from Transplant')
    plt.ylabel('Graft Survival Probability')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
        
    if len(risk_groups) == 2:
        group1_data = df[df[risk_var] == risk_groups[0]]
        group2_data = df[df[risk_var] == risk_groups[1]]
        
        group1_valid = group1_data[(group1_data['graftyrs'] > 0) & (group1_data['graftyrs'].notna())]
        group2_valid = group2_data[(group2_data['graftyrs'] > 0) & (group2_data['graftyrs'].notna())]
        
        if len(group1_valid) > 0 and len(group2_valid) > 0:
            result = logrank_test(group1_valid['graftyrs'], group2_valid['graftyrs'],
                                group1_valid['grf_fail'], group2_valid['grf_fail'])
            print(f"Log-rank test ({risk_groups[0]} vs {risk_groups[1]}): p = {result.p_value:.4f}")
            if result.p_value < 0.05:
                print("Significant difference between risk groups")
            else:
                print("No significant difference between risk groups")
    
    elif len(risk_groups) == 3:

        from itertools import combinations
        
        for group1, group2 in combinations(risk_groups, 2):
            group1_data = df[df[risk_var] == group1]
            group2_data = df[df[risk_var] == group2]
            
            group1_valid = group1_data[(group1_data['graftyrs'] > 0) & (group1_data['graftyrs'].notna())]
            group2_valid = group2_data[(group2_data['graftyrs'] > 0) & (group2_data['graftyrs'].notna())]
            
            if len(group1_valid) > 10 and len(group2_valid) > 10:
                result = logrank_test(group1_valid['graftyrs'], group2_valid['graftyrs'],
                                    group1_valid['grf_fail'], group2_valid['grf_fail'])
                
                sig_marker = "(significant)" if result.p_value < 0.05 else ""
                print(f"        • {group1} vs {group2}: p = {result.p_value:.4f}{sig_marker}")
        
    else:
        print(f"Unexpected number of risk groups: {len(risk_groups)}")
        print(f"Groups found: {risk_groups}")
        
    if len(risk_groups) == 3:
        _recommend_analysis_approach(risk_var, risk_groups)

def _recommend_analysis_approach(risk_var, risk_groups):
    if risk_var in df.columns and 'grf_fail' in df.columns:
        print(f"\nEvent rates:")
        rates = df.groupby(risk_var)['grf_fail'].agg(['count', 'sum', 'mean'])
        
        for group in risk_groups:
            if group in rates.index:
                total = rates.loc[group, 'count']
                events = rates.loc[group, 'sum'] 
                rate = rates.loc[group, 'mean'] * 100
                print(f"           • {group}: {rate:.1f}% failure rate (n={total})")
        
        if 'low' in rates.index and 'interm' in rates.index and 'high' in rates.index:
            low_rate = rates.loc['low', 'mean']
            interm_rate = rates.loc['interm', 'mean'] 
            high_rate = rates.loc['high', 'mean']
            
            if low_rate < interm_rate < high_rate:
                print(f"\nKeep all 3 groups")
            elif abs(interm_rate - low_rate) < abs(interm_rate - high_rate):
                print(f"\nCombine interm + low vs high")
            else:
                print(f"\nCombine high + interm vs low")

if 'graftyrs' in df.columns and 'grf_fail' in df.columns:
    
    if 'Wiebe_Risk' in df.columns:
        _plot_survival_curves('Wiebe_Risk', 'Eplet Risk')
    
    if 'Bin_1_FIBERS_Risk' in df.columns:
        _plot_survival_curves('Bin_1_FIBERS_Risk', 'FIBERS Risk')

# %%

def _create_dummies_with_reference(df, risk_var, reference_category='low'):
    
    df_encoded = pd.get_dummies(df, columns=[risk_var], prefix=risk_var, drop_first=False)
    
    all_dummy_cols = [col for col in df_encoded.columns if col.startswith(risk_var)]
    
    reference_col = f"{risk_var}_{reference_category}"
    
    if reference_col in df_encoded.columns:
        df_encoded = df_encoded.drop(columns=[reference_col])
        reference_used = reference_category
        print(f"Dropped '{reference_col}' - '{reference_category}' is reference group")
    else:
        available_categories = sorted([col.replace(f'{risk_var}_', '') for col in all_dummy_cols])
        print(f"'{reference_col}' not found in data")
        print(f"Available categories: {available_categories}")
        
        fallback_ref = available_categories[0]
        fallback_col = f"{risk_var}_{fallback_ref}"
        df_encoded = df_encoded.drop(columns=[fallback_col])
        reference_used = fallback_ref
        print(f"Using '{fallback_ref}' as reference instead")
        print(f"Dropped '{fallback_col}' - '{fallback_ref}' is reference group")
    
    final_dummy_cols = [col for col in df_encoded.columns if col.startswith(risk_var)]
    print(f"Final dummy columns: {final_dummy_cols}")
    print(f"Reference group: {reference_used} (HR = 1.00)")
    
    return df_encoded, reference_used
    
def _run_cph_separate(df, risk_var, reference_category='low'):

    print(f"\n{risk_var} Cox model:")
    
    model_df = df[['graftyrs', 'grf_fail', risk_var]].dropna()
    
    if len(model_df) == 0:
        print("No valid data for analysis")
        return None
    
    model_df_encoded, reference_used = _create_dummies_with_reference(
        model_df, risk_var, reference_category
    )
    
    print(f"Created dummy variables: {[col for col in model_df_encoded.columns if col.startswith(risk_var)]}")
    print(f"Reference group: {reference_used} (HR = 1.00)")

    cph = CoxPHFitter()
    cph.fit(model_df_encoded, duration_col='graftyrs', event_col='grf_fail')
    
    c_index = cph.concordance_index_
    print(f"C-index: {c_index:.3f}")
    
    print(f"\nHazard Ratios:")
    
    for var_name in cph.summary.index:
        if var_name.startswith(risk_var):
            hr = np.exp(cph.summary.loc[var_name, 'coef'])
            ci_lower = np.exp(cph.summary.loc[var_name, 'coef lower 95%'])
            ci_upper = np.exp(cph.summary.loc[var_name, 'coef upper 95%'])
            p_val = cph.summary.loc[var_name, 'p']
            
            level = var_name.replace(f'{risk_var}_', '')
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
            comparison = f"vs {reference_used}"
            print(f"        • {level} {comparison}: HR = {hr:.2f} (95% CI: {ci_lower:.2f}-{ci_upper:.2f}), p = {p_val:.4f} {sig}")
    
    return {'c_index': c_index, 'model': cph}

def _run_cph_combined(df):
    model_vars = ['graftyrs', 'grf_fail', 'Wiebe_Risk', 'Bin_1_FIBERS_Risk']
    model_df = df[model_vars].dropna()
    
    if len(model_df) < 50:
        print("Insufficient data for combined analysis")
        return
    
    print(f"        Sample size: {len(model_df):,} patients")
    
    wiebe_levels = len(model_df['Wiebe_Risk'].unique()) if 'Wiebe_Risk' in model_df.columns else 0
    fibers_levels = len(model_df['Bin_1_FIBERS_Risk'].unique()) if 'Bin_1_FIBERS_Risk' in model_df.columns else 0
    print(f"        Eplet levels: {wiebe_levels}, FIBERS levels: {fibers_levels}")
    
    model_df_encoded = pd.get_dummies(model_df, columns=['Wiebe_Risk', 'Bin_1_FIBERS_Risk'], 
                                    prefix=['Wiebe_Risk', 'FIBERS_Risk'], drop_first=False)
    
    dropped_refs = []
    if 'Wiebe_Risk_low' in model_df_encoded.columns:
        model_df_encoded = model_df_encoded.drop(columns=['Wiebe_Risk_low'])
        dropped_refs.append('Wiebe: low')
    elif any(col.startswith('Wiebe_Risk_') for col in model_df_encoded.columns):
        wiebe_cols = [col for col in model_df_encoded.columns if col.startswith('Wiebe_Risk_')]
        first_wiebe = sorted(wiebe_cols)[0]
        model_df_encoded = model_df_encoded.drop(columns=[first_wiebe])
        dropped_refs.append(f'Wiebe: {first_wiebe.replace("Wiebe_Risk_", "")}')
    
    if 'FIBERS_Risk_low' in model_df_encoded.columns:
        model_df_encoded = model_df_encoded.drop(columns=['FIBERS_Risk_low'])
        dropped_refs.append('FIBERS: low')
    elif any(col.startswith('FIBERS_Risk_') for col in model_df_encoded.columns):
        fibers_cols = [col for col in model_df_encoded.columns if col.startswith('FIBERS_Risk_')]
        first_fibers = sorted(fibers_cols)[0]
        model_df_encoded = model_df_encoded.drop(columns=[first_fibers])
        dropped_refs.append(f'FIBERS: {first_fibers.replace("FIBERS_Risk_", "")}')
    
    print(f"Reference groups: {', '.join(dropped_refs)}")
    
    cph = CoxPHFitter()
    cph.fit(model_df_encoded, duration_col='graftyrs', event_col='grf_fail')
    
    print(f"Combined C-index: {cph.concordance_index_:.3f}")
    
    print(f"Eplet Risk Categories (vs reference):")
    for var in cph.summary.index:
        if var.startswith('Wiebe_Risk_'):
            hr = np.exp(cph.summary.loc[var, 'coef'])
            ci_lower = np.exp(cph.summary.loc[var, 'coef lower 95%'])
            ci_upper = np.exp(cph.summary.loc[var, 'coef upper 95%'])
            p_val = cph.summary.loc[var, 'p']
            
            level = var.replace('Wiebe_Risk_', '')
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
            print(f"        • {level}: HR = {hr:.2f} (95% CI: {ci_lower:.2f}-{ci_upper:.2f}), p = {p_val:.4f} {sig}")
    
    print(f"FIBERS Risk Categories (vs reference):")
    for var in cph.summary.index:
        if var.startswith('FIBERS_Risk_'):
            hr = np.exp(cph.summary.loc[var, 'coef'])
            ci_lower = np.exp(cph.summary.loc[var, 'coef lower 95%'])
            ci_upper = np.exp(cph.summary.loc[var, 'coef upper 95%'])
            p_val = cph.summary.loc[var, 'p']
            
            level = var.replace('FIBERS_Risk_', '')
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
            print(f"        • {level}: HR = {hr:.2f} (95% CI: {ci_lower:.2f}-{ci_upper:.2f}), p = {p_val:.4f} {sig}")
    
    print(f"\nResults:")
    wiebe_sig_vars = [var for var in cph.summary.index if var.startswith('Wiebe_Risk_') and cph.summary.loc[var, 'p'] < 0.05]
    fibers_sig_vars = [var for var in cph.summary.index if var.startswith('FIBERS_Risk_') and cph.summary.loc[var, 'p'] < 0.05]
    
    if wiebe_sig_vars and fibers_sig_vars:
        print("Both methods contribute independent risk information")
        print("Categories from both methods remain significant when combined")
    elif wiebe_sig_vars:
        print("Only eplet categories remain significant when both included")
    elif fibers_sig_vars:
        print("Only FIBERS categories remain significant when both included")
    else:
        print("No categories significant when both methods included")
        print("Possible multicollinearity or overlapping information")

# %%

if 'graftyrs' in df.columns and 'grf_fail' in df.columns:
    
    analysis_df = df[(df['graftyrs'] > 0) & (df['graftyrs'].notna())].copy()
    
    univariate_results = {}
    
    if 'Wiebe_Risk' in analysis_df.columns:
        print(f"\nEplet Risk only")
        univariate_results['Wiebe'] = _run_cph_separate(analysis_df, 'Wiebe_Risk')
    
    if 'Bin_1_FIBERS_Risk' in analysis_df.columns:
        print(f"\nFIBERS Risk only")
        univariate_results['FIBERS'] = _run_cph_separate(analysis_df, 'Bin_1_FIBERS_Risk')
    
# %%

if 'graftyrs' in df.columns and 'grf_fail' in df.columns:
    
    _run_cph_combined(analysis_df)
    
# %%


def _calculate_c_index_categorical(df, risk_var):
    
    model_df = df[['graftyrs', 'grf_fail', risk_var]].dropna()
    
    if len(model_df) == 0:
        return 0.5
    
    model_df_encoded = pd.get_dummies(model_df, columns=[risk_var], prefix=risk_var, drop_first=False)
    
    reference_col = f"{risk_var}_low"
    if reference_col in model_df_encoded.columns:
        model_df_encoded = model_df_encoded.drop(columns=[reference_col])
    else:
        dummy_cols = [col for col in model_df_encoded.columns if col.startswith(risk_var)]
        if dummy_cols:
            first_col = sorted(dummy_cols)[0]
            model_df_encoded = model_df_encoded.drop(columns=[first_col])
    
    cph = CoxPHFitter()
    cph.fit(model_df_encoded, duration_col='graftyrs', event_col='grf_fail')
    
    return cph.concordance_index_

if 'graftyrs' in df.columns and 'grf_fail' in df.columns:
    
    analysis_df = df[(df['graftyrs'] > 0) & (df['graftyrs'].notna())].copy()
    
    comparison_results = {}
    
    print(f"\nComparison:")
    
    if 'Wiebe_Risk' in analysis_df.columns:
        wiebe_c = _calculate_c_index_categorical(analysis_df, 'Wiebe_Risk')
        comparison_results['Wiebe'] = wiebe_c
        print(f"Eplet C-index: {wiebe_c:.3f}")
    
    if 'Bin_1_FIBERS_Risk' in analysis_df.columns:
        fibers_c = _calculate_c_index_categorical(analysis_df, 'Bin_1_FIBERS_Risk')
        comparison_results['FIBERS'] = fibers_c
        print(f"FIBERS C-index: {fibers_c:.3f}")
    
    if 'Wiebe' in comparison_results and 'FIBERS' in comparison_results:
        diff = comparison_results['FIBERS'] - comparison_results['Wiebe']
        print(f"\nResults:")
        print(f"   • Difference (FIBERS - Wiebe): {diff:+.3f}")
        
        if abs(diff) < 0.02:
            print("Both methods have similar performance")
        elif diff > 0.02:
            print("FIBERS AAMM")
        else:
            print("Eplet MM Wiebe")
        