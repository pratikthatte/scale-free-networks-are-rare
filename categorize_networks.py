import numpy as np
import pandas as pd

def test_strong(rows):
    S1 = False
    S2 = False
    SA = False
    number_of_deg_seqs = len(rows)
    strong = 0
    strong_alone = 0
    for i, row in rows.iterrows():
        if row.power_law_p_value > 0.1 and row.alpha < 3 and row.alpha > 2:
            strong_alone += 1
            if row.exponential_decision > -1 and row.lognormal_decision > -1 and row.truncated_power_law_decision > -1 and row.stretched_exponential_decision > -1:
                strong +=1
    if strong_alone >= 0.9*number_of_deg_seqs:
            SA = True
    if strong >= 0.5*number_of_deg_seqs:
            S2=True
    if SA and strong >= 0.9*number_of_deg_seqs:
            S1=True
    return (S1,S2)

def test_weak(rows):
    W = False
    West = False
    SW = False
    n = len(rows)
    weak = 0
    weakest = 0
    sweak = 0
    for i, row in rows.iterrows():
        if row.power_law_p_value>0.1:
            weakest += 1
            if row.n_tail>=50:
                weak += 1
        if row.exponential_decision != 0 and row.lognormal_decision != 0 and row.truncated_power_law_decision != 0 and row.stretched_exponential_decision != 0:
            sweak += 1
    if weak >= 1:
        W = True
    if weakest >= 1:
        West = True
    if sweak >= 1:
        SW = True
    return (W, West, SW)

def categorize_networks(deg_seq_data_df):
    unique_network_name = np.unique(deg_seq_data_df.network_name)
    print(unique_network_name)
    results = pd.DataFrame(columns = ['Strongest', "Strong", "Weak", "Weakest",
                                       "Super_Weak", "median_alpha",
                                       "median_ntail"], index=unique_network_name )
    for i, network in enumerate(unique_network_name):
        query = "network_name == '%s'" %network
        rows = deg_seq_data_df.query(query)
        [strongest, strong] = test_strong(rows)
        [weak, weakest, superweak] = test_weak(rows)
        results.loc[network]['Strongest'] = strongest
        print("values for network: "+str(network))
        results.loc[network]['Strong'] = strong
        results.loc[network]['Weak'] = weak
        results.loc[network]['Weakest']=weakest
        results.loc[network]['Super_Weak'] = superweak
        results.loc[network]['median_alpha'] = np.median(rows.alpha)
        results.loc[network]['median_ntail']=np.median(rows.n_tail)
    return results