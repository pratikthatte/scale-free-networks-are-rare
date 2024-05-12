import numpy as np
import pandas as pd

def read_deg_seq(deg_seq_file_path):
    deg_seq_df = pd.read_csv(deg_seq_file_path)
    expanded_deg_seq  = np.concatenate([np.array([deg_seq_df.xvalue[i] for repeat in range(deg_seq_df.counts[i])]) for i in range(len(deg_seq_df))])
    return expanded_deg_seq