import networkx as nx
import numpy as np
import pandas as pd
import gml_analysis as gma
import deg_sequence_analysis as dsa
import categorize_networks as cn

if __name__ == "__main__":
    deg_seq_dir = './example_deg_seqs'

    analysis_df = dsa.analyze_deg_sequences(deg_seq_dir)
    result_df = cn.categorize_networks(analysis_df)
    result_df.to_excel('output.xlsx')