import numpy as np
import pandas as pd
import os
import fit
import import_files as im

def analyze_deg_sequences(deg_sequences_path):
    file_path_list = []
    deg_seq_data_df = pd.DataFrame(columns=['network_name', 'alpha','xmin','n_tail','ks','power_law_p_value','exponential_decision','lognormal_decision','truncated_power_law_decision','stretched_exponential_decision'])
    for root, dirs, files in os.walk(deg_sequences_path):
        for file_name in files:
            file_path_list.append(os.path.join(root,file_name))
    for file_path in file_path_list:
        print("File Path is: "+file_path)
        network_name = file_path.split("/")[-1].split(".")[0]
        print("network name is: "+network_name)
        deg_seq = im.read_deg_seq(file_path)
        mean_degree = np.mean(deg_seq)
        if mean_degree<2 or mean_degree>np.sqrt(len(deg_seq)):
            continue
        xmin,alpha,n_tail,ks, exponential_decision, lognormal_decision, stretched_exponential_decision, truncated_power_law_decision = fit.fit_power_law_and_compare_with_alternative(deg_seq)
        p = fit.power_law_pval(deg_seq, alpha, int(xmin),ks)
        deg_seq_data_df.loc[file_path] = ""
        deg_seq_data_df.loc[file_path]['network_name'] = network_name
        print("Value in df: " + deg_seq_data_df.loc[file_path]['network_name'])
        deg_seq_data_df.loc[file_path]['alpha'] = alpha
        deg_seq_data_df.loc[file_path]['xmin'] = xmin
        deg_seq_data_df.loc[file_path]['n_tail'] = n_tail
        deg_seq_data_df.loc[file_path]['power_law_p_value']=p
        deg_seq_data_df.loc[file_path]['exponential_decision']= exponential_decision
        deg_seq_data_df.loc[file_path]['lognormal_decision']= lognormal_decision
        deg_seq_data_df.loc[file_path]['truncated_power_law_decision']= truncated_power_law_decision
        deg_seq_data_df.loc[file_path]['stretched_exponential_decision']= stretched_exponential_decision
    return deg_seq_data_df






