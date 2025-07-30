import os


step=10000
start_index_list=[0,10000,20000,30000,40000,50000,60000,70000,80000,90000]#parallelize
conn_p_list=[0.01]
is_smoothing_list=['yes']

# Female
save_dir='../Raw_module/'#save raw path
DIR_DATA='../data/entity_weight_matchit_0407.csv'#comorbidity and weight
DIR_NET='../data/LCC_0407.csv'#iBKH net

for is_smoothing in is_smoothing_list:
    for conn_p in conn_p_list:
        for start_index in start_index_list:
            cmd=f'python module_generation.py {start_index} {step} {conn_p} {is_smoothing} {save_dir}{is_smoothing}Smooth/ {DIR_DATA} {DIR_NET} &'
            os.system(cmd)
