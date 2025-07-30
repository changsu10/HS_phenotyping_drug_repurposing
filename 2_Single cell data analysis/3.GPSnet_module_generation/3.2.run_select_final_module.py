import os

cutoff_list=[0.005]
trait_list=['v1']
conn_p_list=[0.01]
is_smoothing=['yesSmooth']
steps='0,10000,20000,30000,40000,50000,60000,70000,80000,90000'
summary_df='../summary.tsv'

for t in trait_list:
    for cutoff in cutoff_list:
        for subfolder in is_smoothing:
            save_final_module_path=f'../GPSnet_result/{subfolder}/'
            save_final_module_with_score_path=f'../GPSnet_result_keep_score/{subfolder}/'
            cmd='python select_final_module.py %s %s %s %f %s %s %s' % (t,subfolder,steps,cutoff,save_final_module_path,save_final_module_with_score_path,summary_df)
            os.system(cmd)
    print('finish '+t)

