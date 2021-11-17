import pandas as pd
import numpy as np
import math

def calculate_charge(mut):
  charge = 0
  charge -= mut.count('D')
  charge -= mut.count('E')
  charge += mut.count('K')
  charge += mut.count('R')
  return charge
def file2Dict(file,strain):
    df = pd.read_csv(file, sep='\t')
    df_strain = df[(df.strain == strain)][['ID', 'fit']]
    strain_fit_dic = df2Dict(df_strain, 'ID', 'fit')
    return strain_fit_dic

def df2Dict(DF,KEY,VALUE):
    Dict={}
    for index, row in DF.iterrows():
        ID=row[KEY]
        fit=row[VALUE]
        Dict[ID]=fit
    return Dict
def cal_Enrich_ab(Fa,Fb):
    enri=math.log10(Fb)-math.log10(Fa)
    return  enri
def calcu_trajectory(strain_fit_dict,resi_aa_dict,start_seq,enri=0):
    trajectory_dict={}
    trajectory_history_seq=start_seq.rsplit('-')[:-1]
    wt_seq=start_seq.rsplit('-')[-1]
    for i in range(len(wt_seq)):
        for n in resi_aa_dict[i]:
            if n==wt_seq[i]:continue
            mut_seq = wt_seq[:i] + n + wt_seq[i + 1:]
            if mut_seq in trajectory_history_seq:continue #mutation never go back
            charg = calculate_charge(mut_seq)
            if abs(charg) >= 2:continue # charge limitation
            mut_enri = cal_Enrich_ab(strain_fit_dict[wt_seq], strain_fit_dict[mut_seq])
            enri_final=mut_enri + enri
            seq = start_seq + '-' + mut_seq
            trajectory_dict[seq] = enri_final
    return trajectory_dict

def steps_trajectory_predict(steps,strain,resi_aa_dict,query_seq,output):
    trajectory_dict_total = {}
    trajectory_dict_total2 = {}
    trajectory_dict_total3 = {}
    final_state_dict={}
    all_state_dict={}
    strain_fit_dict = file2Dict('result/NA_compile_results.tsv', strain)
    trajectory_dict = calcu_trajectory(strain_fit_dict, resi_aa_dict, query_seq)
    if steps == 1:
        # trajectory_df = pd.DataFrame.from_dict(trajectory_dict, orient='index', columns=['Trajectory score'])
        final_state_dict = trajectory_dict
        all_state_dict = trajectory_dict
    if steps == 2:
        for id in trajectory_dict.keys():
            trajectory_dict2 = calcu_trajectory(strain_fit_dict, resi_aa_dict, id, trajectory_dict[id])
            trajectory_dict_total.update(trajectory_dict2)

        final_state_dict.update(trajectory_dict_total)
        all_state_dict.update(trajectory_dict)
        all_state_dict.update(trajectory_dict_total)

        # trajectory_df = pd.DataFrame.from_dict(trajectory_dict_total, orient='index', columns=['Trajectory score'])
    if steps == 3:
        for id in trajectory_dict.keys():
            trajectory_dict2 = calcu_trajectory(strain_fit_dict, resi_aa_dict, id, trajectory_dict[id])
            trajectory_dict_total.update(trajectory_dict2)
        for id in trajectory_dict_total.keys():
            trajectory_dict2 = calcu_trajectory(strain_fit_dict, resi_aa_dict, id, trajectory_dict_total[id])
            trajectory_dict_total2.update(trajectory_dict2)
        final_state_dict.update(trajectory_dict_total2)
        all_state_dict.update(trajectory_dict_total)
        all_state_dict.update(trajectory_dict)
        all_state_dict.update(trajectory_dict_total2)
        # trajectory_df = pd.DataFrame.from_dict(trajectory_dict_total2, orient='index', columns=['Trajectory score'])
    if steps == 4:
        for id in trajectory_dict.keys():
            trajectory_dict2 = calcu_trajectory(strain_fit_dict, resi_aa_dict, id, trajectory_dict[id])
            trajectory_dict_total.update(trajectory_dict2)
        for id in trajectory_dict_total.keys():
            trajectory_dict2 = calcu_trajectory(strain_fit_dict, resi_aa_dict, id, trajectory_dict_total[id])
            trajectory_dict_total2.update(trajectory_dict2)
        for id in trajectory_dict_total2.keys():
            trajectory_dict2 = calcu_trajectory(strain_fit_dict, resi_aa_dict, id, trajectory_dict_total2[id])
            trajectory_dict_total3.update(trajectory_dict2)
        final_state_dict.update(trajectory_dict_total3)
        all_state_dict.update(trajectory_dict_total)
        all_state_dict.update(trajectory_dict)
        all_state_dict.update(trajectory_dict_total2)
        all_state_dict.update(trajectory_dict_total3)
        # trajectory_df = pd.DataFrame.from_dict(trajectory_dict_total3, orient='index', columns=['Trajectory score'])
    write_trajectory_files(steps,output,final_state_dict,all_state_dict,strain_fit_dict,resi_aa_dict)

def calculate_p_dict(strain_fit_dict,resi_aa_dict,start_seq):
    P_dict={}
    P_history_seq = start_seq.rsplit('-')[:-2]
    wt_seq = start_seq.rsplit('-')[-2]
    for i in range(len(wt_seq)):
        for n in resi_aa_dict[i]:
            if n == wt_seq[i]: continue
            mut_seq = wt_seq[:i] + n + wt_seq[i + 1:]
            if mut_seq in P_history_seq: continue  # mutation never go back
            charg = calculate_charge(mut_seq)
            if abs(charg) >= 2: continue  # charge limitation
            mut_enri = cal_Enrich_ab(strain_fit_dict[wt_seq], strain_fit_dict[mut_seq])
            seq = wt_seq + '-' + mut_seq
            P_dict[seq]=mut_enri
    P_dict = {k: v / total for total in (sum(P_dict.values()),) for k, v in P_dict.items()}
    return P_dict
def write_trajectory_files(steps,file,final_state_dict,all_state_dict,strain_fit_dict,resi_aa_dict):
    outfile = open(file, 'w')
    outfile.write('step' + "\t" + 'E_score' + "\t" + 'one_step_Probability' + "\t" + "Trajectory" + "\n")
    for key in final_state_dict.keys():
        for x in range(steps + 1):
            if x == 0:
                wt_seq=key[:7]
                wt_fit=strain_fit_dict[wt_seq]
                y = math.log10(wt_fit)
                p = 1
            else:
                trajectory_seq=key[:(8 * (x + 1) - 1)]
                p_dict=calculate_p_dict(strain_fit_dict,resi_aa_dict,trajectory_seq)
                p=p_dict[trajectory_seq[-15:]]
                y = all_state_dict[trajectory_seq]
            outfile.write(str(x) + "\t" + str(y) + "\t" + str(p) +"\t" + key + "\n")
    outfile.close()

def main():

    resi_aa_dict={0:['N', 'K'],
              1:['N', 'S', 'D'],
              2:['K', 'E', 'R'],
              3:['N', 'S', 'G'],
              4:['K', 'E'],
              5:['K', 'E', 'D', 'T'],
              6:['L', 'S']}

    steps_trajectory_predict(1,'HK19',resi_aa_dict,"KSENETS","result/trajectory_prediction_HK19-1step.tsv")
    steps_trajectory_predict(1, 'Vic11', resi_aa_dict, "KNENETS", "result/trajectory_prediction_Vic11-1step.tsv")
    steps_trajectory_predict(3, 'Mos99', resi_aa_dict, "KNESEKL", "result/trajectory_prediction_Mos99-3steps.tsv")
    steps_trajectory_predict(3, 'Bei89', resi_aa_dict, "KNKGEEL", "result/trajectory_prediction_Bei89-3steps.tsv")
    steps_trajectory_predict(2, 'Bk79', resi_aa_dict, "KNKSEES", "result/trajectory_prediction_Bk79-2steps.tsv")
    steps_trajectory_predict(4, 'HK68', resi_aa_dict, "NDRSKDL", "result/trajectory_prediction_HK68-4steps.tsv")






if __name__ == "__main__":
  main()