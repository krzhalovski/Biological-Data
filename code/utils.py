import os
import pprint
import json
import pandas as pd
from math import sqrt

SWISSPROT_SIZE = 563_972
SWISSPROT_N_AA = 203_185_243
SWISSPROT_PATH = "../data/SwissProt/uniprot_sprot.fasta"
GROUND_TRUTH = pd.read_csv('reviewed_domain_positions.csv', index_col=1)

TOOLS_PATH = "../tools"
MODEL_PATH = "../models"
MSA_PATH = "../data/msa"



def build_hmm(msa_name, model_name):
    """
    Builds an hmm model from a t-coffee MSA.
    """
    if os.path.isdir(f'{MODEL_PATH}/hmm/{model_name}'):
        print("A model with this name exists, try again")
        return -1
    
    os.mkdir(f"{MODEL_PATH}/hmm/{model_name}")
    os.system(f"{TOOLS_PATH}/hmmer-3.3.1/src/hmmbuild {MODEL_PATH}/hmm/{model_name}/{model_name}.hmm {MSA_PATH}/{msa_name}")

def search_hmm(model_name):
    placeholder = f"{MODEL_PATH}/hmm/{model_name}/{model_name}"
    
    os.system(f"{TOOLS_PATH}/hmmer-3.3.1/src/hmmsearch --tblout {placeholder}.hmmer_tblout --domtblout {placeholder}.hmmer_domtblout {placeholder}.hmm {SWISSPROT_PATH} > {placeholder}.hmmer_align")

def build_pssm(msa_name, model_name):
    """
    Builds a pssm model from a t-coffee MSA
    """
    if os.path.isdir(f'{MODEL_PATH}/pssm/{model_name}'):
        print("A model with this name exists, try again")
        return -1
    os.mkdir(f"{MODEL_PATH}/pssm/{model_name}")
    os.system(f"{TOOLS_PATH}/ncbi-blast-2.11.0+/bin/psiblast -in_msa {MSA_PATH}/{msa_name} -db {SWISSPROT_PATH} -outfmt 7 -max_target_seqs 1500 -num_iterations 3 -out_pssm {MODEL_PATH}/pssm/{model_name}/{model_name}.pssm -save_pssm_after_last_round > {MODEL_PATH}/pssm/{model_name}/{model_name}.out")

def evaluate_matching_sequences(pr_dom):
    """
    Evaluates a model considering ground truth.
    """
    
    sp_ids = set(GROUND_TRUTH.index.values)
    pr_ids = set(pr_dom['name'].values)
    
    metrics = {
        "TP": 0, 
        "TN": 0,
        "FP": 0,
        "FN": 0,
    }
    
    for name in pr_ids:
        if name in sp_ids:
            metrics['TP'] += 1
        else:
            metrics['FP'] += 1 

    for name in sp_ids:
        if name not in pr_ids:
            metrics['FN'] += 1
        
    metrics['TN'] = SWISSPROT_SIZE - metrics['TP']

    TP, TN, FP, FN = metrics['TP'], metrics['TN'], metrics['FP'], metrics['FN']
    
    metrics['accuracy'] = (TP + TN)/(TP+TN+FP+FN)
    metrics['precision'] = TP/(TP+FP)
    metrics['sensitivity'] = TP/(TP+FN)
    metrics['specificity'] = TN/(TN+FP)
    metrics['MMC'] = (TP*TN + FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    metrics['f1_score'] = 2*TP/(2*TP+FP+FN)
    
    return metrics
    

def evaluate_matching_positions(pr_dom):
    TP = TN = FP = FN = 0
    sp_ids = set(GROUND_TRUTH.index.values)

    for _, x in pr_dom.iterrows():
        name = x['name']
        if name not in sp_ids:
            continue
        
        gt = GROUND_TRUTH.loc[name]
        gt_start, gt_end, hs_start, hs_end = gt['start'], gt['end'], x['start'], x['end']
        length = int(x['length'])

        for i in range(1, length+1):
            in_gt = i >= gt_start and i <= gt_end
            in_hs = i >= hs_start and i <= hs_end

            if in_gt and in_hs:
                TP+=1
            if in_gt and not in_hs:
                FN+=1
            if not in_gt and in_hs:
                FP+=1
    
    TN = SWISSPROT_N_AA - TP - FN - FP
    metrics = {
        "TP": TP,
        "TN": TN,
        "FP": FP,
        "FN": FN,
    }
    
    metrics['accuracy'] = (TP + TN)/(TP+TN+FP+FN)
    metrics['precision'] = TP/(TP+FP)
    metrics['sensitivity'] = TP/(TP+FN)
    metrics['specificity'] = TN/(TN+FP)
    metrics['MMC'] = (TP*TN + FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    metrics['f1_score'] = 2*TP/(2*TP+FP+FN)
    
    return metrics

def parse_domtblout(file_name, E=0.001):
    if not file_name.endswith('domtblout'):
        return -1
    
    entries = []
    with open(file_name) as f:
        for a in f.readlines():
            if a.startswith('sp'):
                splitted = a.split()
                
                e_val_pos = 6
                e_val = float(splitted[e_val_pos])
                
                if e_val >= E:
                    break
                    
                cols_to_save = [2, 5, 17, 18, 19, 20]
                
                values = [int(splitted[i]) for i in cols_to_save]
                name = splitted[0].split('|')[1]
                gt_len = values[0]
                gt_start, gt_end, hs_start, hs_end = values[2:]

                entries.append({
                    "name": name,
                    "length": gt_len,
                    "start": hs_start,
                    "end": hs_end
                })
            
    return pd.DataFrame(entries)

def parse_psiblast_output(file_name, E=0.001):
    if not file_name.endswith('out'):
        return -1
    
    ld = []
    columns = ['identity', 'length', 'mismatches', 'gap_opens', 'qs', 'qe', 'start', 'end', 'eval']
    
    with open(file_name) as f:
        for line in f.readlines()[6:-1]:
            splitted = line.split('\t')
            splitted[-1] = splitted[-1][:-1]

            ld.append({
                columns[i-2]: float(splitted[i]) for i in range(2, 11)
            })
            
            ld[-1]['name'] = splitted[1].split('|')[1]
        
    return pd.DataFrame(ld)

def complete_analysis_pssm(msa_name, model_name, Eval=1e-10, verbose=False):
    if verbose:
        print(f"1. Building PSSM Model from {msa_name}")
    out = build_pssm(msa_name, model_name)
    
    if out == -1:
        print("Canceling")
        return -1
    
    pr_dom = parse_psiblast_output(f"{MODEL_PATH}/pssm/{model_name}/{model_name}.out", E=Eval)
    
    if verbose:
        print(f"2. Evaluating Matching Sequences")
    sequences_metrics = evaluate_matching_sequences(pr_dom)
    
    with open(f'{MODEL_PATH}/pssm/{model_name}/sequences_metrics.json', "w") as f:
        json.dump(sequences_metrics, f)
        
    if verbose:
        print(f"3. Evaluating Domain Positions")
    domain_position_metrics = evaluate_matching_positions(pr_dom)
    
    with open(f'{MODEL_PATH}/pssm/{model_name}/positions_metrics.json', "w") as f:
        json.dump(domain_position_metrics, f)
    
    if verbose:
        print(f"Saved model in {MODEL_PATH}/pssm/{model_name}/")

def complete_analysis_hmm(msa_name, model_name, Eval=1e-10, verbose=False):
    if verbose:
        print(f"1. Building HMM Model from {msa_name}")
    out = build_hmm(msa_name, model_name)
    
    if out == -1:
        print("Canceling")
        return -1
        
    if verbose:
        print(f"2. Calling HMM Search")
    search_hmm(model_name)
    
    pr_dom = parse_domtblout(f"{MODEL_PATH}/hmm/{model_name}/{model_name}.hmmer_domtblout", E=Eval)
    
    if verbose:
        print(f"3. Evaluating Matching Sequences")
    sequences_metrics = evaluate_matching_sequences(pr_dom)
    
    with open(f'{MODEL_PATH}/hmm/{model_name}/sequences_metrics.json', "w") as f:
        json.dump(sequences_metrics, f)
        
    if verbose:
        print(f"4. Evaluating Domain Positions")
    domain_position_metrics = evaluate_matching_positions(pr_dom)
    
    with open(f'{MODEL_PATH}/hmm/{model_name}/positions_metrics.json', "w") as f:
        json.dump(domain_position_metrics, f)
    
    if verbose:
        print(f"Saved model in {MODEL_PATH}/hmm/{model_name}/")