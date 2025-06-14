#!/usr/bin/env python3
import json
import subprocess
import argparse

def bcftools_norm(IN_VCF):
    folder_list = IN_VCF.split('/')
    name = folder_list[-1].replace('.vcf', '')
    split_vcf = name + '_norm.vcf'
    cmd = f'''
    bcftools norm \
    -m -any \
    -Ov -o {split_vcf} \
    {IN_VCF}
    '''
    subprocess.run(cmd, shell=True)
    return split_vcf

def CHECK_ERROR_JSON(IN_JSON):
    check_key = True
    f1 = open(IN_JSON, 'r')
    list1 = f1.read().split('\n')
    f1.close()
    head = list1.pop(0)
    list1.pop(-1)
    tail = list1.pop(-1)
    if head != '{':
        print('first line is not {')
    if tail != '}':
        print('last line is not }')
    body_list = list1
    for i in range(len(body_list)):
        if 'FILTER' not in body_list[i] and 'INFO' not in body_list[i] and 'FORMAT' not in body_list[i]:
            print(f'There is not any column name in {body_list[i]}')
            check_key = False
        elif 'FILTER.' not in body_list[i] and 'INFO.' not in body_list[i] and 'FORMAT.' not in body_list[i]:
            print(f'There is not any dot between column name and tag in {body_list[i]}')
            check_key = False
        elif ': ' not in body_list[i] and ':' in body_list[i]:
            print(f'Space is not exist after ":" in {body_list[i]}')
            check_key = False
        elif '>' not in body_list[i] and '<' not in body_list[i] and '=' not in body_list[i]:
            print(f'numeric comparisons is not exist in {body_list[i]}')
            check_key = False
    return check_key

def json_to_bcftools_expr(IN_JSON):
    with open(IN_JSON) as f:
        rules = json.load(f)
    
    parts = []
    for key, cond in rules.items():
        tag = key.replace('.', '/')
        parts.append(f"{tag}{cond}")
    
    expr = " && ".join(parts)
    
    return expr

def bcftools_filter(expr, IN_VCF, OUT_VCF):
    cmd = f'''
    bcftools filter \
    -i "{expr}" \
    --soft-filter PASS \
    --mode x \
    -Ov -o {OUT_VCF} \
    {IN_VCF}
    '''
    subprocess.run(cmd, shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf')
    parser.add_argument('-j', '--json')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    input_vcf = args.vcf
    input_json = args.json
    output_vcf = args.output
    if CHECK_ERROR_JSON(input_json) == True:
        norm_vcf = bcftools_norm(input_vcf)
        expr = json_to_bcftools_expr(input_json)
        bcftools_filter(expr, norm_vcf, output_vcf)
    print('done')
