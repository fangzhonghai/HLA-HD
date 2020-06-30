# -*- coding:utf-8 -*-
from functools import reduce
import pandas as pd
import subprocess
import optparse
import logging
import yaml
import time
import sys
import os
import re


def print_usage(option, opt, value, parser):
    usage_message = """
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    python3 HLA-HD.py --config hla-hd.yaml -s fq.ls -p /path/to/work --project project_name
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def yaml_read(yaml_file):
    with open(yaml_file, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    return yaml_dic


def get_sample_info(df):
    for sample in df['sample'].values:
        fq1 = df.loc[df['sample'] == sample, 'fq1'].values[0]
        fq2 = df.loc[df['sample'] == sample, 'fq2'].values[0]
        yield sample, fq1, fq2


def job_id_in_sge(command):
    command_status, command_output = subprocess.getstatusoutput(command)
    command_jobid = re.findall(r"Your job (\d+) ", command_output)[0]
    return command_status, command_jobid


def job_num_in_sge():
    command = "qstat | grep `whoami` |wc -l"
    command_status, command_output = subprocess.getstatusoutput(command)
    return command_status, int(command_output)


def job_status_in_sge(jobid):
    command = "qstat | grep " + "\"" + jobid + " " + "\""
    command_status, command_output = subprocess.getstatusoutput(command)
    return command_status, command_output


def jobs(script, flow, parent, resource):
    job_dic = {script: {}}
    job_dic[script]['flow'] = flow
    job_dic[script]['parent'] = parent
    job_dic[script]['child'] = []
    job_dic[script]['resource'] = resource
    job_dic[script]['jobid'] = ''
    if os.path.exists(script + '.check'):
        job_dic[script]['status'] = 'complete'
    else:
        job_dic[script]['status'] = 'incomplete'
    return job_dic


def jobs_map(job_dic) -> dict:
    job = list(filter(lambda x: job_dic[x]['flow'] != 'first_search', job_dic.keys()))
    for i in job:
        for j in job_dic[i]['parent']:
            job_dic[j]['child'].append(i)
    return job_dic


def create_hlahd(sample, fq1, fq2, conf, path):
    min_len = conf['parameters']['min_len_reads']
    mtsize = int(int(min_len)/2)
    n = conf['parameters']['nsize']
    cr = conf['parameters']['cut_rate']
    bowtie = conf['path']['bowtie2']
    hlahd = conf['path']['hlahd']
    threads = conf['resource']['threads']
    hla_gene_split = conf['path']['hla_gene_split']
    freq_dir = conf['path']['freq_dir']
    hla_dictionary = conf['path']['hla_dictionary']
    param_map = conf['parameters']['bowtie2_map']
    param_pmatch = conf['parameters']['bowtie2_pmatch']
    wd = os.path.join(path, sample + '/mapfile')
    wd_sample = os.path.join(path, sample)
    if not os.path.exists(wd):
        os.makedirs(wd)
    mapping = os.path.join(wd_sample, 'hlahd.' + sample + '.mapping')
    p_match = os.path.join(wd_sample, 'hlahd.' + sample + '.perfect_match')
    drop_reads = os.path.join(wd_sample, 'hlahd.' + sample + '.drop_reads')
    map_leaving = os.path.join(wd_sample, 'hlahd.' + sample + '.map_leaving')
    split_pm_read = os.path.join(wd_sample, 'hlahd.' + sample + '.split_pm_read')
    with open(mapping + '.1.sh', 'w') as f:
        m1 = jobs(mapping + '.1.sh', 'first_search', [], conf['resource']['mapping'])
        if not fq1.endswith('.gz'):
            text = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -x {hla_dictionary}/all_exon_intron_N150.fasta -U {fq1} -S {wd}/{sample}_all.R1.sam
if [ $? -ne 0]; then
    exit 1
fi
stfr -L {min_len} {wd}/{sample}_all.R1.sam {fq1} {wd}/{sample}.1
if [ $? -ne 0]; then
    exit 1
fi
mv {wd}/{sample}.1.fastq {wd}/{sample}.R1.fastq
touch {mapping}.1.sh.check
'''.format(**locals())
        else:
            text = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -x {hla_dictionary}/all_exon_intron_N150.fasta -U {fq1} -S {wd}/{sample}_all.R1.sam
if [ $? -ne 0]; then
    exit 1
fi
zcat {fq1} > {wd}/tmp_R1.fastq
if [ $? -ne 0]; then
    exit 1
fi
stfr -L {min_len} {wd}/{sample}_all.R1.sam {wd}/tmp_R1.fastq {wd}/{sample}.1
if [ $? -ne 0]; then
    exit 1
fi
mv {wd}/{sample}.1.fastq {wd}/{sample}.R1.fastq
touch {mapping}.1.sh.check
'''.format(**locals())
        f.write(text)
    with open(mapping + '.2.sh', 'w') as f:
        m2 = jobs(mapping + '.2.sh', 'first_search', [], conf['resource']['mapping'])
        if not fq2.endswith('.gz'):
            text = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -x {hla_dictionary}/all_exon_intron_N150.fasta -U {fq2} -S {wd}/{sample}_all.R2.sam
if [ $? -ne 0]; then
    exit 1
fi
stfr -L {min_len} {wd}/{sample}_all.R2.sam {fq2} {wd}/{sample}.2
if [ $? -ne 0]; then
    exit 1
fi
mv {wd}/{sample}.2.fastq {wd}/{sample}.R2.fastq
touch {mapping}.2.sh.check
'''.format(**locals())
        else:
            text = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -x {hla_dictionary}/all_exon_intron_N150.fasta -U {fq2} -S {wd}/{sample}_all.R2.sam
if [ $? -ne 0]; then
    exit 1
fi
zcat {fq2} > {wd}/tmp_R2.fastq
if [ $? -ne 0]; then
    exit 1
fi
stfr -L {min_len} {wd}/{sample}_all.R2.sam {wd}/tmp_R2.fastq {wd}/{sample}.2
if [ $? -ne 0]; then
    exit 1
fi
mv {wd}/{sample}.2.fastq {wd}/{sample}.R2.fastq
touch {mapping}.2.sh.check
'''.format(**locals())
        f.write(text)

    with open(p_match + '.exon.1.sh', 'w') as f:
        pe1 = jobs(p_match + '.exon.1.sh', 'perfect_match', [mapping + '.1.sh'], conf['resource']['mapping'])
        text1 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} {param_pmatch} -p {threads} -x {hla_dictionary}/all_exon_N150.fasta -U {wd}/{sample}.R1.fastq -S {wd}/{sample}.all_exon.R1.pmap.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 0 --NM 0 {wd}/{sample}.all_exon.R1.pmap.sam 1.0 > {wd}/{sample}.all_exon.R1.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {p_match}.exon.1.sh.check
'''.format(**locals())
        f.write(text1)
    with open(p_match + '.exon.2.sh', 'w') as f:
        pe2 = jobs(p_match + '.exon.2.sh', 'perfect_match', [mapping + '.2.sh'], conf['resource']['mapping'])
        text1 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} {param_pmatch} -p {threads} -x {hla_dictionary}/all_exon_N150.fasta -U {wd}/{sample}.R2.fastq -S {wd}/{sample}.all_exon.R2.pmap.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 0 --NM 0 {wd}/{sample}.all_exon.R2.pmap.sam 1.0 > {wd}/{sample}.all_exon.R2.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {p_match}.exon.2.sh.check
'''.format(**locals())
        f.write(text1)
    with open(p_match + '.intron.1.sh', 'w') as f:
        pi1 = jobs(p_match + '.intron.1.sh', 'perfect_match', [mapping + '.1.sh'], conf['resource']['mapping'])
        text1 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} {param_pmatch} -p {threads} -x {hla_dictionary}/all_intron_N150.fasta -U {wd}/{sample}.R1.fastq -S {wd}/{sample}.all_intron.R1.pmap.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 0 --NM 0 {wd}/{sample}.all_intron.R1.pmap.sam 1.0 > {wd}/{sample}.all_intron.R1.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {p_match}.intron.1.sh.check
'''.format(**locals())
        f.write(text1)
    with open(p_match + '.intron.2.sh', 'w') as f:
        pi2 = jobs(p_match + '.intron.2.sh', 'perfect_match', [mapping + '.2.sh'], conf['resource']['mapping'])
        text1 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} {param_pmatch} -p {threads} -x {hla_dictionary}/all_intron_N150.fasta -U {wd}/{sample}.R2.fastq -S {wd}/{sample}.all_intron.R2.pmap.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 0 --NM 0 {wd}/{sample}.all_intron.R2.pmap.sam 1.0 > {wd}/{sample}.all_intron.R2.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {p_match}.intron.2.sh.check
'''.format(**locals())
        f.write(text1)
    with open(p_match + '.merge.1.sh', 'w') as f:
        pm1 = jobs(p_match + '.merge.1.sh', 'perfect_match', [p_match + '.exon.1.sh', p_match + '.intron.1.sh'], conf['resource']['fq_merge'])
        text1 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
cat {wd}/{sample}.all_exon.R1.pmap.NM.sam {wd}/{sample}.all_intron.R1.pmap.NM.sam > {wd}/{sample}.all.R1.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
stfr {wd}/{sample}.all.R1.pmap.NM.sam {wd}/{sample}.R1.fastq {wd}/{sample}.R1.pm
if [ $? -ne 0]; then
    exit 1
fi
get_diff_fasta {wd}/{sample}.R1.fastq {wd}/{sample}.R1.pm.fastq > {wd}/{sample}.R1.diff.fastq
if [ $? -ne 0]; then
    exit 1
fi
touch {p_match}.merge.1.sh.check
'''.format(**locals())
        f.write(text1)
    with open(p_match + '.merge.2.sh', 'w') as f:
        pm2 = jobs(p_match + '.merge.2.sh', 'perfect_match', [p_match + '.exon.2.sh', p_match + '.intron.2.sh'], conf['resource']['fq_merge'])
        text1 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
cat {wd}/{sample}.all_exon.R2.pmap.NM.sam {wd}/{sample}.all_intron.R2.pmap.NM.sam > {wd}/{sample}.all.R2.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
stfr {wd}/{sample}.all.R2.pmap.NM.sam {wd}/{sample}.R2.fastq {wd}/{sample}.R2.pm
if [ $? -ne 0]; then
    exit 1
fi
get_diff_fasta {wd}/{sample}.R2.fastq {wd}/{sample}.R2.pm.fastq > {wd}/{sample}.R2.diff.fastq
if [ $? -ne 0]; then
    exit 1
fi
touch {p_match}.merge.2.sh.check
'''.format(**locals())
        f.write(text1)

    with open(drop_reads + '.exon.1.sh', 'w') as f:
        de1 = jobs(drop_reads + '.exon.1.sh', 'drop_reads', [p_match + '.merge.1.sh'], conf['resource']['mapping'])
        text2 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -x {hla_dictionary}/all_exon_N150.fasta -U {wd}/{sample}.R1.diff.fastq -S {wd}/{sample}.all_exon.R1.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 0 --NM {n} {wd}/{sample}.all_exon.R1.sam {cr} > {wd}/{sample}.all_exon.R1.diff.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {drop_reads}.exon.1.sh.check
'''.format(**locals())
        f.write(text2)
    with open(drop_reads + '.exon.2.sh', 'w') as f:
        de2 = jobs(drop_reads + '.exon.2.sh', 'drop_reads', [p_match + '.merge.2.sh'], conf['resource']['mapping'])
        text2 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -x {hla_dictionary}/all_exon_N150.fasta -U {wd}/{sample}.R2.diff.fastq -S {wd}/{sample}.all_exon.R2.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 0 --NM {n} {wd}/{sample}.all_exon.R2.sam {cr} > {wd}/{sample}.all_exon.R2.diff.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {drop_reads}.exon.2.sh.check
'''.format(**locals())
        f.write(text2)
    with open(drop_reads + '.intron.1.sh', 'w') as f:
        di1 = jobs(drop_reads + '.intron.1.sh', 'drop_reads', [p_match + '.merge.1.sh'], conf['resource']['mapping'])
        text2 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -x {hla_dictionary}/all_intron_N150.fasta -U {wd}/{sample}.R1.diff.fastq -S {wd}/{sample}.all_intron.R1.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 2 --NM {n} {wd}/{sample}.all_intron.R1.sam {cr} > {wd}/{sample}.all_intron.R1.diff.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {drop_reads}.intron.1.sh.check
'''.format(**locals())
        f.write(text2)
    with open(drop_reads + '.intron.2.sh', 'w') as f:
        di2 = jobs(drop_reads + '.intron.2.sh', 'drop_reads', [p_match + '.merge.2.sh'], conf['resource']['mapping'])
        text2 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -x {hla_dictionary}/all_intron_N150.fasta -U {wd}/{sample}.R2.diff.fastq -S {wd}/{sample}.all_intron.R2.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 2 --NM {n} {wd}/{sample}.all_intron.R2.sam {cr} > {wd}/{sample}.all_intron.R2.diff.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {drop_reads}.intron.2.sh.check
'''.format(**locals())
        f.write(text2)
    with open(drop_reads + '.merge.1.sh', 'w') as f:
        dm1 = jobs(drop_reads + '.merge.1.sh', 'drop_reads', [drop_reads + '.exon.1.sh', drop_reads + '.intron.1.sh'], conf['resource']['fq_merge'])
        text2 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
cat {wd}/{sample}.all_exon.R1.diff.pmap.NM.sam {wd}/{sample}.all_intron.R1.diff.pmap.NM.sam > {wd}/{sample}.all.R1.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
stfr {wd}/{sample}.all.R1.pmap.NM.sam {wd}/{sample}.R1.diff.fastq {wd}/{sample}.R1.diff2
if [ $? -ne 0]; then
    exit 1
fi
mv {wd}/{sample}.R1.diff2.fastq {wd}/{sample}.R1.diff.fastq
touch {drop_reads}.merge.1.sh.check
'''.format(**locals())
        f.write(text2)
    with open(drop_reads + '.merge.2.sh', 'w') as f:
        dm2 = jobs(drop_reads + '.merge.2.sh', 'drop_reads', [drop_reads + '.exon.2.sh', drop_reads + '.intron.2.sh'], conf['resource']['fq_merge'])
        text2 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
cat {wd}/{sample}.all_exon.R2.diff.pmap.NM.sam {wd}/{sample}.all_intron.R2.diff.pmap.NM.sam > {wd}/{sample}.all.R2.pmap.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
stfr {wd}/{sample}.all.R2.pmap.NM.sam {wd}/{sample}.R2.diff.fastq {wd}/{sample}.R2.diff2
if [ $? -ne 0]; then
    exit 1
fi
mv {wd}/{sample}.R2.diff2.fastq {wd}/{sample}.R2.diff.fastq
touch {drop_reads}.merge.2.sh.check
'''.format(**locals())
        f.write(text2)

    with open(map_leaving + '.exon.1.sh', 'w') as f:
        le1 = jobs(map_leaving + '.exon.1.sh', 'map_leaving', [drop_reads + '.merge.1.sh'], conf['resource']['mapping'])
        text3 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -a -x {hla_dictionary}/all_exon_N150.fasta -U {wd}/{sample}.R1.diff.fastq -S {wd}/{sample}.all_exon.R1.sam
if [ $? -ne 0]; then
    exit 1
fi
grep -v ^@ {wd}/{sample}.all_exon.R1.pmap.NM.sam >> {wd}/{sample}.all_exon.R1.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 0 --NM {n} {wd}/{sample}.all_exon.R1.sam {cr} > {wd}/{sample}.all_exon.R1.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {map_leaving}.exon.1.sh.check
'''.format(**locals())
        f.write(text3)
    with open(map_leaving + '.exon.2.sh', 'w') as f:
        le2 = jobs(map_leaving + '.exon.2.sh', 'map_leaving', [drop_reads + '.merge.2.sh'], conf['resource']['mapping'])
        text3 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -a -x {hla_dictionary}/all_exon_N150.fasta -U {wd}/{sample}.R2.diff.fastq -S {wd}/{sample}.all_exon.R2.sam
if [ $? -ne 0]; then
    exit 1
fi
grep -v ^@ {wd}/{sample}.all_exon.R2.pmap.NM.sam >> {wd}/{sample}.all_exon.R2.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 0 --NM {n} {wd}/{sample}.all_exon.R2.sam {cr} > {wd}/{sample}.all_exon.R2.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {map_leaving}.exon.2.sh.check
'''.format(**locals())
        f.write(text3)
    with open(map_leaving + '.intron.1.sh', 'w') as f:
        li1 = jobs(map_leaving + '.intron.1.sh', 'map_leaving', [drop_reads + '.merge.1.sh'], conf['resource']['mapping'])
        text3 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -a -x {hla_dictionary}/all_intron_N150.fasta -U {wd}/{sample}.R1.fastq -S {wd}/{sample}.all_intron.R1.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 2 --NM {n} {wd}/{sample}.all_intron.R1.sam {cr} > {wd}/{sample}.all_intron.R1.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {map_leaving}.intron.1.sh.check
'''.format(**locals())
        f.write(text3)
    with open(map_leaving + '.intron.2.sh', 'w') as f:
        li2 = jobs(map_leaving + '.intron.2.sh', 'map_leaving', [drop_reads + '.merge.2.sh'], conf['resource']['mapping'])
        text3 = '''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
bowtie2 {param_map} -p {threads} -a -x {hla_dictionary}/all_intron_N150.fasta -U {wd}/{sample}.R2.fastq -S {wd}/{sample}.all_intron.R2.sam
if [ $? -ne 0]; then
    exit 1
fi
pm_extract --MP 2 --NM {n} {wd}/{sample}.all_intron.R2.sam {cr} > {wd}/{sample}.all_intron.R2.NM.sam
if [ $? -ne 0]; then
    exit 1
fi
touch {map_leaving}.intron.2.sh.check
'''.format(**locals())
        f.write(text3)
    with open(split_pm_read + '.sh', 'w') as f:
        s = jobs(split_pm_read + '.sh', 'split_pm_read', [map_leaving + '.exon.1.sh', map_leaving + '.intron.1.sh', map_leaving + '.exon.2.sh', map_leaving + '.intron.2.sh'], conf['resource']['split_pm_read'])
        text4 = r'''#!/bin/bash
export PATH=$PATH:{bowtie}
export PATH=$PATH:{hlahd}
split_pm_read -f {freq_dir} -t {min_len} -m {mtsize} {wd}/{sample}.all_exon.R1.NM.sam {wd}/{sample}.all_exon.R2.NM.sam \
{wd}/{sample}.all_intron.R1.NM.sam {wd}/{sample}.all_intron.R2.NM.sam {hla_gene_split} \
{sample} {hla_dictionary} 150 {path}
if [ $? -ne 0]; then
    exit 1
fi
split_shell {wd_sample}/estimation.sh 29
if [ $? -ne 0]; then
    exit 1
fi
touch {split_pm_read}.sh.check
'''.format(**locals())
        f.write(text4)
    est = []
    est_scripts = []
    for i in range(29):
        est.append(jobs(wd_sample+'/estimation.sh.'+str(i), 'estimation', [split_pm_read + '.sh'], conf['resource']['estimation']))
        est_scripts.append(wd_sample+'/estimation.sh.'+str(i))
    pick_up = jobs(wd_sample+'/pickup.sh', 'pickup', est_scripts, conf['resource']['estimation'])
    jobs_list = [m1, m2, pe1, pe2, pi1, pi2, pm1, pm2, de1, de2, di1, di2, dm1, dm2, le1, le2, li1, li2, s, pick_up] + est
    jobs_dic = reduce(lambda x, y: dict(x, **y), jobs_list)
    return jobs_dic


def job_submit(job, job_resource, pro):
    job_num_status, job_num = job_num_in_sge()
    if job_num_status != 0:
        logger_main.error('qstat Failed!')
        sys.exit(1)
    while job_num >= 1999:
        time.sleep(600)
        job_num_status, job_num = job_num_in_sge()
        if job_num_status != 0:
            logger_main.error('qstat Failed!')
            sys.exit(1)
    job_path = os.path.dirname(job)
    command = "qsub -wd " + job_path + " -P " + pro + " " + job_resource + " " + job
    command_status, command_jobid = job_id_in_sge(command)
    return command_status, command_jobid


def add_env_in_head(script, conf):
    env = conf['path']['hlahd']
    text = '''#!/bin/bash
export PATH=$PATH:{env}
'''.format(**locals())
    with open(script, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(text + content)


def add_flag_in_script(script):
    text = '''if [ $? -ne 0]; then
exit 1
fi
touch {script}.sh.check
'''.format(**locals())
    with open(script, 'a+') as f:
        f.write(text)


def work_flow(jobs_dic, pro, conf):
    queue = list(filter(lambda x: jobs_dic[x]['flow'] == 'first_search', jobs_dic.keys()))
    while len(queue) > 0:
        queue_add = []
        queue_remove = []
        for q in queue:
            if os.path.exists(q + '.check'):
                jobs_dic[q]['status'] = 'complete'
            if jobs_dic[q]['status'] == 'incomplete':
                if jobs_dic[q]['jobid'] == '':
                    if len(jobs_dic[q]['parent']) == 0:
                        command_status, command_jobid = job_submit(q, jobs_dic[q]['resource'], pro)
                        if command_status != 0:
                            logger_main.error(q + ' Submit Failed!')
                            sys.exit(1)
                        jobs_dic[q]['jobid'] = command_jobid
                        logger_main.info(q + ' Submit Success!')
                    else:
                        status_list = [jobs_dic[p]['status'] for p in jobs_dic[q]['parent']]
                        if 'incomplete' not in status_list:
                            if q.__contains__('/estimation.sh.') or q.__contains__('/pickup.sh'):
                                add_env_in_head(q, conf)
                                add_flag_in_script(q)
                            command_status, command_jobid = job_submit(q, jobs_dic[q]['resource'], pro)
                            if command_status != 0:
                                logger_main.error(q + ' Submit Failed!')
                                sys.exit(1)
                            jobs_dic[q]['jobid'] = command_jobid
                            logger_main.info(q + ' Submit Success!')
                else:
                    command_status, command_output = job_status_in_sge(jobs_dic[q]['jobid'])
                    if command_status != 0 and command_output == '' and not os.path.exists(q + '.check'):
                        logger_main.error(q + ' Run Failed!')
                        sys.exit(1)
            else:
                queue_remove.append(q)
                queue_add += jobs_dic[q]['child']
                logger_main.info(q + ' Finished!')
        if len(queue_add) == 0:
            time.sleep(60)
        queue = list(set(list(set(queue) - set(queue_remove)) + queue_add))


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-s', dest='info', default=None, metavar='file')
    parser.add_option('-p', '--pwd', dest='pwd', default=None, metavar='string')
    parser.add_option('--project', dest='project', help='project queue in sge', default=None, metavar='string')
    parser.add_option('--config', dest='config', default=None, metavar='file')
    (opts, args) = parser.parse_args()
    pwd = opts.pwd
    info = opts.info
    config = opts.config
    project = opts.project
    config_dic = yaml_read(config)
    fh = logging.FileHandler(pwd + '/main.log', 'w', encoding='utf-8')
    logger_main = logging.getLogger()
    logger_main.setLevel(logging.INFO)
    fm = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    logger_main.addHandler(fh)
    fh.setFormatter(fm)
    logger_main.info('HLA-HD Pipeline Create Jobs Start!')
    info_df = pd.read_csv(info, sep='\t', header=None)
    info_df.columns = ['sample', 'fq1', 'fq2']
    hla = []
    for sampleid, fastq1, fastq2 in get_sample_info(info_df):
        hla.append(create_hlahd(sampleid, fastq1, fastq2, config_dic, pwd))
    hla_dic = reduce(lambda x, y: dict(x, **y), hla)
    hla_dic_final = jobs_map(hla_dic)
    with open(pwd + '/jobs_map.yaml', 'w') as fp:
        yaml.dump(hla_dic_final, fp, default_flow_style=False)
        logger_main.info('HLA-HD Pipeline Create Jobs & Map Finish!')
    work_flow(hla_dic_final, project, config_dic)
    logger_main.info('All Jobs Finished!')
