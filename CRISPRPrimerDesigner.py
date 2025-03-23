#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import regex as re
import os
import subprocess
import argparse
import primer3
import pandas as pd
import time

def get_header():
	return (r'''
			CRISPRPrimerDesigner
--A script that designs PCR or RPA primers for CRISPR nucleic acid detection.
Version v1.0.0
Last Revised:2025/3/23
''')

def Fasta_reverse(sequence):
	complement={'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return ''.join([complement.get(base, 'N') for base in reversed(sequence.upper())])

def get_sg(genome_file):
	rre=re.compile(r'TTT[ACG].{23}')
	rrv=re.compile(r'.{23}[TGC]AAA')
	seq={}
	with open(genome_file) as f:
		for line in f:
			line=line.strip()
			if line.startswith('>'):
				name=line.replace('>','').split()[0]
				seq[name]=[]
			else:
				seq[name].append(line.upper())

	with open('tmp.txt','w') as out:
		for key,value in seq.items():
			count=0
			value=''.join(value) 
			sgRNA_f=rre.finditer(value,overlapped=True)
			sgRNA_r=rrv.finditer(value,overlapped=True)
			for i in sgRNA_f:
				sss=i.group()
				count+=1
				start=i.start()
				end=i.end()
				name=key+'__F__'+str(count)
				out.write(f'{name}\t{start}\t{end}\t{sss}\n')
			count=0
			for i in sgRNA_r:
				sss=Fasta_reverse(i.group())
				count+=1
				start=i.start()
				end=i.end()
				name=key+'__R__'+str(count)
				out.write(f'{name}\t{start}\t{end}\t{sss}\n')
			count=0
	with open('tmp.txt') as ff:
		dic={}
		with open('tmpp.txt','w') as ot:
			for line in ff:
				if line.split('\t')[-1][:20] not in dic.keys():
					dic[line.split('\t')[-1][:20]]=line
			for value in dic.values():
				ot.write(value)
	os.remove('tmp.txt')
	with open('tmpp.txt') as ff:
		dic={}
		with open('sgRNA.txt','w') as ot:
			for line in ff:
				name=line.split('\t')[0].split('__')[0]+'__'+line.split('\t')[0].split('__')[1]
				if name not in dic.keys():
					dic[name]=[]
				dic[name].append(line)
			for x,y in dic.items():
				n=0
				for i in y:
					n+=1
					lst=i.split('\t')
					del(lst[0])
					word=x+'__'+str(n)+'\t'+'\t'.join(lst)
					ot.write(word)
	os.remove('tmpp.txt')

def get_casoffinder_input(genome_file):
	with open('sgRNA.txt') as f,open('casoffinder_input.txt','w') as ot:
		ot.write(genome_file+'\n'+'TTTNNNNNNNNNNNNNNNNNNNNNNNN')
		for line in f:
			line=line.strip().split()
			name=line[0]
			sg=line[-1]
			new_sg='NNNN'+sg[4:]
			ot.write('\n'+new_sg+' '+'5'+' '+name)

def run_casoffinder():
	cmd = ['./cas-offinder','casoffinder_input.txt','G0','casoffinder_output.txt']
	process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	stdout, stderr = process.communicate()
	if process.returncode != 0:
		raise RuntimeError(f"cas-offinder failed: {stderr.decode()}")

def get_index():
	ff=open('casoffinder_input.txt')
	f=ff.readlines()[2:]
	ff.close()
	dic={}
	mis_dic={}
	for line in f:
		sp=line.strip().split()
		seq=sp[0]
		num=int(sp[1])
		id_=sp[-1]
		dic[seq]=id_
		mis_dic[id_]={str(i):0 for i in range(num+1)}
	return dic,mis_dic

def count():
	with open('casoffinder_output.txt') as ff:
		for line in ff:
			sp=line.strip().split()
			seq=sp[0]
			mis=sp[-2]
			id_=sg_index_dict[seq]
			mismatch_dict[id_][mis]+=1

def make_result():
	out=open('calc_mismatch.txt','w')
	for id_,mis_dic in mismatch_dict.items():
		out.write(id_+'\t')
		values=list(mis_dic.values())
		values=[str(i) for i in values]
		out.write('\t'.join(values)+'\t')
		total=0
		for i in values[1:]:
			total+=int(i)
		out.write(str(total)+'\n')	
	out.close()

def work():
	global sg_index_dict
	global mismatch_dict
	sg_index_dict,mismatch_dict=get_index()
	count()
	make_result()

def combine_files():
	dic={}
	with open('calc_mismatch.txt') as f:
		for line in f:
			sp=line.strip().split()
			sg_id=sp[0]
			mismatch=','.join(sp[1:])
			dic[sg_id]=mismatch

	with open('sgRNA.txt') as ff,open('sg_offtarget.txt','w') as ot:
		for line in ff:
			sp=line.strip().split()
			sg_id=sp[0]
			ot.write(line.strip()+'\t'+dic[sg_id]+'\n')

def choose_low_offtarget_sg():
	with open('sg_low_offtarget.txt','w') as ot:
		set1=set()
		with open('sg_offtarget.txt') as f:
			for line in f:
				sp=line.strip().split()
				N23=sp[3][4:]
				mis=sp[-1].strip().split(',')
				mis0=mis[0]
				mistotal=mis[-1]
				if mis0=='1' and mistotal=='0':
					if N23 not in set1:
						set1.add(N23)
						ot.write(line)

def get_N227_seq(genome_file):
	seq={}
	with open(genome_file) as f:
		for line in f:
			line=line.strip()
			if line.startswith('>'):
				name=line.replace('>','').split()[0]
				seq[name]=[]
			else:
				seq[name].append(line.upper())

	segments={}
	with open('sg_low_offtarget.txt') as f:
		for line in f:
			sp=line.strip().split()
			seq_idd=sp[0]
			seq_id=sp[0].split('__')[0]
			ori=sp[0].split('__')[1]
			start=sp[1]
			end=sp[2]
			sg=sp[3]
			new_start=int(start)-200
			new_end=int(end)+200
			DNA_string=''.join(seq[seq_id])
			if new_start>0 and new_end<len(DNA_string):
				N427=DNA_string[new_start:new_end]
				gc_count=N427.count('G')+N427.count('C')
				gc_content=f'{gc_count/427:.2f}'
				if 0.4<=float(gc_content)<=0.6:
					segments[seq_idd]=N427
	return segments

def primer_design(genome_file,amplify_mode):
	ot=open('data_primer.txt','w')
	seq_args={}
	segments=get_N227_seq(genome_file)
	for key,value in segments.items():
		seq_args['SEQUENCE_ID']=key
		seq_args['SEQUENCE_TEMPLATE']=value
		seq_args['SEQUENCE_INCLUDED_REGION']=[0,len(value)]

		if amplify_mode=='PCR':
			global_args={
				'PRIMER_OPT_SIZE': 20,
				'PRIMER_MIN_SIZE': 18,
				'PRIMER_MAX_SIZE': 25,
				'PRIMER_OPT_TM': 60.0,
				'PRIMER_MIN_TM': 57.0,
				'PRIMER_MAX_TM': 63.0,
				'PRIMER_MIN_GC': 40.0,
				'PRIMER_MAX_GC': 60.0,
				'PRIMER_MAX_POLY_X': 100,
				'PRIMER_INTERNAL_MAX_POLY_X': 4,
				'PRIMER_SALT_MONOVALENT': 50.0,
				'PRIMER_DNA_CONC': 50.0,
				'PRIMER_MAX_NS_ACCEPTED': 0,
				'PRIMER_MAX_SELF_ANY': 8,
				'PRIMER_MAX_SELF_END': 4,
				'PRIMER_PAIR_MAX_COMPL_ANY': 8,
				'PRIMER_PAIR_MAX_COMPL_END': 4,
				'PRIMER_PRODUCT_SIZE_RANGE': [150,250],
				'PRIMER_NUM_RETURN':1,
				'PRIMER_MAX_HAIRPIN_TH':20,
			}
			primer3_result=primer3.bindings.design_primers(seq_args, global_args)

			if 'PRIMER_LEFT_0_SEQUENCE' in primer3_result:
				seq_idd=key
				pl_seq=primer3_result['PRIMER_LEFT_0_SEQUENCE']
				pl_coor_len=primer3_result['PRIMER_LEFT_0']
				pl_start=int(pl_coor_len[0])
				pl_tm=f"{primer3_result['PRIMER_LEFT_0_TM']:.2f}"
				pl_gc=f"{primer3_result['PRIMER_LEFT_0_GC_PERCENT']:.2f}"
				pr_seq=primer3_result['PRIMER_RIGHT_0_SEQUENCE']
				pr_coor_len=primer3_result['PRIMER_RIGHT_0']
				pr_start=int(pr_coor_len[0])
				pr_tm=f"{primer3_result['PRIMER_RIGHT_0_TM']:.2f}"
				pr_gc=f"{primer3_result['PRIMER_RIGHT_0_GC_PERCENT']:.2f}"
				product_size=primer3_result['PRIMER_PAIR_0_PRODUCT_SIZE']
				if 80<pl_start<190 or 230<pr_start<370:
					ot.write(f'{seq_idd}\t{pl_seq}\t{pl_coor_len}\t{pl_tm}\t{pl_gc}\t{pr_seq}\t{pr_coor_len}\t{pr_tm}\t{pr_gc}\t{product_size}\n')
			seq_args={}

		elif amplify_mode=='RPA':
			global_args={
				'PRIMER_OPT_SIZE': 30,
				'PRIMER_MIN_SIZE': 30,
				'PRIMER_MAX_SIZE': 35,
				'PRIMER_OPT_TM': 65.0,
				'PRIMER_MIN_TM': 60.0,
				'PRIMER_MAX_TM': 70.0,
				'PRIMER_MIN_GC': 40.0,
				'PRIMER_MAX_GC': 60.0,
				'PRIMER_MAX_POLY_X': 100,
				'PRIMER_INTERNAL_MAX_POLY_X': 4,
				'PRIMER_SALT_MONOVALENT': 50.0,
				'PRIMER_DNA_CONC': 50.0,
				'PRIMER_MAX_NS_ACCEPTED': 0,
				'PRIMER_MAX_SELF_ANY': 8,
				'PRIMER_MAX_SELF_END': 4,
				'PRIMER_PAIR_MAX_COMPL_ANY': 8,
				'PRIMER_PAIR_MAX_COMPL_END': 4,
				'PRIMER_PRODUCT_SIZE_RANGE': [150,250],
				'PRIMER_NUM_RETURN':1,
				'PRIMER_MAX_HAIRPIN_TH':20,
			}
			primer3_result = primer3.bindings.design_primers(seq_args, global_args)

			if 'PRIMER_LEFT_0_SEQUENCE' in primer3_result:
				seq_idd=key
				pl_seq=primer3_result['PRIMER_LEFT_0_SEQUENCE']
				pl_coor_len=primer3_result['PRIMER_LEFT_0']
				pl_start=int(pl_coor_len[0])
				pl_tm=f"{primer3_result['PRIMER_LEFT_0_TM']:.2f}"
				pl_gc=f"{primer3_result['PRIMER_LEFT_0_GC_PERCENT']:.2f}"
				pr_seq=primer3_result['PRIMER_RIGHT_0_SEQUENCE']
				pr_coor_len=primer3_result['PRIMER_RIGHT_0']
				pr_start=int(pr_coor_len[0])
				pr_tm=f"{primer3_result['PRIMER_RIGHT_0_TM']:.2f}"
				pr_gc=f"{primer3_result['PRIMER_RIGHT_0_GC_PERCENT']:.2f}"
				product_size=primer3_result['PRIMER_PAIR_0_PRODUCT_SIZE']
				if 80<pl_start<190 or 230<pr_start<370:
					ot.write(f'{seq_idd}\t{pl_seq}\t{pl_coor_len}\t{pl_tm}\t{pl_gc}\t{pr_seq}\t{pr_coor_len}\t{pr_tm}\t{pr_gc}\t{product_size}\n')
			seq_args={}
	ot.close()

def get_primer_output(output):
	sg_dic={}
	with open('sgRNA.txt') as f:
		for line in f:
			sp=line.strip().split()
			idd=sp[0]
			iddd=sp[0].split('__')[0]
			start=sp[1]
			end=sp[2]
			sg=sp[3]
			info=f'{sg}\t{iddd}:{start}-{end}'
			sg_dic[idd]=info
	with open('data_primer.txt') as ff,open(output,'w') as out:
		out.write(f'index\tsgRNA\tseqid:start-end\tstrand\tleft_primer\t[primer_start,length]\tTM\tGC%\tright_primer\t[primer_start,length]\tTM\tGC%\tproduct_size\n')
		index=0
		for line in ff:
			index+=1
			sp=line.strip().split('\t')
			idd=sp[0]
			ori=sp[0].split('__')[1]
			pl_seq=sp[1]
			pl_coor_len=sp[2].replace(' ','')
			pl_tm=sp[3]
			pl_gc=sp[4]
			pr_seq=sp[5]
			pr_coor_len=sp[6].replace(' ','')
			pr_tm=sp[7]
			pr_gc=sp[8]
			product_size=sp[9]
			if ori=='F':
				out.write(f'{index}\t{sg_dic[idd]}\t+\t{pl_seq}\t{pl_coor_len}\t{pl_tm}\t{pl_gc}\t{pr_seq}\t{pr_coor_len}\t{pr_tm}\t{pr_gc}\t{product_size}\n')
			else:
				out.write(f'{index}\t{sg_dic[idd]}\t-\t{pl_seq}\t{pl_coor_len}\t{pl_tm}\t{pl_gc}\t{pr_seq}\t{pr_coor_len}\t{pr_tm}\t{pr_gc}\t{product_size}\n')

def main():
	start_time = time.time()
	parser=argparse.ArgumentParser(description='CRISPRPrimerDesigner Parameters',formatter_class=argparse.ArgumentDefaultsHelpFormatter,
		usage='python CRISPRPrimerDesigner_v1.py -g <genome_file> [-a <choose amplify mode: PCR or RPA>] [-o <output file>]')
	parser.add_argument('-g','--genome',type=str,required=True,help='genome fasta file path')
	parser.add_argument('-a','--amplify_mode',type=str,choices=['PCR', 'RPA'],required=False,help='amplify mode choose (PCR/RPA)',default='PCR')
	parser.add_argument('-o','--output',type=str,required=False,help='output text file path',default='primer_output.txt')
	args=parser.parse_args()

	genome_file,amplify_mode,output=args.genome,args.amplify_mode,args.output
	print(get_header())
	get_sg(genome_file)
	get_casoffinder_input(genome_file)
	run_casoffinder()
	work()
	combine_files()
	choose_low_offtarget_sg()
	primer_design(genome_file,amplify_mode)
	get_primer_output(output)

	os.remove('sgRNA.txt')
	os.remove('casoffinder_input.txt')
	os.remove('casoffinder_output.txt')
	os.remove('calc_mismatch.txt')
	os.remove('sg_offtarget.txt')
	os.remove('sg_low_offtarget.txt')
	os.remove('data_primer.txt')

	end_time=time.time()
	total_seconds=end_time - start_time
	minutes, seconds = divmod(total_seconds, 60)
	print('Task done!')
	print(f'\nTotal runtime: {int(minutes)}m {seconds:.2f}s')


if __name__=='__main__':
	main()
