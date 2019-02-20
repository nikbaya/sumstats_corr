#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 13:01:15 2019

Creates correlation matrices from both sumstats and genotype matrix tables.

Most recent version is listed first.

@author: nbaya
"""

import hail as hl
import numpy as np
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hail.linalg import BlockMatrix

"""
╔═══════════╗
║ Version 2 ║
╚═══════════╝
Use ld_matrix for both correlation matrices. Keep at ~13m variants instead of
filtering to HM3 variants.
"""
variant_set = 'qc_pos'
variants = hl.read_table('gs://nbaya/split/'+variant_set+'_variants.ht') # for hm3: import table hapmap3_variants.tsv'

def sumstats_correlation(chr_list):
    """
    Get LD correlation matrix using genotypes of white British, using
    variant_set variants.
    """
    starttime = datetime.datetime.now()
    
    mt0 = hl.read_matrix_table('gs://phenotype_31063/hail/gwas.imputed_v3.both_sexes.annotated.mt')
    mt1 = mt0.filter_rows(hl.is_defined(variants[mt0.locus,mt0.alleles])) #filter to variant_set variants
    mt1.describe()
    
    for ch in chr_list:
        mt_chr = hl.filter_intervals(mt1,[hl.parse_locus_interval(str(ch))])
        print(mt_chr.count_rows())
        ld = hl.ld_matrix(mt_chr.beta, mt_chr.locus, radius=3e7)
        ld_sparse = ld.sparsify_band(lower=1, upper = 1000)
        ld_sparse.write('gs://nbaya/sumstats_corr/hm3_ss_correlation_chr'+str(ch)+'.bm',overwrite=True)
    
    endtime = datetime.datetime.now()
    elapsed = endtime-starttime
    print('\n####################')
    print('Iteration time: '+str(round(elapsed.seconds/60, 2))+' minutes')
    print('####################')

def genotype_correlation(chr_list):
    """
    Get classic LD correlation matrix from genotypes of white British, using
    variant_set variants.
    """
    starttime = datetime.datetime.now()
    
    gt1 = hl.read_matrix_table('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_repart.mt')
        
    print(gt1.count())
    print(gt1.describe())
    
    for ch in chr_list:
        mt_chr = hl.filter_intervals(gt1,[hl.parse_locus_interval(str(ch))])
        print(mt_chr.count_rows())
        print(mt_chr.describe())
        ld = hl.ld_matrix(mt_chr.dosage, mt_chr.locus, radius=3e7)
        ld_sparse = ld.sparsify_band(lower=1, upper = 1000)
        ld_sparse.write('gs://nbaya/sumstats_corr/hm3_gt_correlation_chr'+str(ch)+'.bm',overwrite=True)
    
    endtime = datetime.datetime.now()
    elapsed = endtime-starttime
    print('\n####################')
    print('Iteration time: '+str(round(elapsed.seconds/60, 2))+' minutes')
    print('####################')

def plot_correlation_matrices(chr_list):
    """
    Plot combined correlation matrices for genotype-correlation and 
    sumstats-correlation matrices
    """
    for ch in chr_list:
        ss_ch = BlockMatrix.read('gs://nbaya/sumstats_corr/'+variant_set+'_ss_correlation_chr{}.bm/'.format(ch))
        gt_ch = BlockMatrix.read('gs://nbaya/sumstats_corr/'+variant_set+'_gt_correlation_chr{}.bm/'.format(ch))
        M_max = int(1e4)                             #max number of variants to be taken from the block matrices (suggested: 2e4)
        M = ss_ch.shape[0]                      #dimension of block matrix
#        for idx in range(int(M/M_max)+1):       #index of which disjoint window we are looking at in the block matrix
        for idx in range(0,int(M/M_max)+1):       #index of which disjoint window we are looking at in the block matrix
            M0 = M_max*(idx)               #start variant index for block matrix filtering
            M1 = min(M_max*(idx+1),M)      #stop variant index for block matrix filtering
            ss_np = ss_ch[M0:M1,M0:M1].to_numpy()
            gt_np = gt_ch[M0:M1,M0:M1].to_numpy()
            print('\nStarting variant window: ['+str(M0)+','+str(M1)+']')
            w = int(5e3)                            #window width of variants for correlation matrix (suggested: 2e3)
            for i in range(int((M1-M0-1)/w)+1):
                w0 = w*i                        #start variant index for window of correlation matrix
                w1 = min(w*(i+1),M1-M0)         #stop variant index for window of correlation matrix
                full = (ss_np[w0:w1,w0:w1]+
                        gt_np[w0:w1,w0:w1].T)
                np.fill_diagonal(full,1)
                fig,ax = plt.subplots()
                ax.imshow(full,cmap='bwr')
                ax.plot([0, w],[0, w],'k--',alpha=0.5,lw=2)
                plt.xlim([0,w])
                plt.ylim([w,0])
                ax.text(w*0.83,w*0.1,"SS",fontsize=60,alpha=0.5)
                ax.text(w*0.02,w*0.97,"GT",fontsize=60,alpha=0.5)
                plt.title('chr'+str(ch)+' '+variant_set+' variants ('+str(M0+w0)+'-'+str(M0+w1)+')')
                fig=plt.gcf()
                fig.set_size_inches(10,10)
                path=('gs://nbaya/sumstats_corr/plots/chr'+str(ch)+'_'+variant_set+
                      '_'+str(M0+w0).zfill(len(str(M)))+'-'+str(M0+w1).zfill(len(str(M)))+'.png')
                with hl.hadoop_open(path, 'wb') as f: 
                    fig.savefig(f,dpi=600)
                plt.close()
            print('\nFinished variant window: ['+str(M0)+','+str(M1)+']')
        
if __name__ == "__main__":
    chr_list = [22]
    #sumstats_correlation(chr_list)
    #genotype_correlation(chr_list)
    plot_correlation_matrices(chr_list)









"""
╔═══════════╗
║ Version 1 ║
╚═══════════╝
"""

        ##mt0 = hl.read_matrix_table('gs://phenotype_31063/hail/gwas.imputed_v3.both_sexes.annotated.mt')
        #variant_set = 'hm3'
        ##variants = hl.read_table('gs://nbaya/split/'+variant_set+'_variants.ht') # for hm3: import table hapmap3_variants.tsv'
        ##mt = mt0.filter_rows(hl.is_defined(variants[mt0.locus,mt0.alleles])) #filter to variant_set variants
        #
        ##chr_idx0 = [0]*22
        ##for ch in range(1,23):
        ##    mt_chr = hl.filter_intervals(mt,[hl.parse_locus_interval(str(ch))])
        ##    M = mt_chr.count_rows()
        ##    print('chr '+str(ch)+': '+str(M))
        ##    chr_idx0[ch-1] = M
        ##    
        ##chr_idx = np.cumsum(chr_idx0).tolist()
        ##chr_idx.insert(0,0)
        #
        #
        #
        ##if variant_set == 'hm3':
        ##chr_idx = [0, 90487, 182237, 258449, 326550, 396844, 462960, 523307, 582695,
        ##           632758, 691334, 747508, 801921, 843314, 879600, 912540, 946268,
        ##           976018, 1008524, 1028993, 1057641, 1073259, 1089172]
        ##
        ##X_hm3 = BlockMatrix.read('gs://nbaya/sumstats_corr/sumstats_hm3.bm')
        ##X_hm3 = X_hm3.T
        ##print(X_hm3.shape)
        ##
        ##for ch in range(1,23):
        ##    chr_bm = X_hm3[:,chr_idx[ch-1]:chr_idx[ch]]
        ##    M = chr_bm.shape[1]
        ##    print('chr '+str(ch)+': '+str(M))
        ##    X = chr_bm
        ##    CX = X-X.sum(axis=0)/X.shape[0] #mean-centered
        ##    cov = (1/X.shape[0])*CX.T@CX
        ##    D = (cov**(-1/2)).sparsify_band(lower=0,upper=0)
        ##    X_s = CX@D
        ##    R = (1/X.shape[0])*X_s.T@X_s
        ##    R.sparsify_band(lower=0,upper=1000).write('gs://nbaya/sumstats_corr/hm3_correlation_chr'+str(ch)+'.bm',overwrite=True)
        #
        ################################################################################
        #"""
        #Calculate LD matrix from genotypes
        #"""
        #
        ##gt0 = hl.read_matrix_table('gs://phenotype_31063/hail/imputed/ukb31063.GT.autosomes.mt/')
        #variant_set = 'hm3'
        ##variants = hl.read_table('gs://nbaya/split/'+variant_set+'_variants.ht') # for hm3: import table hapmap3_variants.tsv'
        ##gt = gt0.filter_rows(hl.is_defined(variants[gt0.locus,gt0.alleles])) #filter to variant_set variants
        #gt = hl.read_matrix_table('gs://nbaya/split/ukb31063.'+variant_set+'_variants.gwas_samples_repart.mt')
        #
        ##for ch in range(1,23)[::-1]:
        #for ch in range(1,22)[::-1]:
        #    print('starting chr '+str(ch))
        #    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
        #    gt_chr = hl.filter_intervals(gt,[hl.parse_locus_interval(str(ch))])
        #    ld = hl.ld_matrix(gt_chr.dosage, gt_chr.locus, radius=1e6)
        #    print('shape: '+str(ld.shape))
        #    ld_sparse = ld.sparsify_band(lower=0, upper = 1000)
        #    ld_sparse.write('gs://nbaya/sumstats_corr/hm3_dosage_correlation_chr'+str(ch)+'.bm',overwrite=True)
        #    print('finished chr '+str(ch))
        #    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now())) 
        #
        ##for ch in range(1,23):
        ##    print('starting chr '+str(ch))
        ##    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now()))
        ##    chr_mt = hl.filter_intervals(gt,[hl.parse_locus_interval(str(ch))])
        ##    chr_bm = BlockMatrix.from_entry_expr(chr_mt.dosage)
        ##    chr_bm = chr_bm.T
        ##    M = chr_bm.shape[1]
        ##    print('chr '+str(ch)+': '+str(M))
        ##    X = chr_bm
        ##    CX = X-X.sum(axis=0)/X.shape[0] #mean-centered
        ##    cov = (1/X.shape[0])*CX.T@CX
        ##    D = (cov**(-1/2)).sparsify_band(lower=0,upper=0)
        ##    X_s = CX@D
        ##    R = (1/X.shape[0])*X_s.T@X_s
        ##    R.sparsify_band(lower=1,upper=1000).write('gs://nbaya/sumstats_corr/hm3_correlation_chr'+str(ch)+'.bm',overwrite=True)
        ##    print('finished chr '+str(ch))
        ##    print('Time: {:%H:%M:%S (%Y-%b-%d)}'.format(datetime.datetime.now())) 



