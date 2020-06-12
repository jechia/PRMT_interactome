#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2016-6-5
E-mail: fanxiaojuan@picb.ac.cn
Description: Count z-score of motif
"""

"""
input file format:
#------------------------------------------------------------------------------ 
motif    count    frequency
GCGTT    292    0.0546714098483
AAATG    1374    0.257255195656
GCCCG    459    0.0859389627411
GCCCA    1019    0.190788241902
AAATC    1061    0.198651937839
"""

"""
formula
#------------------------------------------------------------------------------ 
z_score = (f(H)1-f(H)2)/((N1-N2)*p(1-p))**(1/2)
p = (f(H)1*N1-f(H)2*N2)/(N1+N2)
"""

import argparse
import math
import os

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to count z-score of motif.'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-p', '--input sample 1', dest='p', default='F:/experiment/WeiHuanhuan/XiaoJuan/PRMTs_293T/arginine_frequency/', help='input sample path')
    parser.add_argument('--i-s2', '--input sample 2', dest='fnIn_s2', default='F:/experiment/WeiHuanhuan/XiaoJuan/motifs_batch2/motif_files/uniprot-human-reviewed-sp-20180830-flat_frequency.txt', help='input sample 2 as background')
    op=parser.parse_args()
    return op

if __name__ == '__main__':
    op = createHelp()
    
    print ('Build motif dictionary in background dataset.')
    
    background = open(op.fnIn_s2).readlines()
    bg_dic = {}
    for i in range(0,len(background)):
        words_bg = background[i].strip().split('\t')
        motif_bg = words_bg[0]
        count_bg = int(words_bg[1])
        freq_bg = float(words_bg[2])
        bg_dic[motif_bg] = [count_bg,freq_bg]
        
    #print bg_dic
    print ('Background dic done!')
    
    print ('Start to calculate z-score...')
    for dirNA in os.listdir(op.p):
        if 'mock_ctl_specific_proteins_frequency.txt' in dirNA:
            print (dirNA + ' start...')
            file_prefix = dirNA.split('.')[0]
            freq_file = open(op.p+dirNA).readlines()
            
            final_list = []
            for i in range(0,len(freq_file)):
                words_cir = freq_file[i].strip().split('\t')
                motif_cir = words_cir[0]
                count_cir = int(words_cir[1])
                freq_cir = float(words_cir[2])
                if motif_cir in bg_dic:
                    if freq_cir != 0:
                        p = (count_cir + bg_dic[motif_cir][0])/((count_cir / freq_cir + bg_dic[motif_cir][0] / bg_dic[motif_cir][1]))
                        #print p
                        z_score = (freq_cir - bg_dic[motif_cir][1])/math.sqrt((1/(count_cir/freq_cir) + 1/(bg_dic[motif_cir][0]/bg_dic[motif_cir][1])) * p * (1-p))
                        #print z_score
                        fold_change = freq_cir / bg_dic[motif_cir][1]
                        final_list.append('\t'.join(words_cir)+'\t'+str(bg_dic[motif_cir][0])+'\t'+str(bg_dic[motif_cir][1])+'\t'+str(fold_change)+'\t'+str(z_score))
                else:
                    #print motif_cir
                    final_list.append('\t'.join(words_cir))
                    
            #print final_list
            print ('Done! Start to write to output file...')
            output = open(op.p+file_prefix+'_score.txt','w')
        #     final_list.sort(key = lambda l:(l[1]),reverse=True)
            output.write('\n'.join(final_list)+'\n')
        
    print ('Finish!')
