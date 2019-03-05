import sys
import os
import re
import time
import numpy as np
import Bio.SeqIO as seqio
#import seqio
from optparse import OptionParser,OptionGroup
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParam
from Bio.SeqUtils import GC123

#return_orf,start_coordinate,strand_direction,orf_integrity,orf_len, orf_frame_score
class FindCDS:
    '''
    Find the most like CDS in a given sequence 
    The most like CDS is the longest ORF found in the sequence
    When having same length, the upstream ORF is printed
    modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar
    '''
    def __init__(self,seq):
        self.seq = seq
        self.result = (0,0,0,0,0)
        self.longest = 0
        self.basepair = {"A":"T","T":"A","U":"A","C":"G","G":"C","N":"N","X":"X"}

    def _reversecompliment(self):
        return "".join(self.basepair[base] for base in self.seq)[::-1]

    def get_codons(self,frame_number):
        '''
        Record every nucleotide triplet and its coordinate position for input sequence in one frame
        '''
        coordinate = frame_number
        while coordinate + 3 <= len(self.seq):
            yield (self.seq[coordinate:coordinate+3], coordinate)
            coordinate += 3 
    
    def find_longest_in_one(self,myframe,direction,start_codon,stop_codon):
        '''
        find the longest ORF in one reading myframe
        '''
        triplet_got = self.get_codons(myframe)
        starts = start_codon
        stops = stop_codon
        '''
        Extend sequence by triplet after start codon encountered
        End ORF extension when stop codon encountered
        '''
        while True:
            try: 
                codon,index = triplet_got.next()
            except StopIteration:
                break 
            if codon in starts and codon not in stops:
                '''
                find the ORF start
                '''
                orf_start = index
                end_extension = False
                while True:
                    try: 
                        codon,index = triplet_got.next()
                    except StopIteration:
                        end_extension = True
                        integrity = -1
                    if codon in stops:
                        integrity = 1
                        end_extension = True
                    if end_extension:
                        orf_end = index + 3
                        Length = (orf_end - orf_start)
                        if Length > self.longest:
                            self.longest = Length
                            self.result = [direction,orf_start,orf_end,Length,integrity]
                        if Length == self.longest and orf_start < self.result[1]:
                            '''
                            if ORFs have same length, return the one that if upstream
                            '''
                            self.result = [direction,orf_start,orf_end,Length,integrity]
                        break

    def longest_orf(self,direction,start_codon={"ATG":None}, stop_codon={"TAG":None,"TAA":None,"TGA":None}):
        return_orf = ""
        #orf_frame_score = orf_frame_score(self,direction,start_codon={"ATG":None}, stop_codon={"TAG":None,"TAA":None,"TGA":None})
        for frame in range(3):
            self.find_longest_in_one(frame,"+",start_codon,stop_codon)
        return_orf = self.seq[self.result[1]:self.result[2]][:]
        start_coordinate = self.result[1]
        strand_direction = "+"
        orf_integrity = self.result[4]
        #orf_len = len(return_orf)
        '''
        Also check reverse chain if -r is chosen
        '''
        if direction == "-":
            self.seq = self._reversecompliment()
            for frame in range(3):
                self.find_longest_in_one(frame,"-",start_codon,stop_codon)
            if self.result[0] == "-":
                return_orf = self.seq[self.result[1]:self.result[2]][:]
                start_coordinate = self.result[1]
                strand_direction = "-"
                orf_integrity = self.result[4]
                #orf_len = len(return_orf)
        
        return return_orf,start_coordinate,strand_direction,orf_integrity


#fickett_score
class Fickett:
    '''
    calculate Fickett TESTCODE for full sequence
    NAR 10(17) 5303-531
    modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar 
    '''
    def __init__(self):
        '''new compiled Fickett look-up table'''
        self.position_parameter  = [1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,0.0]
        self.content_parameter  = [0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.19,0.17,0]
        self.position_probability = {
            "A":[0.51,0.55,0.57,0.52,0.48,0.58,0.57,0.54,0.50,0.36],
            "C":[0.29,0.44,0.55,0.49,0.52,0.60,0.60,0.56,0.51,0.38],
            "G":[0.62,0.67,0.74,0.65,0.61,0.62,0.52,0.41,0.31,0.17],
            "T":[0.51,0.60,0.69,0.64,0.62,0.67,0.58,0.48,0.39,0.24],
            }
        self.position_weight = {"A":0.062,"C":0.093,"G":0.205,"T":0.154}
        self.content_probability = {
            "A":[0.40,0.55,0.58,0.58,0.52,0.48,0.45,0.45,0.38,0.19],
            "C":[0.50,0.63,0.59,0.50,0.46,0.45,0.47,0.56,0.59,0.33],
            "G":[0.21,0.40,0.47,0.50,0.52,0.56,0.57,0.52,0.44,0.23],
            "T":[0.30,0.49,0.56,0.53,0.48,0.48,0.52,0.57,0.60,0.51]
            }
        self.content_weight = {"A":0.084,"C":0.076,"G":0.081,"T":0.055}


    def look_up_position_probability(self,value, base):
        '''
        look up positional probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx,val in enumerate (self.position_parameter):
            if (float(value) >= val):
                return float(self.position_probability[base][idx]) * float(self.position_weight[base])

    def look_up_content_probability(self,value, base):
        '''
        look up content probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx,val in enumerate (self.content_parameter):
            if (float(value) >= val):
                return float(self.content_probability[base][idx]) * float(self.content_weight[base])

    def fickett_value(self,dna):
        '''
        calculate Fickett value from full RNA transcript sequence
        '''
        if len(dna) < 2:
            return 0
        fickett_score=0
        dna=dna
        total_base = len(dna)
        A_content = float(dna.count("A"))/total_base
        C_content = float(dna.count("C"))/total_base
        G_content = float(dna.count("G"))/total_base
        T_content = float(dna.count("T"))/total_base

        phase_0 = dna[::3]
        phase_1 = dna[1::3]
        phase_2 = dna[2::3]

        phase_0_A = phase_0.count("A")
        phase_1_A = phase_1.count("A")
        phase_2_A = phase_2.count("A")
        phase_0_C = phase_0.count("C")
        phase_1_C = phase_1.count("C")
        phase_2_C = phase_2.count("C")
        phase_0_G = phase_0.count("G")
        phase_1_G = phase_1.count("G")
        phase_2_G = phase_2.count("G")
        phase_0_T = phase_0.count("T")
        phase_1_T = phase_1.count("T")
        phase_2_T = phase_2.count("T")

        A_content = float(phase_0_A + phase_1_A + phase_2_A)/total_base
        C_content = float(phase_0_C + phase_1_C + phase_2_C)/total_base
        G_content = float(phase_0_G + phase_1_G + phase_2_G)/total_base
        T_content = float(phase_0_T + phase_1_T + phase_2_T)/total_base
        A_position= np.max([phase_0_A,phase_1_A,phase_2_A])/(np.min([phase_0_A,phase_1_A,phase_2_A]) +1.0)
        C_position= np.max([phase_0_C,phase_1_C,phase_2_C])/(np.min([phase_0_C,phase_1_C,phase_2_C]) +1.0)
        G_position= np.max([phase_0_G,phase_1_G,phase_2_G])/(np.min([phase_0_G,phase_1_G,phase_2_G]) +1.0)
        T_position= np.max([phase_0_T,phase_1_T,phase_2_T])/(np.min([phase_0_T,phase_1_T,phase_2_T]) +1.0)

        fickett_score += self.look_up_content_probability(A_content,"A")
        fickett_score += self.look_up_content_probability(C_content,"C")
        fickett_score += self.look_up_content_probability(G_content,"G")
        fickett_score += self.look_up_content_probability(T_content,"T")

        fickett_score += self.look_up_position_probability(A_position,"A")
        fickett_score += self.look_up_position_probability(C_position,"C")
        fickett_score += self.look_up_position_probability(G_position,"G")
        fickett_score += self.look_up_position_probability(T_position,"T")

        return fickett_score

#from Bio.Seq import translate
def mRNA_translate(mRNA):
    return Seq(mRNA).translate()
    #return translate(mRNA)

def GCs(seqRNA):
    gc, gc1, gc2, gc3 = GC123(seqRNA)
    return gc1/100, gc2/100, gc3/100

#stop_codon_count
def stop_codon_func(translate_prot):
    stop_num = translate_prot.count("*")
    return stop_num

#orf_coverage
def orf_coverage(orf_len, trans_len):
    orf_coverage = float(orf_len)/trans_len
    return orf_coverage

#isoelectric_point
def protein_param(putative_seqprot):
    pi = putative_seqprot.isoelectric_point()
    return pi

#Mw
#translate_seqprot = mRNA_translate(seqRNA)
#from Bio.SeqUtils import molecular_weight
def Mw(translated_seqprot):
    if(len(translated_seqprot)==0):
        mw = 0.0
        return mw
    else:
        
        #translated_seqprot = translated_seqprot.strip("*")
        translate_seqprot = ProtParam.ProteinAnalysis(str(translated_seqprot.strip("*")))

        mw = translate_seqprot.molecular_weight()
        return mw 


def pi_mw_ratio(pi, mw):
    if(pi == 0):
        pi_mw = 0.0
        return pi_mw
    else:
        pi_mw = np.log10((float(mw)/pi) + 1)
        return pi_mw

def gravy(translated_seqprot):
    if(len(translated_seqprot)==0):
        Gravy = 4.0
        return Gravy
    else:

        #print("translated seqprot : " , translated_seqprot)
        translated_seqprot = ProtParam.ProteinAnalysis(str(translated_seqprot.strip("*")))
        #print("* baad : ", translated_seqprot)
        #instab = translated_seqprot.instability_index()
        Gravy = translated_seqprot.gravy()
        return Gravy

def instb_ind(translated_seqprot):
    if(len(translated_seqprot)==0):
        instab = 150.0
        return instab
    #print("translated seqprot : " , translated_seqprot)
    else:
        #print("len : ", len(translated_seqprot))
        #translated_seqprot = ProtParam.ProteinAnalysis(str(translated_seqprot))
        translated_seqprot = ProtParam.ProteinAnalysis(str(translated_seqprot.strip("*")))
        #print("trans_seq_prot : ", translated_seqprot)

        #print("* baad : ", translated_seqprot)
        instab = translated_seqprot.instability_index()
        #print("instab : ", instab)
        return instab

def __main():
    start_time = time.time()
    usage = "usage: %prog [options] -i input.fasta -o output_file"
    description = "Contact: Kang Yujian <kangyj@mail.cbi.pku.edu.cn>"
    parser = OptionParser(usage,version="%prog 0.1",description = description)
    Common_group = OptionGroup(parser,"Common Options")
    Common_group.add_option("-i",dest="fasta",help="input sequence in fasta format [Required]",metavar="FILE",type="string",default=None)
    Common_group.add_option("-o",dest="outfile",help="output file [Default: cpc2output.txt]",metavar="FILE",type="string",default="cpc2output.txt")
    Common_group.add_option("-r",dest="reverse",help="also check the reverse strand [Default: FALSE]",action="store_true")
    parser.add_option_group(Common_group)
    (options, args) = parser.parse_args()
    if options.fasta == None:
        parser.print_help()
        return -1
    else:
        if not os.path.isfile(options.fasta):
            sys.stderr.write("[ERROR] %s is not a file\n"%options.fasta)
            return -1
    if options.reverse:
        strand = "-"
    else:
        strand = "+"
    if make_proc_file(options.fasta,strand,options.outfile):
        return 1
    sys.stderr.write("[INFO] cost time: %ds\n"%(time.time()-start_time))
    return 0

def make_proc_file(fasta,strand,outfile):
    '''
    Calculate three features: putative peptide length,pI and Fickett
    And assess coding potential based on SVM model
    '''
    strinfoAmbiguous = re.compile("X|B|Z|J|U",re.I)
    ptU = re.compile("U",re.I)
    ftmp_feat = file(outfile + ".feat","w")
    ftmp_feat.write("\t".join(map(str,["#ID","transcript_length","peptide_length","Fickett_score","ORF_integrity","ORF_coverage", 
                                         "GC1", "GC2", "GC3","stop_codon_num","instability", "gravy", "pI", "Mw", "PW", "#label"]))+"\n")
    fickett_obj = Fickett()
    for seq in seqio.parse(fasta, "fasta"):
    #for seq in seqio.fasta_read(fasta):
        #id
        seqid = seq.id
        seqRNA = ptU.sub("T",str(seq.seq).strip())
        '''seqRNA:transcript full sequence'''
        seqRNA = seqRNA.upper()
        seqCDS,start_pos,orf_strand,orf_fullness = FindCDS(seqRNA).longest_orf(strand)
        seq_gc1, seq_gc2, seq_gc3 = GCs(seq.seq)
        seq_orf_coverage = orf_coverage(len(seqCDS), len(seq))
        
        '''seqCDS:longest ORF'''
        seq_dna_prot = mRNA_translate((seqRNA))
        #print("********")
        #print("seqCDS : ", seqCDS)
        seqprot = mRNA_translate(seqCDS)
        #print("seqprot : ", seqprot)
        #print("********")
        #print("....................")
        #print("longest orf translated_protein : ", seqprot)
        seq_prot_stop_codon_num = stop_codon_func(seq_dna_prot)
        #print("stop_codon_num : ", seq_prot_stop_codon_num)
        #print("....................")
        print("seqid : ", seqid)
        seq_prot_mw = Mw(seqprot)
        seq_prot_instb = instb_ind(seqprot)
        seq_prot_gravy = gravy(seqprot)
        
        pep_len = len(seqprot) #pep_len = len(seqprot.strip("*"))
        newseqprot = strinfoAmbiguous.sub("",str(seqprot))
        '''exclude ambiguous amino acid X, B, Z, J, Y in peptide sequence'''
        
        fickett_score = fickett_obj.fickett_value(seqRNA)
        
        protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot.strip("*")))
        
        #iso = protparam_obj.isoelectric_point()
        
        if pep_len > 0:
            #fickett_score = fickett_obj.fickett_value(seqCDS)
            isoelectric_point = protein_param(protparam_obj)
            #isoelectric_point = iso

				        	
        else:
            #fickett_score = 0.0
            orf_fullness = -1
            isoelectric_point = 0.0
        seq_prot_pi_mw_ratio = pi_mw_ratio(isoelectric_point, seq_prot_mw) 
        
        ftmp_feat.write("\t".join(map(str,[seqid, len(seqRNA), pep_len, fickett_score, orf_fullness, seq_orf_coverage, seq_gc1, seq_gc2, seq_gc3, seq_prot_stop_codon_num, seq_prot_instb, seq_prot_gravy, isoelectric_point, seq_prot_mw, seq_prot_pi_mw_ratio, "noncoding"]))+"\n")
        #ftmp_feat.write("\t".join(map(str,[seqid, len(seqRNA), pep_len, fickett_score, orf_fullness, seq_orf_coverage, seq_gc1, seq_gc2, seq_gc3, seq_prot_stop_codon_num, seq_prot_instb, seq_prot_gravy, isoelectric_point, seq_prot_mw, seq_prot_pi_mw_ratio, "coding"]))+"\n")   
        #ftmp_svm.write("".join(map(str,["999"," 1:",pep_len," 2:",fickett_score," 3:",isoelectric_point," 4:",orf_fullness]))+"\n")
    ftmp_feat.close()
    
if __name__ == "__main__":
    sys.exit(__main())

