import numpy as np
import pandas as pd
import argparse
from itertools import combinations
from itertools import product
from time import sleep

class DNA:
    '''
    A = np.matrix([[0,0],[0,0]])
    T = np.matrix([[1,0],[0,1]])
    C = np.matrix([[0,1],[1,1]])
    G = np.matrix([[1,1],[1,0]])
    '''
    # Generator_matrix=[
        #                 [ 'T', 'A','A','A','A','A','T','T','G','C','T','T'],
        #                 [ 'A', 'T','A','A','A','A','C','G','T','T','T','T'],
        #                 [ 'A', 'A','T','A','A','A','G','C','G','A','T','C'],
        #                 [ 'A', 'A','A','T','A','A','C','A','C','G','T','G'],
        #                 [ 'A', 'A','A','A','T','A','G','C','A','C','G','T'],
        #                 [ 'A', 'A','A','A','A','T','A','G','C','G','C','T']]
    Gen_mat_12_6_6=[
                        [ 1,0 , 0,0 ,0,0 ,0,0 ,0,0 ,0,0 ,1,0 ,1,0 ,1,1 ,0,1 ,1,0 ,1,0 ],
                        [ 0,1 , 0,0 ,0,0 ,0,0 ,0,0 ,0,0 ,0,1 ,0,1 ,1,0 ,1,1 ,0,1 ,0,1 ],
                        [ 0,0 , 1,0 ,0,0 ,0,0 ,0,0 ,0,0 ,0,1 ,1,1 ,1,0 ,1,0 ,1,0 ,1,0 ],
                        [ 0,0 , 0,1 ,0,0 ,0,0 ,0,0 ,0,0 ,1,1 ,1,0 ,0,1 ,0,1 ,0,1 ,0,1 ],
                        [ 0,0 , 0,0 ,1,0 ,0,0 ,0,0 ,0,0 ,1,1 ,0,1 ,1,1 ,0,0 ,1,0 ,0,1 ],
                        [ 0,0 , 0,0 ,0,1 ,0,0 ,0,0 ,0,0 ,1,0 ,1,1 ,1,0 ,0,0 ,0,1 ,1,1 ],
                        [ 0,0 , 0,0 ,0,0 ,1,0 ,0,0 ,0,0 ,0,1 ,0,0 ,0,1 ,1,1 ,1,0 ,1,1 ],
                        [ 0,0 , 0,0 ,0,0 ,0,1 ,0,0 ,0,0 ,1,1 ,0,0 ,1,1 ,1,0 ,0,1 ,1,0 ],
                        [ 0,0 , 0,0 ,0,0 ,0,0 ,1,0 ,0,0 ,1,1 ,0,1 ,0,0 ,0,1 ,1,1 ,1,0 ],
                        [ 0,0 , 0,0 ,0,0 ,0,0 ,0,1 ,0,0 ,1,0 ,1,1 ,0,0 ,1,1 ,1,0 ,0,1 ],
                        [ 0,0 , 0,0 ,0,0 ,0,0 ,0,0 ,1,0 ,0,0 ,1,1 ,0,1 ,1,1 ,0,1 ,1,0 ],
                        [ 0,0 , 0,0 ,0,0 ,0,0 ,0,0 ,0,1 ,0,0 ,1,0 ,1,1 ,1,0 ,1,1 ,0,1 ]]
                        
    # Generator_matrix=[
    #                 [ 'T', 'A','A','C','T','C'],
    #                 [ 'A', 'T','A','G','G','T'],
    #                 [ 'A', 'A','T','A','G','C'],
    
    Gen_mat_6_3_3=[
                        [ 1,0 , 0,0 ,0,0 ,0,1 ,1,0 ,0,1 ],
                        [ 0,1 , 0,0 ,0,0 ,1,1 ,0,1 ,1,1 ],
                        [ 0,0 , 1,0 ,0,0 ,1,1 ,1,1 ,1,0 ],
                        [ 0,0 , 0,1 ,0,0 ,1,0 ,1,0 ,0,1 ],
                        [ 0,0 , 0,0 ,1,0 ,0,0 ,1,1 ,0,1 ],
                        [ 0,0 , 0,0 ,0,1 ,0,0 ,1,0 ,1,1 ]]
    
        # Generator_matrix=[
    #                 [ 'T', 'A','C','G'],
    #                 [ 'A', 'T','G','C'],
    
    Gen_mat_4_2_3=[
                        [ 1,0 , 0,0 ,1,1 ,0,1],
                        [ 0,1 , 0,0 ,1,0 ,1,1],
                        [ 0,0 , 1,0 ,0,1 ,1,1],
                        [ 0,0 , 0,1 ,1,1 ,1,0]]
    
    Generator_matrix=Gen_mat_4_2_3

    header_length = 15
    Generator_matrix=np.matrix(Generator_matrix)    
    dna_code_length = int((Generator_matrix.shape[1])/2)
    
    bin_data = ""
    data_byte = "" # ascii code "##" -> 010101011
    matrix = np.empty([2,2])
    matrix_bin = np.empty([2,2])
    original_text = ""
    dict_bin = {'00' : 'A',
            '11' : 'T',
            '10' : 'C',
            '01' : 'G'
        }
    dict_DNA = {'A' : np.matrix([[0,0],[0,0]]),
            'T' : np.matrix([[1,0],[0,1]]),
            'C' : np.matrix([[0,1],[1,1]]),
            'G' : np.matrix([[1,1],[1,0]])
    }
    coset_leaders_dict=dict()
    


    def __init__(self) :
        None
    def readfile(self, args_file) :
        with open(args_file, 'r') as f:
                    while True :
                        line = f.readline()
                        for word in line:
                            self.bin_data = self.bin_data + bin(ord(word)).replace("0b","").zfill(8) #zfill을 이용해서 무조건 8자리로 만듬
                        if not line:
                            break
        with open(str(args_file)+"bin_data","w") as g:
            g.write(self.bin_data)
        g.close()
        return self.bin_data

    def DNA_to_matrix(self,DNA_seq) :
        matrix_ = np.empty([2,2])
        for nucleotide in DNA_seq :
            matrix_ = np.concatenate((matrix_,self.dict_DNA[nucleotide]), axis =1 )
        matrix_=matrix_[0:,2:]
        return matrix_

    def matrix_to_DNA(self, mat) :
        DNA_seq = ""
        iter_nums = int((mat.shape[1])/2)
        j=0
        for n in range(iter_nums) :
            DNA_seq += str(self.return_key_DNA(mat[:,j:j+2]))
            j+=2
        return DNA_seq

    def bin_to_DNA_seq(self) :
        DNA_seq = ""
        for a,b in zip(self.bin_data[0::2], self.bin_data[1::2]) :
            c = a+b
            DNA_seq += self.dict_bin[c]
        return DNA_seq

    def header(self):
        header_str = bin(len(self.bin_data)).replace("0b","").zfill(self.header_length)
        num_zeros = self.dna_code_length*2-(len(self.bin_data)+self.header_length)%(self.dna_code_length*2)
        for i in range(num_zeros):
            header_str +='0'
        self.bin_data = header_str+self.bin_data
        return self.bin_data, header_str

    def matrix_multiply(self,Mat) :
        result = Mat*self.Generator_matrix
        result = result%2
        return result

    def encoding_to_mat_seq(self, Mat) :
        """
        self.dna_code_length = 12
        self.header_length  = 20
        self.matrix
        """
        block_count = Mat.shape[1]
        for i in range(int(block_count/self.dna_code_length)) :
            sub_matrix = self.matrix_multiply(Mat[0:,i*self.dna_code_length:(i+1)*self.dna_code_length])
            self.matrix_bin = np.concatenate((self.matrix_bin, sub_matrix), axis = 1 )
        self.matrix_bin = self.matrix_bin[0:,2:].astype(int)
        return self.matrix_bin
        
    def Codeword_seq_to_vec_seq(self, codeword_seq) :
        cnt_vec = int(len(codeword_seq)/self.dna_code_length)
        if len(codeword_seq)%self.dna_code_length !=0 :
            print("Length of Codeword sequence is not multiple of {}".format(self.dna_code_length))
            return
        vec_seq = ""
        for i in range(cnt_vec) :
            codeword = codeword_seq[(i)*(self.dna_code_length):(i+1)*(self.dna_code_length)]
            codeword = self.DNA_to_matrix(codeword)
            error_corrected_codeword=self.Syndrome_decoding(codeword,i)
            error_corrected_codeword=self.matrix_to_DNA(error_corrected_codeword)
            vec = error_corrected_codeword[0:int(self.dna_code_length/2)]
            vec_seq += vec
        return vec_seq
    

    def coset_leaders(self) :
        identity_mat=np.identity(self.dna_code_length*2,dtype=int)  #binary identiti matrix
        canonical_bases=[identity_mat[2*i:2*(i+1),:self.dna_code_length*2] for i in range(self.dna_code_length)]  

        #non zero scalars
        T_=DNA.dict_DNA['T']
        C_=DNA.dict_DNA['C']
        G_=DNA.dict_DNA['G']
        scalars=[T_,C_,G_]

        bases=[]  
        rows=(self.Generator_matrix.T).shape[0] 
        cols=(self.Generator_matrix.T).shape[1] 
        cnt_basis=int(rows/2)

        item=list(range(cnt_basis))

        for i in range(cnt_basis):
            bases.append((self.Generator_matrix.T)[2*i:2*(i+1),:cols])

        for k in range(2):
            comb=combinations(item, k+1)
            coeff=product(scalars, repeat=k+1)
            comb_coeff=product(comb, coeff)
 
            for tup in comb_coeff:
                vecs=[tup[1][i]*bases[tup[0][i]] for i in range(k+1)]
                error_vecs=[tup[1][i]*canonical_bases[tup[0][i]] for i in range(k+1)]
                coset_leader=sum(vecs)%2
                errors=sum(error_vecs)%2
                coset_leader_DNA=self.matrix_to_DNA(coset_leader)
                error_vec_DNA=self.matrix_to_DNA(errors)
                if coset_leader_DNA in self.coset_leaders_dict.keys():
                    if error_vec_DNA not in self.coset_leaders_dict[coset_leader_DNA]:
                        self.coset_leaders_dict[coset_leader_DNA].append(error_vec_DNA)
                else:
                    self.coset_leaders_dict.update({coset_leader_DNA:[error_vec_DNA]})



    def Syndrome_decoding(self, codeword,i_th) :
        syndrome=(codeword*(self.Generator_matrix.T))%2
        #print(syndrome.shape)
        is_error=False
        if np.array_equal(syndrome, np.matrix(np.zeros((2,self.dna_code_length)))):
            decoded=codeword
        else:
            is_error=True 
            syndrome_DNA = self.matrix_to_DNA(syndrome)
            if syndrome_DNA in self.coset_leaders_dict.keys():
                print(f"The {i_th}-th codeword is corrected")
                error_vec_DNA = self.coset_leaders_dict[syndrome_DNA][0]
                decoded = codeword+self.DNA_to_matrix(error_vec_DNA)
                decoded=decoded%2
            else:
                print(f"The {i_th}-th codeword has too many errors to be corrected")
                decoded=codeword

        return decoded

    def DNA_seq_to_bin(self, DNA) :
        bin_seq = ""
        for i in DNA :
            bin_seq += self.return_key_bin(i)
        return bin_seq

    def bin_to_txt(self, bin_seq) :
        txt=""
        header_vec = bin_seq[0:self.header_length]
        len_binary_seq = int(header_vec,2)
        bin_txt=bin_seq[-len_binary_seq:]
        for i in range(len(bin_txt)//8):
            ord_bin = bin_txt[i*8:(i+1)*8]
            txt += chr(int(str(ord_bin),2))
        return txt
        
    def return_key_bin(self, val):
        for key, value in self.dict_bin.items():
            if value==val :
                return key
    def return_key_DNA(self, val):
        for key, value in self.dict_DNA.items():
            if (value==val).all():
                return key


    def indel_correcting_encoding(self, codeword_seq) :
        cnt_vec = int(len(codeword_seq)/self.dna_code_length)
        poly_T_seq_mat = self.DNA_to_matrix("".join(['T' for i in range(self.dna_code_length+2)]))

        if len(codeword_seq)%self.dna_code_length !=0 :
            print("Length of Codeword sequence is not multiple of {}".format(self.dna_code_length))
            return
        indel_seq = ""
        for i in range(cnt_vec) :
            codeword = codeword_seq[(i)*(self.dna_code_length):(i+1)*(self.dna_code_length)]
            if codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='AA':
                codeword=codeword[:int(self.dna_code_length/2)]+'CC'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='AT':
                codeword=codeword[:int(self.dna_code_length/2)]+'GG'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='AC':
                codeword=codeword[:int(self.dna_code_length/2)]+'TT'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='AG':
                codeword=codeword[:int(self.dna_code_length/2)]+'TT'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='TA':
                codeword=codeword[:int(self.dna_code_length/2)]+'GG'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='TT':
                codeword=codeword[:int(self.dna_code_length/2)]+'CC'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='TC':
                codeword=codeword[:int(self.dna_code_length/2)]+'GG'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='TG':
                codeword=codeword[:int(self.dna_code_length/2)]+'CC'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='CA':
                codeword=codeword[:int(self.dna_code_length/2)]+'TT'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='CT':
                codeword=codeword[:int(self.dna_code_length/2)]+'GG'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='CC':
                codeword=codeword[:int(self.dna_code_length/2)]+'AA'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='CG':
                codeword=codeword[:int(self.dna_code_length/2)]+'TT'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='GA':
                codeword=codeword[:int(self.dna_code_length/2)]+'TT'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='GT':
                codeword=codeword[:int(self.dna_code_length/2)]+'CC'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='GC':
                codeword=codeword[:int(self.dna_code_length/2)]+'TT'+codeword[int(self.dna_code_length/2):]
            elif codeword[int(self.dna_code_length/2)-1:int(self.dna_code_length/2)+1]=='GG':
                codeword=codeword[:int(self.dna_code_length/2)]+'AA'+codeword[int(self.dna_code_length/2):] 
            if i>0 and indel_seq[-1] == codeword[0]:
                complementary_codeword=self.matrix_to_DNA((self.DNA_to_matrix(codeword)+poly_T_seq_mat)%2)
                codeword=complementary_codeword
            indel_seq += codeword
        return indel_seq


    def is_complementary(self, x,y,c_psi):
        if x=='A' and y=='A':
            if c_psi =="C":
                complementary_check=False
            elif c_psi =="G":
                complementary_check=True
            else:
                complementary_check="Error"                                
        elif x=='A' and y=='T':
            if c_psi =="G":
                complementary_check=False
            elif c_psi =="C":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='A' and y=='C':
            if c_psi =="T":
                complementary_check=False
            elif c_psi =="G":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='A' and y=='G':
            if c_psi =="T":
                complementary_check=False
            elif c_psi =="C":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='T' and y=='A':
            if c_psi =="G":
                complementary_check=False
            elif c_psi =="C":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='T' and y=='T':
            if c_psi =="C":
                complementary_check=False
            elif c_psi =="G":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='T' and y=='C':
            if c_psi =="G":
                complementary_check=False
            elif c_psi =="A":
                complementary_check=True
            else:
                complementary_check="Error"        
        elif x=='T' and y=='G':
            if c_psi =="C":
                complementary_check=False
            elif c_psi =="A":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='C' and y=='A':
            if c_psi =="T":
                complementary_check=False
            elif c_psi =="G":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='C' and y=='T':
            if c_psi =="G":
                complementary_check=False
            elif c_psi =="A":
                complementary_check=True
            else:
                complementary_check="Error"        
        elif x=='C' and y=='C':
            if c_psi =="A":
                complementary_check=False
            elif c_psi =="T":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='C' and y=='G':
            if c_psi =="T":
                complementary_check=False
            elif c_psi =="A":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='G' and y=='A':
            if c_psi =="T":
                complementary_check=False
            elif c_psi =="C":
                complementary_check=True
            else:
                complementary_check="Error"        
        elif x=='G' and y=='T':
            if c_psi =="C":
                complementary_check=False
            elif c_psi =="A":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='G' and y=='C':
            if c_psi =="T":
                complementary_check=False
            elif c_psi =="A":
                complementary_check=True
            else:
                complementary_check="Error"
        elif x=='G' and y=='G':
            if c_psi =="A":
                complementary_check=False
            elif c_psi =="T":
                complementary_check=True
            else:
                complementary_check="Error"
        return complementary_check

    def deletion_correction(self, codeword_seq) :
        seq_len = len(codeword_seq)
        indel_corrected_seq = ""
        seq_pos=0
        block_cnt=0
        poly_T_seq_mat = self.DNA_to_matrix("".join(['T' for i in range(self.dna_code_length)]))
        while seq_pos <= seq_len-int(self.dna_code_length):
            if seq_len-seq_pos==int(self.dna_code_length)+2: 
                block_codeword=codeword_seq[seq_pos:]                
                block_cnt+=1
                x=block_codeword[int(self.dna_code_length/2)-1] ;y=block_codeword[int(self.dna_code_length/2)+2] ; c_psi = block_codeword[int(self.dna_code_length/2)]
                information_block = block_codeword[:int(self.dna_code_length/2)]
                if self.is_complementary(x,y,c_psi)==True:
                    print(information_block)
                    print(poly_T_seq_mat)
                    indel_corrected_block=self.matrix_to_DNA((self.matrix_multiply(self.DNA_to_matrix(information_block))+poly_T_seq_mat)%2)
                else:
                    indel_corrected_block=self.matrix_to_DNA(self.matrix_multiply(self.DNA_to_matrix(information_block)))
                indel_corrected_seq += indel_corrected_block[:int(self.dna_code_length)]
                Lprint = str(block_cnt) + " block " + str(block_codeword) + ' has no deletion'
                print(Lprint, ', decoded as ', indel_corrected_block[:int(self.dna_code_length)])
                break
            elif seq_len-seq_pos==int(self.dna_code_length)+1:
                block_codeword=codeword_seq[seq_pos:]
                block_cnt+=1
                Lprint=str(block_cnt)+" block "+str(block_codeword)+' has a deletion'
                c_psi = block_codeword[int(self.dna_code_length/2)]
                if block_codeword[int(self.dna_code_length/2)]==block_codeword[int(self.dna_code_length/2)+1]: 
                    information_block = block_codeword[:int(self.dna_code_length/2)]
                    indel_corrected_block=self.matrix_to_DNA(self.matrix_multiply(self.DNA_to_matrix(information_block)))
                    x=indel_corrected_block[int(self.dna_code_length/2)-1] ;y=indel_corrected_block[int(self.dna_code_length/2)]
                    if self.is_complementary(x,y,c_psi)==True:
                        indel_corrected_block=self.matrix_to_DNA((self.DNA_to_matrix(indel_corrected_block)+poly_T_seq_mat)%2)
                elif block_codeword[int(self.dna_code_length/2)-1]==block_codeword[int(self.dna_code_length/2)]:
                    block_codeword=block_codeword[int(self.dna_code_length/2)+1:int(self.dna_code_length)+1]
                    block_codeword=block_codeword[::-1]
                    indel_corrected_block=self.matrix_to_DNA(self.matrix_multiply(self.DNA_to_matrix(block_codeword)))[::-1]
                    x=indel_corrected_block[int(self.dna_code_length/2)-1] ;y=indel_corrected_block[int(self.dna_code_length/2)]
                    if self.is_complementary(x,y,c_psi)==True:
                        indel_corrected_block=self.matrix_to_DNA((self.DNA_to_matrix(indel_corrected_block)+poly_T_seq_mat)%2)
                else: 
                    block_codeword=block_codeword[:int(self.dna_code_length/2)] 
                    indel_corrected_block=self.matrix_to_DNA(self.matrix_multiply(self.DNA_to_matrix(block_codeword)))
                    x=indel_corrected_block[int(self.dna_code_length/2)-1] ;y=indel_corrected_block[int(self.dna_code_length/2)]
                    if self.is_complementary(x,y,c_psi)==True:
                        indel_corrected_block=self.matrix_to_DNA((self.DNA_to_matrix(indel_corrected_block)+poly_T_seq_mat)%2)
                indel_corrected_seq += indel_corrected_block[:int(self.dna_code_length)]
                Lprint = str(block_cnt) + " block " + str(block_codeword) + ' has 1 deletion'
                print(Lprint, ', decoded as ', indel_corrected_block[:int(self.dna_code_length)])
                break
            else: 
                block_codeword=codeword_seq[seq_pos:seq_pos+self.dna_code_length+2] 
                block_cnt+=1
                c_psi = block_codeword[int(self.dna_code_length/2)]
                if block_codeword[int(self.dna_code_length/2)]==block_codeword[int(self.dna_code_length/2)+1]:  
                    information_block = block_codeword[:int(self.dna_code_length/2)]
                    pre_indel_corrected_block=self.matrix_to_DNA(self.matrix_multiply(self.DNA_to_matrix(information_block)))
                    x=pre_indel_corrected_block[int(self.dna_code_length/2)-1] ;y=pre_indel_corrected_block[int(self.dna_code_length/2)]
                    indel_corrected_block=self.matrix_to_DNA(self.DNA_to_matrix(pre_indel_corrected_block))
                    if indel_corrected_block[int(self.dna_code_length/2):]==block_codeword[int(self.dna_code_length/2)+2:] :
                        seq_pos+=int(self.dna_code_length)+2 
                        if self.is_complementary(x,y,c_psi)==True :
                            indel_corrected_block=self.matrix_to_DNA((self.DNA_to_matrix(pre_indel_corrected_block)+poly_T_seq_mat)%2)

                    else:
                        seq_pos+=int(self.dna_code_length)+1 
                        if self.is_complementary(x,y,c_psi)==True :
                            indel_corrected_block=self.matrix_to_DNA((self.DNA_to_matrix(pre_indel_corrected_block)+poly_T_seq_mat)%2)
                        Lprint=str(block_cnt)+" block "+str(block_codeword)+' has a deletion in back'
                        print(Lprint,', decoded as ', indel_corrected_block[:int(self.dna_code_length)])
                else: 
                    Lprint=str(block_cnt)+" block "+str(block_codeword)+' has a deletion in front' 
                    block_codeword=block_codeword[int(self.dna_code_length/2)+1:int(self.dna_code_length)+1]
                    block_codeword=block_codeword[::-1]
                    indel_corrected_block=self.matrix_to_DNA(self.matrix_multiply(self.DNA_to_matrix(block_codeword)))[::-1]
                    x=indel_corrected_block[int(self.dna_code_length/2)-1] ;y=indel_corrected_block[int(self.dna_code_length/2)]
                    if self.is_complementary(x,y,c_psi)==True:
                        indel_corrected_block=self.matrix_to_DNA((self.DNA_to_matrix(indel_corrected_block)+poly_T_seq_mat)%2)
                    seq_pos+=int(self.dna_code_length)+1 
                    print(Lprint,', decoded as ', indel_corrected_block[:int(self.dna_code_length)])
                indel_corrected_seq += indel_corrected_block[:int(self.dna_code_length)]
        return indel_corrected_seq

def Encode(dna, file) :
    print("Encoding Process Starting...")
    result = dna.readfile(file)    
    if input("Show input text? [y/n]") == "y":
        print(dna.original_text)

    bin_data, header_str = dna.header()
    if input("Show header bit information? [y/n]") == "y":
        print(header_str)

    DNA_seq = dna.bin_to_DNA_seq()    
    if input("show DNA Seq? [y/n]") == "y" :
        print(DNA_seq)

    DNA_matrix = dna.DNA_to_matrix(DNA_seq)
    if input("show DNA Matrix? [y/n]") == "y" :
        print(DNA_matrix)
    
    DNA_Generator_multiplied = dna.encoding_to_mat_seq(DNA_matrix)
    DNA_Seq = dna.matrix_to_DNA(DNA_Generator_multiplied)
    print("Generator Matrix multiplied.")

    DNA_Encoded_Seq = dna.indel_correcting_encoding(DNA_Seq)

    with open("DNA_Encoded_"+str(file),"w") as f:
        f.write(DNA_Encoded_Seq)
    f.close()

    if input("Show Encoded DNA Matrix? [y/n]")== "y" :
        print(DNA_Encoded_Seq)
    print(f'length: {len(DNA_Encoded_Seq)}')
        

def Decode(dna, file) :
    print("Decoding Process Starting...")
    
    DNA_Encoded_seq = ""
    with open(file, "r") as f:
        DNA_Encoded_seq = f.read()
    f.close()

    print(f'DNA Encoding Length : {len(DNA_Encoded_seq)}')
    if input("Show Encoded DNA Seq? [y/n]") == "y" :
        print(DNA_Encoded_seq)
    
    print("Coset leader generating.")

    dna.coset_leaders()

    print("Coset leader generated!")

    Del_corrected=dna.deletion_correction(DNA_Encoded_seq) 
    vec_seq = dna.Codeword_seq_to_vec_seq(Del_corrected) # syndrome
    print("Syndrome Decoding.")

    if input("Show Fixed DNA Seq? [y/n]") == "y":
        print(vec_seq)


    bin_seq = dna.DNA_seq_to_bin(vec_seq)
    txt_data = dna.bin_to_txt(bin_seq)
    if input("Show Original txt? [y/n]") == "y":
        print(txt_data) 
    with open("Retrieved_" +str(file), "w") as f:
        f.write(txt_data)
    f.close()


    return
    

def getArgumentParser() :
    parser = argparse.ArgumentParser(description = "Ultimate Matrix Encoder/Decoder for DNA", prog = 'dna')
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-e","--encode", help = "file to be encoded")
    group.add_argument("-d","--decode", help = "file to be decoded")
    args = parser.parse_args()
    return args

def main() :
    args = getArgumentParser()
    dna = DNA()
    
    if args.encode :
        Encoded_text = Encode(dna, args.encode)
    elif args.decode:
        Decoded_text = Decode(dna, args.decode)

if __name__ == "__main__" :
    main()
