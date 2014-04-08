import numpy as np
import subprocess
import os
import sys
import StringIO


def __index(file, separator='###*'):
    index_list = []
    if os.path.exists(file):
        result=subprocess.check_output(['grep', '-ob', separator, file])
        for item in result[:-1].split('\n'):
            index_list.append(int(item.split(':')[0]))
        #print "Block Number: {}".format(len(index_list)), index_list
        #result=subprocess.check_output(['wc','-c',file])
        #length=int(result.split(' ')[0])
        length=os.path.getsize(file)
        for i in range(0, len(index_list)-1):
            index_list[i]=(index_list[i], index_list[i+1])
        index_list[-1]=(index_list[-1], length)
        # print index_list

    else:
        print "Data file '{}' does not exist!".format(file)
        sys.exit()

    #index_list.append(os.path.getsize(file))
    #the last item is the end of the file

    return index_list, length

def __parse_scope(index_list, length,scope=':'):
    scope_bound_str=scope.split(':')
    if len(scope_bound_str)>2:
        print "Only 2 parameters are allowed in range: {}".format(scope)
        sys.exit()
    else:
        try:
            if scope_bound_str[0]!='':
                header=int(scope_bound_str[0])
            else:
                header=0

            if scope_bound_str[1]!='':
                scope_bound=index_list[header:int(scope_bound_str[1])]
            else:
                scope_bound=index_list[header:]

        except:
            print "Please check your range: {}".format(scope)
            sys.exit()

    return scope_bound

def __parse_dim(dimstr):
    dim=[]
    dim_name=[]
    strbuff=dimstr.strip(' \t\n\r')[1:].split(',')
    for i in range(0,len(strbuff)):
        elem=strbuff[i].strip(' \t\n\r').split(':')
        if len(elem)==1:
            dim.append(int(elem[0]))
            if len(strbuff)<=3:
                if i==0:
                    dim_name.append("X")
                elif i==1:
                    dim_name.append("Y")
                else:
                    dim_name.append("Z")
            else:
                dim_name.append('X'+str(i-1))

        elif len(elem)==2:
            dim.append(int(elem[1]))
            if elem[0].strip(' \t\n\r')=='':
                dim_name.append('X'+str(i-1))
            else:
                dim_name.append(elem[0])
        else:
            print "Illegal dimensional information: {}".format(dimstr)
    if len(strbuff)==1:
        dim_name[0]="X"
    if len(strbuff)==3:
        dim_name[0]="X"
        dim_name[1]="Y"
        dim_name[2]="Z"

    return dim, dim_name

def read_array(file, scope='-1:', comment='#', separator='###*'):
    '''
    Return: a tuple as (list, list, list)
    ------------------------
        -list: a list contains blocks of data
        -list: a list conatins dimensional information, e.g. [4,3]
        -list: a list contains the name of each dimension, e.g. ['X','Y']

    Parameters:
    -----------
    file : string
        The data file path and name

    scope : string with format 'int:int', default='-1:'
        The blocks range needed to read, as the index in list

    comment : string, default='#'
        The comment characters in the data file

    separator : string, default='###*'
        The string to separate different blocks in the data file
    '''
    index_list, length=__index(file, separator)
    scope_bound=__parse_scope(index_list, length, scope)
    #print scope_bound
    data_block=[]
    f=open(file,'r')
    #print len(scope_bound)
    for i in range(0,len(scope_bound)):
        f.seek(scope_bound[i][0])
        buff=StringIO.StringIO(f.read(scope_bound[i][1]-scope_bound[i][0]))
        buff.readline()
        dim, dim_name=__parse_dim(buff.readline())
        # strbuff=[e.strip() for e in buff.readline()[1:].split(',')]
        a=np.loadtxt(buff,comments=comment)
        if a.shape[1]>2:
            print "At most two numbers in each row!"
            sys.exit()
        elif a.shape[1]==2:
            # print a[0:3,:]
            a.dtype=complex
            #print a[0:3]
        data_block.append(a.reshape(dim,order='F'))
    f.close()
    #print len(data_block)
    if len(dim)==2:
        dim_name=["X","Y"]
    else:
        dim_name=[]
    return data_block, dim, dim_name
    # offset=scope_bound[0]
    # for line in f:
    #if offset<scope_bound[1]:
    #         print line
    #     else:
    #         break
    #     offset+=len(line)
    #a=np.genfromtxt(file,comments='#',skip_header=scope_bound[0],
    #   skip_footer=length-scope_bound[1])
    #print a


if __name__=="__main__":
    read_array('0_0.50_1_Gam_matrix.dat')
