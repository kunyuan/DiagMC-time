import numpy as np
import subprocess
import os
import sys
import StringIO


def index(file, separator='###*'):
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

def parse_scope(index_list, length,scope=':'):
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
    index_list, length=index(file, separator)
    scope_bound=parse_scope(index_list, length, scope)
    #print scope_bound
    data_block=[]
    f=open(file,'r')
    #print len(scope_bound)
    for i in range(0,len(scope_bound)):
        f.seek(scope_bound[i][0])
        buff=StringIO.StringIO(f.read(scope_bound[i][1]-scope_bound[i][0]))
        print buff.readline()
        # strbuff=[e.strip() for e in buff.readline()[1:].split(',')]
        strbuff=buff.readline().strip(' \t\n\r')[1:]
        dim=tuple(map(int, strbuff.split(',')))
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
