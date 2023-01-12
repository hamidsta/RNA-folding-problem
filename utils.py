import copy
import os
import pandas as pd
import collections
from collections import defaultdict
import math
from pathlib import Path
import sys
import os.path
import sys
from urllib.request import urlopen
import matplotlib.pyplot as plt


def plot_script_2(nested_dict,direct):
    for key, value in nested_dict.items():
        print('plotting:' , key)
        #for k2,v2 in value.items():
        names = list(value.keys())
        values = list(value.values())
        fig=plt.plot(names, values, **{'color': 'firebrick', 'marker': 'o'})
        plt.title(key)
        plt.xlabel('distance')
        plt.ylabel('pseudo energie score')
        plt.savefig(direct+'/'+key+'.png')
        plt.close()



def save_file ( directory,parent,file):
    #try :
        # Directory
        dir = directory
        # Parent Directory path
        par = parent
        # Path
        path = os.path.join(parent, directory)
        x=os.mkdir(path)
        os.chdir(path)
        for key,value in file.items():
            with open(key+'.txt', 'w') as f:
                for k2, v2 in value.items():
                    f.write('%s:%s\n' % (k2, v2))
        f.close()
    #except FileExistsError :
       #return print('the file name : {} already exist in {}, put another file name '.format(directory,parent))




def all_in_one(directory):
    d = dict()
    d_no_round=dict()
    if os.path.isdir(directory):
        print('*** analyse of directory {}'.format(directory))
        for path in os.listdir(directory):
            full_path = os.path.join(directory, path)
            f = open(full_path, 'r')
            text = f.read().split('\n')  # List containing  header
            header = text[0][10:13]
            if header == 'RNA':
                print('*** processing: {} '.format(path))
                A = readFasta(full_path)
                for i in range(0, len(A) - 1):
                    for y in range(1, len(A)):
                        if (float(A[y][5]) - float(A[i][5])) > 3 and A[i][4] == A[y][4]:
                            if distance_no_round(A,i, y) <= 20:
                                add_values_in_dict( A, d, (A[i][3]) + ((A[y][3])), distance(A ,i, y))
                                add_values_in_dict( A, d_no_round, (A[i][3]) + ((A[y][3])), distance_no_round(A ,i, y))
                            else:
                                pass
            else:
                return print( '*** the pdb file : {} , is not RNA ! '.format(full_path) )
    elif os.path.isfile(directory):
        # print('analyse of file {}'.format(directory))
        f = open(directory, 'r')
        text = f.read().split('\n')  # Liste contenant le header et les lignes de la séquence
        header = text[0][10:13]
        if header == 'RNA':
            print(' *** processing: {} '.format(os.path.basename(directory)))
            A = readFasta(directory)
            for i in range(0, len(A) - 1):
                for y in range(1, len(A)):
                    if (float(A[y][5]) - float(A[i][5])) > 3 and A[i][4] == A[y][4]:
                        if distance_no_round(A, i, y) <= 20:
                            add_values_in_dict(A, d, (A[i][3]) + ((A[y][3])), distance(A, i, y))
                            add_values_in_dict(A, d_no_round, (A[i][3]) + ((A[y][3])), distance_no_round(A, i, y))
                        else:
                            pass
        else:
            return print('*** the pdb file : {} , is not RNA ! '.format(directory))
    else :
        return print('no')

    q = {k: {i: v.count(i) for i in range(0, 21)} for k, v in d.items()}
    for k, v in q.items():
        q[k]['sum'] = sum(v.values())

    observed = nested_dic_observed_frequency(q)
    reference = nested_dic_reference_frequency(q)
    pseudo_NRJ = score_pseudo_energy(observed, reference)

    return pseudo_NRJ , d_no_round


def readFasta(seqFile):
    if os.path.isfile(seqFile):
        storage=[]
        with open(seqFile) as ifile:
            for ligne in ifile:
                if ligne[0:4] == 'ATOM' :
                    df_list=[ligne[:6], ligne[6:11], ligne[12:16], ligne[17:20], ligne[21], ligne[22:26], ligne[30:38],
                             ligne[38:46],ligne[46:54]]
                    delete_empty = [x.strip(' ') for x in df_list ]
                    if delete_empty [2] == "C3'":
                        storage.append(delete_empty)
            return storage
    else:
        print(" *** Erreur : le chemin donne '{}' n'est pas valide.".format(
            seqFile))  # if  file is not found  error message




def add_values_in_dict( file_name, sample_dict, key, list_of_values):
    ''' Append multiple values to a key in
        the given dictionary '''
    A = file_name
    # sample_dict=d
    #key = (A[i][3]) + ((A[y][3]))
    # list_of_values=distance(i, y)
    if key not in sample_dict:
        if key[::-1] not in sample_dict:
            sample_dict[key] = list()
        else:
            key = key[::-1]
    sample_dict[key].append(list_of_values)
    return sample_dict


def distance(file_name,i, y):
    A = file_name
    x_dist = (float(A[i][6]) - float(A[y][6])) ** 2
    y_dist = (float(A[i][7]) - float(A[y][7])) ** 2
    z_dist = (float(A[i][8]) - float(A[y][8])) ** 2
    return round(math.sqrt(x_dist + y_dist + z_dist))


def distance_no_round(file_name,i, y):
    A = file_name
    x_dist = (float(A[i][6]) - float(A[y][6])) ** 2
    y_dist = (float(A[i][7]) - float(A[y][7])) ** 2
    z_dist = (float(A[i][8]) - float(A[y][8])) ** 2
    return math.sqrt(x_dist + y_dist + z_dist)


def nested_dic_observed_frequency(nested_dic):
    observed = copy.deepcopy(
        nested_dic)  # need to  create deepcopy otherwise global change occur to the dic argument
    for k, v in observed.items():
        for ki, vi in v.items():
            if ki in range(0, 21):
                v[ki] /= observed[k]['sum']
    return observed


def nested_dic_reference_frequency(dic):
    previous_sums = defaultdict(int)
    dict3 = copy.deepcopy(dic)
    for key, value in dict3.items():
        for inner_key in value:
            value[inner_key] += previous_sums[inner_key]
            previous_sums[inner_key] = value[inner_key]
    reference = {key: value / previous_sums['sum'] for key, value in previous_sums.items()}
    return reference


def score_pseudo_energy(observed, reference):
    dict3 = copy.deepcopy(observed)
    for key, value in dict3.items():
        for innerkey in value:
            if innerkey in reference:
                try:
                    value[innerkey] = - math.log(value[innerkey]/reference[innerkey])
                except ZeroDivisionError:
                    value[innerkey] = 10
                except ValueError:
                    value[innerkey] = 10
        dict3[key].pop('sum', None)
    return dict3


def new_pseudo_linear_interpolation(d_no_round,pseudo_NRJ):
    #calculate score pseudo NRJ #
    new_pseudo=dict()
    score_new_pseudo=0
    # check if d key is not  the same as pseudo_NRJ key
    for key, value in d_no_round.items():
        if key not in pseudo_NRJ:
            new_key = key[::-1]
            d_no_round[new_key] = d_no_round.pop(key)
    for new_key, value in d_no_round.items():
        new_pseudo[new_key] = list()
        for i in range ( len (value)):
            x = value[i]
            x1=math.floor(value[i])  # tail and floor
            x2=math.ceil(value[i])
            y1=pseudo_NRJ[new_key][x1]
            y2=pseudo_NRJ[new_key][x2]
            interpolation_lineaire=y1 + (x - x1) * ((y2 - y1) / (x2 - x1))
            new_pseudo[new_key].append(interpolation_lineaire)
    for key, value in new_pseudo.items():
        score_new_pseudo += sum(new_pseudo[key])
    return  score_new_pseudo

