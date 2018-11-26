import re
import numpy as np
import os
class PdbReader():

    def __init__(self, pdb_file, model_num=None, x_start_col=None, element_col=None):
        self.pdb_file = pdb_file
        self.model_num = model_num
        self.x_start_col = x_start_col
        self.element_col = element_col
        self.file_name, self.file_ext = os.path.splitext(str(self.pdb_file))

    def get_ATOM(self):
        f = open(self.pdb_file, "r")
        each_model_list = []
        all_model_list = []
        if self.model_num == None:
            for line in f:
                if line.startswith("ATOM"):
                    each_model_list.append(str(line.strip()))
            f.close()
            MODEL = each_model_list
        else:
            for line in f:
                if line.startswith("MODEL"):
                    each_model_list = []

                elif line.startswith("ATOM"):
                    each_model_list.append(line)

                elif line.startswith("ENDMDL"):
                    all_model_list.append(each_model_list)
            MODEL = all_model_list[int(self.model_num) - 1]
            f.close()
        self.MODEL = MODEL
        return(MODEL, all_model_list)

    def move_center(self, MODEL):
        q_array = np.array([])
        weight_q_array = np.array([], dtype = 'float64')
        moment_array = np.array([], dtype = 'float64')
        Total_weight = 0
        for line in MODEL:
            if self.x_start_col == None:
                x = line.split()[6]
                y = line.split()[7]
                z = line.split()[8]
            else:
                x = line.split()[self.x_start_col]
                y = line.split()[self.x_start_col+1]
                z = line.split()[self.x_start_col+2]
            if self.element_col == None:
                element= line.split()[11]
            else:
                element= line.split()[self.element_col]

            q = np.array([x, y, z], dtype = 'float64')
            q_array = np.append(q_array, q)

            if element == 'H':
                weight = 1
            elif element == 'C':
                weight = 12
            elif element == 'N':
                weight = 14
            elif element == 'O':
                weight = 16
            elif element == "S":
                weight = 32
            else:
                print("[[[WARNING]]] Please write the weight of Element:" + element)
            moment = q * weight
            moment_array = np.append(moment_array, moment)
            Total_weight += weight
            weight_q = q * (weight**(1/2))
            weight_q_array = np.append(weight_q_array, weight_q)
        q_matrix = q_array.reshape(int(q_array.shape[0]/3), 3)
        weight_q_matrix=  weight_q_array.reshape(int(weight_q_array.shape[0]/3), 3)
        moment_matrix = moment_array.reshape(int(moment_array.shape[0]/3), 3)
        center = np.array([sum(moment_matrix[:, 0]), sum(moment_matrix[:, 1]), sum(moment_matrix[:, 2])]) / Total_weight
        q_matrix_from_center = q_matrix - center
        self.q_matrix_from_center = q_matrix_from_center
        return(q_matrix_from_center)

    def file_changer(self):
        count = -1
        f = open(self.pdb_file, "r")
        new_file = open(self.file_name + '_from_center' + self.file_ext, 'w')
        if self.model_num == None:
            for line in f:
                if line.startswith("ATOM"):
                    count += 1
                    if self.x_start_col == None:
                        x = line.split()[6]
                        y = line.split()[7]
                        z = line.split()[8]
                    else:
                        x = line.split()[self.x_start_col]
                        y = line.split()[self.x_start_col+1]
                        z = line.split()[self.x_start_col+2]
                    n = 0
                    new_line = line
                    for be in str(x), str(y), str(z):
                        len_be = len(be)
                        af = str(round(self.q_matrix_from_center[count][n], 3))
                        len_af = len(af)
                        if len_be == len_af:
                            new_line = new_line.replace(be, af)
                        elif len_be > len_af:
                            new_line = new_line.replace(be, ' '*(len_be-len_af)+af)
                        elif len_be < len_af:
                            new_line = new_line.replace(' '*(len_af-len_be) + be, af)
                        n += 1
                else:
                    new_line = line
                print(new_line.strip())
                new_file.write(new_line)
        f.close()
        new_file.close()

# Rotate PdbClass1 matrix to Pdbclass2 and output the rotated matrix to an file.
def Rotate(PdbReaderClass1, PdbReaderClass2):
    D = []
    for i in [0, 1, 2]:
        for j in [0, 1,2 ]:
            D.append(np.sum(PdbReaderClass1.q_matrix_from_center[:, i] * PdbReaderClass2.q_matrix_from_center[:, j]))
    D = np.array(D)
    D_matrix = D.reshape(3, 3)
    U, S, V = np.linalg.svd(D_matrix, full_matrices=False)
    R = np.dot(U, V)
    Rotated = np.dot(PdbReaderClass1.q_matrix_from_center, R)
    PdbReaderClass1.file_changer(Rotated, "_rotated")