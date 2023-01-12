from utils import *
import os
import sys
import argparse
import textwrap


parser = argparse.ArgumentParser(
      prog='Main.py',
      formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog=textwrap.dedent('''\
        information:
            .Main.py should be run in the same folder as utils.py. Only basics python library have been used, please check requirements.txt for additional information
            .Main.py allow an interactive utilisation. The script take into account different constraint : 
                - Only pdb file 
                - only RNA file
                - File of testing and training shouldn't be the same 
                - directory / file path should be valid
                
            .Dataset should be contained inside a folder . Make sure that this folder ONLY contain PDB RNA file 
            .For training : When PATH is asking , put the right already existing PATH where file(s) are located (ex: D:\PDB\storage_data ) 
            .When Name is asking , put just the name ( ex : myfolder ) , The name should be unique inside the path
            .For testing  : FILE PATH should be use , and not directory (ex: D:\PDB\storage_data\data1.pdb ) 
            .For testing , the final value will be stored inside training path choosen
            .I strongly advise to put testing data into an other directory as the training data, since it shouldn't be the same
            .To run the script effectively , you can just run 1 , as it will train / plot proposal / test proposal . 
            .If you only want to plot , go to 2 
            .Otherwise 3 for exit  
            .the script allow the user to retest if file / path given is wrong or non pdb/rna for every step
         '''))

args = parser.parse_args()


if __name__=="__main__":
    while True :
        xt=input("""
        1) Training script and testing script 
        2) only Plot training script
        3) Leave

        CHOISIR Function >>>>""")

        # Partie 1 #########################################################################
        if xt =='1':
            scan_file = 0
            list_of_valid_files = []

            try:
                print("--- Training script :")
                directory = str(input('*** put the directory PATH to analyse, it should contain ONLY PDB files :'))
                if os.path.isdir(directory):
                    for path in os.listdir(directory):
                        full_path = os.path.join(directory, path)
                        if full_path.endswith(".pdb"):
                            print('*** {} is valid, and can be analysed : '.format(path))
                            list_of_valid_files.append(path)
                        else:
                            scan_file += 1
                            print('*** ERROR  remove non PDB file  {} '.format(path))
                    if scan_file > 0:
                        pass
                    else:
                        name_directory = str(input('*** put PATH of  EXISTING DIRECTORY where file will be stored : '))
                        nom_folder = str(input('*** put NAME of the FOLDER  where file will be stored ,folder will be created; folder name should be UNIQUE : '))
                        print('test')
                        score = all_in_one(directory)[0]
                        save_file(nom_folder, name_directory, score)
                        print('*** the file have been saved inside{} '.format(name_directory + '/'+nom_folder))

                        yesno = str(input('*** Do you want to plot the files ? Type yes/y : '))
                        if yesno == 'yes' or yesno =='y':
                            again = str(input('*** do you want to save the plot in the same folder as used for saving the file ? Type yes/y : '))
                            if again == 'yes' or again == 'y':
                                plot_script_2(score, name_directory + '/' + nom_folder)
                            else:
                                plot_loop = 0
                                while plot_loop == 0 :
                                    try:
                                        fullpath = str(input('*** Put the PATH where plot will be saved : '))
                                        plot_script_2(score, fullpath)
                                        plot_loop = 1
                                    except FileNotFoundError:
                                        print('directory not found ')
                                        retry = str(input(' do you wanna retry ? Type yes/y : '))
                                        if retry !='yes' and retry !='y':
                                            print('No plotting ')
                                            plot_loop = 1
                        else:
                            print('No plotting')

                        print('--- Testing script : ')
                        test_script= str(input('*** do you want to perform testing script ? Types yes/y : '))
                        if test_script == 'yes' or test_script == 'y':
                            path_test_file = str(input('*** Put file path to analyse , file should be different than the file(s) use for the training script ! : '))
                            user_multiple_test=0
                            while user_multiple_test == 0:
                                if os.path.isfile(path_test_file):
                                    if path_test_file.endswith(".pdb") and os.path.basename(path_test_file) not in list_of_valid_files :
                                        print('*** {} is a valid pdb file  , and can be analysed '.format(os.path.basename(path_test_file)))
                                        dic_test = all_in_one(path_test_file)[1]
                                        pseudo_NRJ = all_in_one(path_test_file)[0]
                                        testing=new_pseudo_linear_interpolation(dic_test, pseudo_NRJ)
                                        print('*** value of testing script ',testing)
                                        print('*** value testing is store inside folder of training {}' . format(name_directory + '/'+nom_folder + '/' + os.path.basename(path_test_file)+'_testing'+'.txt'))
                                        # open text file
                                        text_file = open(name_directory + '/'+nom_folder + '/' + os.path.basename(path_test_file)+'_testing'+'.txt' , "w")
                                        # write string to file
                                        text_file.write(str(testing))
                                        # close file
                                        text_file.close()
                                        user_multiple_test = 1
                                    else:
                                        print('*** {}  might not be a PDB file or the file have already been used in the training script !'.format(path_test_file))
                                        test_again_again =str(input('*** do you want to test again for a valid file ? Type yes/y : '))
                                        if test_again_again =='yes' or test_again_again=='y':
                                            path_test_file=str(input('*** Put file path to analyse , file should be different than the file(s) use for the training script and a PDB ! : '))
                                        else:
                                            print('*** stop')
                                            user_multiple_test = 1
                                else:
                                    print('*** can only accept valid and existing file path not directory')
                                    test_again_again = str(input('*** do you want to test again for a valid file ? Type yes/y : '))
                                    if test_again_again == 'yes' or test_again_again == 'y':
                                        path_test_file = str(input(
                                            '*** Put file path to analyse , file should be different than the file(s) use for the training script and a PDB ! : '))
                                    else:
                                        print('*** stop')
                                        user_multiple_test = 1

                else:
                    print('only accept valid directory as input')
            except FileNotFoundError:
                print('directory name not found ')
            except TypeError:
                pass
            except FileExistsError:
                print('file name {} already exist inside directory {}'.format(nom_folder, name_directory))
            except AttributeError:
                print('The script accept only RNA file ')
            except NameError:
                pass
            # except OSError :
            # print('put a valid name ')


        if xt == '2' :
            scan_file = 0
            try:
                print("ONLY PLOTTING :")
                directory = str(input('*** put the directory path to analyse , it should contain ONLY PDB files :'))
                if os.path.isdir(directory):
                    for path in os.listdir(directory):
                        full_path = os.path.join(directory, path)
                        if full_path.endswith(".pdb"):
                            print('*** {} is valid, and can be analysed '.format(path))
                        else:
                            scan_file += 1
                            print('*** ERROR  remove non PDB file  {} '.format(path))
                    if scan_file > 0:
                        pass
                    else:
                        score = all_in_one(directory)[0]
                        plot_storage = str(input('*** put the path to store the plot'))
                        if os.path.exists(plot_storage) :
                            plot_script_2(score, plot_storage)
                        else:
                            print('*** path no existing ')
                else:
                    print('*** directory name non valid ')
            except FileNotFoundError:
                print('*** directory name not found ! ')
            except FileExistsError:
                print('*** file name {} already exist inside directory {}'.format(nom_folder, name_directory))
            # except OSError :
            # print('put a valid name ')

        if xt== '3':
            break







