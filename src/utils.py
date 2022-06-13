# -*- coding: utf-8 -*-

'''
utils.py -- Definition of useful basic functions

Written by Moisès Bernabeu <mosies.bernabeu.sci@gmail.com>
2022
'''

import subprocess as sp
import sys
import os


# Running bash commands
def run_cmd(cmd, ommit=False):
    '''
    Running a bash command from python

    This function gets the command to execute and tries to execute it, if ommit
    is tre, it communicates you to ommit or not, while if it is false, it will
    ommit in case of it is not running.

    Args:
        cmd (str): string of the command
        ommit (boolean): indicate whether ommit or not.

    Returns:
        int: 0
    '''

    if ommit:
        try:
            process = sp.Popen(cmd, shell=True)
        except:
            pass
        process.communicate("Y\n")
        if process.wait() != 0:
            print("Error occurred, but you chose to ommit it")
    else:
        try:
            process = sp.Popen(cmd, shell=True)
        except OSError:
            sys.exit("Error: Execution cmd failed")
        process.communicate("Y\n")
        if process.wait() != 0:
            sys.exit("ERROR: Execution cmd failed")

    return 0


def run_cmd_return(cmd):
    '''
    Running a bash command from python returning the output

    This function gets the command to execute and tries to execute it. The
    output is stored and returned by the function.

    Args:
        cmd (str): string of the command

    Returns:
        str: output
    '''

    try:
        process = sp.Popen(cmd, shell=True, stdout=sp.PIPE).stdout
    except:
        sys.exit()

    return process.readlines()


# Creating folders
def create_folder(name):
    '''
    This function checks the existance of a folder and if it is not created
    it creates

    Args:
        name (str): string with the folder name

    Returns:
        int: 0
    '''

    if not os.path.exists(name):
        try:
            os.mkdir(name)
        except Exception:
            print('Unable to create the directory: %s' % name)

    return 0


def rmkey(dict, key):
    '''
    This function removes an element of a dictionary based on the key

    Args:
        dict (dict): dictionary
        key (str): string of the dictionary key

    Returns:
        dict: dictionary without the removed key
    '''

    del dict[key]

    return dict


def rmbwbr(string, upper_del, lower_del):
    '''
    Removes elements between upper and lower delimiters

    Gets the string and counts the number of times that the upper string
    appears then if it only appears once, it cuts there and in the lower, and
    retrieves the string containing the elements outside the delimiters. if
    there are more than 1, it returns the element after last upper delimiter

    Args:
        string (str): string of interest containing upper and lower delimiters
        upper_del (str): string with the upper delimiter
        lower_del (str): string with the lower delimiter

    Returns:
        str: string with all the elements between the lower
             and upper delimiters
    '''

    nsets = string.count(upper_del)
    ostr = ''
    if nsets == 1:
        ostr = string.split(upper_del)[0]
        ostr += string.split(lower_del)[1]
    else:
        parts = string.split(upper_del)
        ostrl = list()
        for i, item in enumerate(parts):
            if not i % 2:
                ostrl.append(item.replace(lower_del, ''))
        ostr = ' '.join(ostrl)
        ostr.replace('  ', ' ').replace('..', '.').replace(',,', ',')

    return ostr


def csv_to_dict(file, sep):
    '''
    This function converts a 2 columns data frame in a dictionary

    It iterates over the file lines and stores the first column as the index
    and the second as the value of a dictionary.

    Args:
        file (str): path to the 2 columns table file
        sep (str): the separator of the table file

    Returns:
        dictionary: dictionary which keys are the first column and the values
                    are the second column
    '''

    odict = dict()
    for line in open(file):
        # Checking whether the line is commented
        if line != '' and '#' not in line:
            # Getting the elements
            elements = line.replace('\n', '').split(sep)
            # Writing to the dictionary
            odict[elements[0]] = elements[1]

    return odict


def file_exists(filename):
    '''
    This function checks the file existance

    The function gets the path of a file and checks whether the file exists
    and if it has more than a byte

    Args:
        filename (str): filename path

    Returns:
        boolean: true for the existance and non-emptyness of the file
    '''

    # Check file existence
    if os.path.isfile(filename):
        # Check the file is not empty
        if os.stat(filename).st_size > 1:
            return True
    else:
        return False


if __name__ == '__main__':
    istr = 'asjkdh $ abc cda : miau'
    print(rmbwbr(istr, '$', ':'))
