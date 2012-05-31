# -*- coding:utf-8 -*-
"""
Created on 2010-1-9

@author: madlee
"""

try:
    from settings import KUAI_JOB_EXE
except:
    KUAI_JOB_EXE = r"D:\Madlee\Job\Madlee\kuai\bin32w\kuaijob.exe"

from subprocess import Popen, PIPE, STDOUT
from kuai.iotool import check_error

import subprocess
kwargs = {}
if subprocess.mswindows:
    su = subprocess.STARTUPINFO()
    su.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    su.wShowWindow = subprocess.SW_HIDE
    kwargs['startupinfo'] = su

KUAI_JOB_PIPE = Popen("kuaijob.exe", 100*1024, KUAI_JOB_EXE, stdin = PIPE, stdout = PIPE, stderr = STDOUT, **kwargs)
KUAI_JOB_PIPE.stdin.write("job=inchi\n")

def inchi(mol):
    mol = str(mol)
    KUAI_JOB_PIPE.stdin.write(".mol=embed\n")
    KUAI_JOB_PIPE.stdin.write(mol)
    KUAI_JOB_PIPE.stdin.write("end\n")
    KUAI_JOB_PIPE.stdin.flush()
    
    line = KUAI_JOB_PIPE.stdout.readline()
    while line:
        check_error(line, KUAI_JOB_PIPE.stdout)
        if line.startswith('Structure:'):
            break
        line = KUAI_JOB_PIPE.stdout.readline()
    
    inchi_code=KUAI_JOB_PIPE.stdout.readline()
    check_error(inchi_code, KUAI_JOB_PIPE.stdout)
    inchi_code = inchi_code.rstrip()
    
    aux_code = KUAI_JOB_PIPE.stdout.readline()
    check_error(aux_code, KUAI_JOB_PIPE.stdout)
    aux_code = aux_code.rstrip()
    
    inchi_key = KUAI_JOB_PIPE.stdout.readline()
    check_error(inchi_key, KUAI_JOB_PIPE.stdout)
    assert inchi_key.startswith('InChIKey=')
    inchi_key = inchi_key[9:].rstrip()
    
    return inchi_code, aux_code, inchi_key
    
