'''
Created on 2010-4-9

@author: Madlee
'''
import os.path as path

MAX_EXT_NAME_LENGTH = 10

def extname(file):
    name = path.basename(file)
    index = name.rfind('.') 
    if index != -1:
        result = name[index+1:]
        if len(result) < MAX_EXT_NAME_LENGTH:
            return result.lower()
    return ''

def remove_extname(file):
    folder, name = path.split(file)
    index = name.rfind('.') 
    if index != -1:
        name = name[0:index]
    return path.join(folder, name)

def create_array(file, typename):
    result = []
    fields = [i.strip() for i in file.readline().split(',')]
    line = file.readline()
    while line is not None and len(line) > 0:
        data = eval(line)
        assert len(data) == len(fields)
        record = typename()
        for i in range(len(data)):
            setattr(record, fields[i], data[i])
        result.append(record)
        line = file.readline()
    return result

def dump_string(filename, content):
    with open(filename, 'wt') as file:
        file.write(content)

def load_string(filename):
    with open(filename) as file:
        return file.read()
    
def load_value(filename):
    return eval(load_string(filename))
