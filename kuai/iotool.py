# -*- coding:utf-8 -*-
"""
Created on 2010-1-2

@author: madlee
"""

from cStringIO import StringIO
from os.path import abspath
from kuai.strtool import extname

class Reader:
    def read(self, fileobj):
        raise NotImplementedError()
    
    def read_file(self, filename):
        with open(filename) as file:
            result = self.read(file)
            result.filename = abspath(filename)
            return result

    def parse(self, data):
        file = StringIO(data)
        try:
            result = self.read(file)
            return result
        finally:
            file.close()
        
class ReaderManager:
    def __init__(self):
        self.readers = {}
    
    def read(self, file, type):
        reader = self.get_reader(type)
        return reader.read(file)
    
    def read_file(self, filename, type=None):
        if type is None:
            type = extname(filename)
        reader = self.get_reader(type)
        return reader.read_file(filename)
    
    def parse(self, data, type):
        reader = self.get_reader(type)
        return reader.parse(data)

    def get_reader(self, type):
        if type in self.readers:
            return self.readers[type]
        elif '' in self.readers:
            return self.readers['']
        else:
            raise KeyError('Can not find a reader for %s file' % type)


class Writer:
    def write(self, fileobj, var):
        raise NotImplementedError()
    
    def write_file(self, filename, var):
        with open(filename, 'wt') as file:
            self.write(file, var)

    def format(self, var):
        file = StringIO()
        try:
            self.write(file, var)
            result = file.getvalue()
            return result
        finally:
            file.close()


class WriterManager:
    def __init__(self):
        self.writers = {}
    
    def write(self, file, var, type):
        writer = self.get_writer(type)
        return writer.write(file, var)
    
    def write_file(self, filename, var, type=None):
        if type is None:
            type = extname(filename)
        writer = self.get_writer(type)
        return writer.write_file(filename, var)
    
    def format(self, var, type):
        writer = self.get_writer(type)
        return writer.format(var)
    
    def get_writer(self, type):
        if type in self.writers:
            return self.writers[type]
        elif '' in self.writers:
            return self.writers['']
        else:
            raise KeyError('Can not find a writer for %s file' % type)

def check_error(line, file):
    if line.startswith('#ERROR:'):
        error = [line[7:]]
        line = file.readline()
        while line.strip() != '':
            error.join(line)
            line = file.readline()
        raise RuntimeError(''.join(error))
