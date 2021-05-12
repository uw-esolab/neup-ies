#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 10:15:51 2021

@author: gabrielsoto
"""

import os, re
from util.FileMethods import FileMethods

class FileParser(object):
    
    def find_first_string_instance(fpath, search_str):
        """ Method to parse through file and locate string
        
        This method parses through a given file found at the input filepath
        and searches for the first instance of a given string. It then returns
        the line number and the full line string where it was found. If not found
        it returns None for both outputs. 
        
        Inputs:
            fpath (str)      : full path to relevant file
            search_str (str) : string to search
        Outputs:
            line_number (str)   : line number where search_str was found (None if not found)
            line_contents (str) : full line string where search_str was found (None if not found)
            
        """
        
        fo = open(fpath)
        # line-by-line search
        for line_number,line_contents in enumerate(fo):
            found = search_str in line_contents
            # found the search string in the file!
            if found:
                print(fpath, "[", line_number, "] ")
                print(line_contents)
                #close file
                fo.close()
                return line_number,line_contents
        # close file
        fo.close()
        return None,None