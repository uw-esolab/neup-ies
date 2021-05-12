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
                print("Line: ", line_contents)
                #close file
                fo.close()
                return line_number,line_contents
        # close file
        fo.close()
        return None,None
    

    def grab_lines_between_strings(fpath, str_start, str_end ):
        """ Method to return all lines between two strings in a file
        
        This method parses through a given file found at the input filepath
        and returns a list of lines in between two book-end strings: str_start
        and str_end. 
        
        Inputs:
            fpath (str)      : full path to relevant file
            str_start (str)  : first string to search
            str_end (str)    : last string to search
        Outputs:
            lines (list of str)   : list of lines between book-end strings
            
        """
        
        fo = open(fpath)
        # find starting and ending lines
        i_start, L_start = FileParser.find_first_string_instance(fpath, str_start)
        i_end,   L_end   = FileParser.find_first_string_instance(fpath, str_end)
        # create empty list for lines
        lines = []
        # if that the start and end strings were actually found, otherwise return empty list
        if i_start is not None and i_end is not None:
            # line-by-line search
            for i,l in enumerate(fo):
                if i>i_start and i<i_end:
                    lines.append(l)
        # close file
        fo.close()
        return lines
    
    