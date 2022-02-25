#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 10:15:51 2021

@author: gabrielsoto
"""

import os, re, pandas, openpyxl
import numpy as np
from util.FileMethods import FileMethods
import pyomo as pyo

class FileParser(object):
    """
    The FileParser class is a util class used generally to manipulate text data
    in files. Some main features is the ability to search existing Python class
    files for particular strings. Currently using regex to search for complex
    strings, don't yell at me please.
    """
    
    # Regex rules for parsing through dispatch files (e.g. GeneralDispatch, NuclearDispatch)
    param_latex_regex  = r'\s#([\w+\d+\{\}\^\\\,-]+):\s'           # Matches parameter names in LaTeX form
    param_txt_regex    = r'self.model.(\w+)\s?='                   # Matches parameter names in txt form
    detail_latex_regex = r':\s([\w+\s+\d+\{\}\[\]\^\\\,\$/-]+)'    # Matches explanation of parameters in LaTeX form
    section_regex      = r'###\s([\w+\s+]+)\s###'                  # Matches section titles in txt form
    subsection_regex   = r'#-------\s([\w+\s+]+)\s---------'       # Matches subsection titles in txt form

    
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
                #close file
                fo.close()
                return line_number,line_contents
        # close file
        fo.close()
        return None,None
    

    def get_lines_between_strings(fpath, str_start, str_end ):
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
    
    
    def get_dispatch_string_properties(dispatch_name,return_case=0):
        """ Method to return relevant dispatch string properties
        
        This method parses through a given file found at the input filepath
        and returns a dictionary with string values and line numbers for different
        parameter properties in the dispatcher. For now, it can sort through the 
        Params and Variables of the Dispatch files. It returns latex+text names for all
        entries, latex strings containing parameter details, and section+subsection
        titles.
        
        NOTE: this method is tied to the very specific aesthetic format of the Dispatch
        files, should they change this method will no longer work.
        
        Inputs:
            dispatch_name (str)  : name of Dispatch class to document
            return_case (int)    : 0=return Param properties, 1=return Var properties
        Outputs:
            P (dict)   : dictionary with entries for line number and str contents
            
        """
        
        #TODO: check that dispatch name exists
        filepath  = os.path.join( FileMethods.samsim_dir, 'dispatch', dispatch_name)
        filepath += '.py'
        
        # methods to be used as book-ends
        method_names = ['def generate_params',
                        'def generate_variables',
                        'def add_objective']
        
        # define book-end strings that define where all parameters are located
        param_method_start = method_names[0] if return_case == 0 else method_names[1]
        param_method_end   = method_names[1] if return_case == 0 else method_names[2]
    
        # get list of relevant lines to parse through
        lines_list = FileParser.get_lines_between_strings(filepath, param_method_start, param_method_end)
        
        # define regex (regular expression) rules for matching strings
        params_latex_rule   = re.compile( FileParser.param_latex_regex )
        params_txt_rule     = re.compile( FileParser.param_txt_regex )
        detail_latex_rule   = re.compile( FileParser.detail_latex_regex )
        section_rule        = re.compile( FileParser.section_regex )
        subsection_rule     = re.compile( FileParser.subsection_regex )
        
        # define empty lists for all properties we want
        params_latex   = []
        params_txt     = []
        sections = []
        ssection = []
        detail   = []
        
        # define empty lists for line nu where all properties are located
        params_latex_lineNo  = []
        params_txt_lineNo    = []
        section_lineNo  = []
        ssection_lineNo = []
        detail_lineNo   = []
        
        # define empty dict to store and return everything in
        P = {}
        
        # begin search
        count = 0
        for line in lines_list:
            
            # find strings per the regex rules
            param_latex_str   = params_latex_rule.findall(line)
            param_txt_str     = params_txt_rule.findall(line)
            section_str       = section_rule.findall(line)
            subsection_str    = subsection_rule.findall(line)
            detail_str        = detail_latex_rule.findall(line)
            
            # append strs and line numbers to lists if found
            if len(param_latex_str) != 0:
                params_latex.append(param_latex_str[0])
                params_latex_lineNo.append(count)

            if len(param_txt_str) != 0:
                params_txt.append(param_txt_str[0])
                params_txt_lineNo.append(count)
            
            if len(section_str) != 0:        
                sections.append(section_str[0])
                section_lineNo.append(count)
                
            if len(subsection_str) != 0:        
                ssection.append(subsection_str[0])
                ssection_lineNo.append(count)
            
            if len(detail_str) != 0:
                # getting rid of extra \n's at the end of the lines
                if detail_str[0].endswith('\n'):
                    detail.append(  detail_str[0].split('\n')[0]  )
                else:
                    detail.append(detail_str[0])
                detail_lineNo.append(count)
                
            count += 1
        
        P['params_latex'] = params_latex
        P['params_txt']   = params_txt
        P['sections']     = sections
        P['ssection']     = ssection
        P['detail']       = detail

        P['params_latex_lineNo'] = params_latex_lineNo
        P['params_txt_lineNo']   = params_txt_lineNo
        P['section_lineNo']      = section_lineNo
        P['ssection_lineNo']     = ssection_lineNo
        P['detail_lineNo']       = detail_lineNo
        
        P['N_lines'] = count
        
        return P
    
    #TODO: add file to write LaTeX scripts with all parameter/variable names?
    
    
    def get_single_pyomo_model_data(model, params_list, ind):
        """ Method to return single parameter data from a Pyomo model
        
        This method parses through a given Pyomo model and returns a value for
        a parameter from params_list specified by the index ind. 
        
        Inputs:
            model (PyomoModel) : Pyomo model after execution
            params_list (str)  : list of strings representing parameter names
            ind (int)          : index of params_list to return 
        Outputs:
            modelData (int or float)   : value of requested Pyomo parameter
            
        """
        OrderedSimpleSet = pyo.core.base.set.OrderedSimpleSet
        SimpleParam  = pyo.core.base.param.SimpleParam
        ScalarParam  = pyo.core.base.param.ScalarParam
        IndexedParam = pyo.core.base.param.IndexedParam
        
        # check that the model has the Parameter defined
        if hasattr(model, params_list[ind]):
            # get the model Parameter
            modelParam = getattr(model, params_list[ind])
            
            # depending on how each parameter was initialized, there are 
            #    different ways to retrieve them
            if type(modelParam) is OrderedSimpleSet:
                modelData = modelParam.data()
            elif type(modelParam) is SimpleParam or type(modelParam) is ScalarParam:
                modelData = modelParam.value
            elif type(modelParam) is IndexedParam:
                ParamDict = modelParam.extract_values()
                modelData = [ParamDict[x] for x in ParamDict.keys()]
            else:
                # Parameter couldn't be found, return None for now
                #TODO: raise some sort of error here
                modelData = []
        
            return modelData  
    
        else:
            return []
    
    
    def get_list_of_pyomo_data(model, params_list):
        """ Method to return all parameter data from a Pyomo model
        
        This method parses through a given Pyomo model and returns a value for
        each parameter given in the params_list. Typically should parse through
        the DispatchModel first using self.get_dispatch_string_properties() to
        get params_list.
        
        Inputs:
            model (PyomoModel) : Pyomo model after execution
            params_list (str)  : list of strings representing parameter names
        Outputs:
            data_list (list of int or float)   : list of data values for each Pyomo parameter
            
        """
        N = len(params_list)
        data_list = []
        
        for n in range(N):
            single_data_list = FileParser.get_single_pyomo_model_data(model, params_list, n)
            data_list.append(single_data_list)
        
        return data_list


    def write_pyomo_params_to_excel(model, dispatch_name, xls_file_name):
        """ Method to write all Pyomo parameters to an Excel sheet
        
        This method parses through a given Pyomo model and returns a value for
        each parameter given in the params_list. Then it saves the data to an
        Excel sheet called "Raw Data" in an xlsx file specified by the user input.
        Output file will be written in the simulations/outputs folder. Excel sheet 
        will have three columns: Parameter name, its value, and a txt description.
        
        TODO: perhaps this method belongs in FileMethods?
        
        Inputs:
            model (PyomoModel)  : Pyomo model after execution
            dispatch_name (str) : name of Dispatch model in the simulations/dispatch folder
            xls_file_name (str) : name of desired output file (just filename, no directories)

        """
        
        # parse through Dispatch class and retrieve string properties of Parameters
        P = FileParser.get_dispatch_string_properties(dispatch_name,0)
        
        # get the text name for each parameter and details
        params_txt = P['params_txt']
        param_details = P['detail']
        
        # lambda function to facilitate getting the pyomo parameter data
        get_data = lambda x: FileParser.get_single_pyomo_model_data(model, params_txt, x)
            
        # create dataframe with parameter text name, data value, and text detail of that parameter
        dataframe_list = [ [params_txt[n], get_data(n), param_details[n+3]] for n in range(len(params_txt))  ]
    
        # create pandas DataFrame with dataframe list and column names
        params_dataframe = pandas.DataFrame(dataframe_list,columns=['Parameter', 'Value', 'Description'])
        
        # define the output file destination in the simulation/outputs folder
        xls_path = os.path.join(FileMethods.samsim_dir, 'outputs', xls_file_name)
        sheetname = 'Raw Data' # define worksheet name for dumping data
        
        # if the workbook doesn't already exist with given filename, create it
        if not os.path.exists(xls_path):
            wb = openpyxl.Workbook()
            wb.save(filename=xls_path)
            wb.close()
        
        # open workbook with openpyxl
        with pandas.ExcelWriter(xls_path, engine='openpyxl', mode='a') as writer: 
            workBook = writer.book
            # trying here to overwrite sheet if it exists, otherwise create and write data to it
            try:
                workBook.remove(workBook[sheetname])
            except:
                print("Worksheet does not exist")
            finally:
                params_dataframe.to_excel(writer, sheet_name=sheetname,index=False)
            # save the excel file and we're done!
            writer.save()


    def write_pyomo_params_to_text(model, dispatch_name, txt_file_name, param_or_var=0):
        """ Method to write all Pyomo parameters to an text file with LaTeX
        
        This method parses through a given Pyomo model and returns parameter names 
        and information. In a text file specified by the user, it writes LaTeX formatted 
        sections and subsections of the parameter listings. 
        
        Inputs:
            model (PyomoModel)  : Pyomo model after execution
            dispatch_name (str) : name of Dispatch model in the simulations/dispatch folder
            txt_file_name (str) : name of desired output file (just filename, no directories)
            param_or_var (int)  : select either Parameters or Variables to write 

        """
        # get dictionary of parsed parameter string properties
        P = FileParser.get_dispatch_string_properties(dispatch_name,param_or_var)
        
        # assign line contents for each type of line
        params_latex = P['params_latex'] 
        sections = P['sections']
        subsection = P['ssection'] 
        detail = P['detail']
        
        # assign line numbers where each type of line shows up
        params_latex_lineNo = P['params_latex_lineNo'] 
        section_lineNo    = P['section_lineNo']  
        subsection_lineNo = P['ssection_lineNo']

        # total amount of lines
        n_lines = P['N_lines']
        
        # defining LaTeX commands for beginning and ending sections/alignings
        begin_section    = '\subsection{'
        begin_subsection = '\subsubsection{'
        begin_align      = '\\begin{flalign}'
        end_align        = '\end{flalign}'
        
        # start counting at 0
        doc_lines = []
        count_sec  = 0
        count_subsec = 0
        count_params = 0
        
        for n in range(n_lines):
            # writing section titles
            if n in section_lineNo:
                tmp_section = begin_section + sections[count_sec] + '}'
                doc_lines.append(tmp_section)
                count_sec += 1
                
                # if params in next line, write a \begin{align} statement
                if n+1 in params_latex_lineNo:
                    doc_lines.append(begin_align)
            
            # writing subsection titles
            if n in subsection_lineNo:
                tmp_ssection = begin_subsection + subsection[count_subsec] + '}'
                doc_lines.append(tmp_ssection)
                count_subsec += 1
                
                # write a \begin{align} statement for params next line
                doc_lines.append(begin_align)
            
            # writing LaTeX representation of params with text details
            if n in params_latex_lineNo:
                tmp_params = '&' + params_latex[count_params] + ' &&: \\text{' + detail[count_params] + '}'
                
                # if there's another parameter next, end the line
                if n+1 in params_latex_lineNo:
                    tmp_params += ' \\\\ '
                    doc_lines.append(tmp_params)
                # if no parameters up next, end the align brackets
                else:
                    doc_lines.append(tmp_params)
                    doc_lines.append(end_align + '\n')
           
                count_params += 1
        
        # print and/or save text file
        txt_path = os.path.join(FileMethods.samsim_dir, 'outputs', txt_file_name)
        text_file = open(txt_path,'w')
        text_file.write(' \n'.join(doc_lines))
        text_file.close()


    def write_pyomo_params_to_table(model, dispatch_name, txt_file_name, param_or_var=0):
        """ Method to write all Pyomo parameters to an text file with LaTeX
        
        This method parses through a given Pyomo model and returns parameter names 
        and information. In a text file specified by the user, it writes LaTeX formatted 
        sections and subsections of the parameter listings. 
        
        Inputs:
            model (PyomoModel)  : Pyomo model after execution
            dispatch_name (str) : name of Dispatch model in the simulations/dispatch folder
            txt_file_name (str) : name of desired output file (just filename, no directories)
            param_or_var (int)  : select either Parameters or Variables to write 

        """

        # get dictionary of parsed parameter string properties
        P = FileParser.get_dispatch_string_properties(dispatch_name, param_or_var)
        
        # assign line contents for each type of line
        params_latex = P['params_latex'] 
        sections = P['sections']
        subsection = P['ssection'] 
        detail = P['detail']
        
        # assign line numbers where each type of line shows up
        params_latex_lineNo = P['params_latex_lineNo'] 
        section_lineNo    = P['section_lineNo']  
        subsection_lineNo = P['ssection_lineNo']
        params_detail_lineNo = P['detail_lineNo'] 

        # total amount of lines
        n_lines = P['N_lines']
        
        # start counting at 0
        doc_lines = []
        count_sec  = 0
        count_subsec = 0
        count_params = 0
        
        for n in range(n_lines):
            # writing section titles
            if n in section_lineNo:
                tmp_section = '& &\\ \'' + sections[count_sec] + '\''
                doc_lines.append(tmp_section)
                count_sec += 1
                
                # if params in next line, write a \begin{align} statement
                if n+1 in params_latex_lineNo:
                    doc_lines.append(' ')
                

            # writing subsection titles
            if n in subsection_lineNo:
                tmp_ssection = '& &\\ \'' + subsection[count_subsec]  + '\''
                doc_lines.append(tmp_ssection)
                count_subsec += 1
                
                doc_lines.append(' ')
                

            # writing LaTeX representation of params with text details
            if n in params_latex_lineNo:
                detail_line_match = np.where( np.asarray(params_detail_lineNo) == params_latex_lineNo[count_params])[0][0]
                tmp_detail = detail[ int(detail_line_match)  ]
                
                # if there's another parameter next, end the line
                if n+1 in params_latex_lineNo:
                    if '[' in tmp_detail:
                        detail_txt, unit_txt = tmp_detail.split('[')
                        tmp_params = '${0}$ & ${1}$ & {2} \\\\'.format(
                            params_latex[count_params], unit_txt[:-1], detail_txt )
                    else:
                        tmp_params = '&' + params_latex[count_params] + ' &&: \\text{' + detail[count_params] + '}'
                        
                    doc_lines.append(tmp_params)
                # if no parameters up next, end the align brackets
                else:
                    doc_lines.append(tmp_params)
                    doc_lines.append(' ')

                count_params += 1
        
        # print and/or save text file
        txt_path = os.path.join(FileMethods.samsim_dir, 'outputs', txt_file_name)
        text_file = open(txt_path,'w')
        text_file.write(' \n'.join(doc_lines))
        text_file.close()