#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import wx
import numpy as np
import pandas as pd
import os
from TAME_GUI import mainWindow, newProjectDialog, fixedLevelsDialog
import os.path
import bisect
import configparser
import subprocess
from ObjectListView import ColumnDefn, OLVEvent
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
import glob

import warnings  # only here to stop deprecation warning of objectlistview from clogging up terminal
warnings.filterwarnings("ignore", category=DeprecationWarning)


class MyFrame(mainWindow):
    def __init__(self, *args, **kwds):
        mainWindow.__init__(self, *args, **kwds)
        
        self.set_constants()
        
        self.load_main_config()
        self.load_project() 
        
        self.configure_listviews()
        self.configure_layout()
        
        
    def set_constants(self):
        self.cm_1 = 'cm\u207B\u00B9'
        self.blank_strans_lev = {'label': '', 'j':0.0 , 'energy':0.0 , 'parity':0}
   
        
    def configure_layout(self):
                    
        self.window_2.SetSashPosition(780)
        self.window_1.SetSashPosition(380)

        self.toolbar = NavigationToolbar(self.matplotlib_canvas)
        self.toolbar.Realize()
        self.matplotlib_sizer.Add(self.toolbar, 0, wx.ALIGN_CENTRE, border=5)
        self.toolbar.update()
        
        self.lopt_level_panel.Hide()
        self.lopt_line_panel.Hide()
        self.sizer_8.Layout()
        
        
    def configure_listviews(self):
        self.strans_lev_ojlv.SetColumns([
        ColumnDefn("Level", "left", 100, 'label', isEditable=True),
        ColumnDefn("J", "left", 50, 'j', stringConverter="%.1f", isEditable=True),
        ColumnDefn(f"Energy ({self.cm_1})", "left", 120, 'energy', stringConverter="%.4f", isEditable=True),           
        ColumnDefn("Parity", "left", 50, 'parity', stringConverter="%d", isSpaceFilling=True, isEditable=True)])
    
        self.strans_lev_ojlv.cellEditMode = self.strans_lev_ojlv.CELLEDIT_DOUBLECLICK
        self.strans_lev_ojlv.Bind(OLVEvent.EVT_CELL_EDIT_FINISHED, self.on_strans_lev_edit)
        
        self.strans_lines_ojlv.SetColumns([
            ColumnDefn(f'Wavenumber ({self.cm_1})', 'left', 150, 'wavenumber', stringConverter="%.4f"),
            ColumnDefn('SNR', 'left', 40, 'peak', stringConverter="%d"),
            ColumnDefn('FHWM (mK)', 'left', 90, 'width', stringConverter="%d"),
            ColumnDefn('log(Eq. Width)', 'left', 110, 'eq width', stringConverter="%.2f"),
            ColumnDefn('Fit', 'left', 30, 'tags'),
            ColumnDefn(f'Unc. ({self.cm_1})', 'left', 90, 'unc', stringConverter="%.4f"),
            ColumnDefn('Main Element Transitions', 'left', 500, 'main_desig'),
            ColumnDefn('Other Element Transitions', 'left', 500, 'other_desig', isSpaceFilling=True)])
        
        self.strans_lines_ojlv.SetEmptyListMsg("Run STRANS first")
        
        self.group_column = ColumnDefn(f'Level ({self.cm_1})', 'left', 100, 'main_level', stringConverter="%.4f", groupKeyConverter=self.loptGroupKeyConverter)       
        self.lopt_lev_ojlv.SetColumns([
            self.group_column,
            ColumnDefn('', 'left', 20, 'star', stringConverter=self.star_to_asterix),
            ColumnDefn('Fit', 'left', 30, 'tags', stringConverter=self.fit_converter),
            ColumnDefn('Intensity', 'left', 70, 'log_ew', stringConverter=self.log_ew_converter),
            ColumnDefn('SNR', 'left', 50, 'peak', stringConverter=self.snr_converter),
            ColumnDefn(f'Wn ({self.cm_1})', 'left', 100, 'wavenumber', stringConverter=self.wn_converter),
            ColumnDefn(f'Unc. ({self.cm_1})', 'left', 90, 'uncW_o', stringConverter=self.wn_converter),
            ColumnDefn(f'Obs-Calc ({self.cm_1})', 'left', 120, 'dWO-C', stringConverter=self.neg_num_str),
            ColumnDefn('Level', 'left', 100, 'other_level'),
            ColumnDefn('Tags', 'left', 40, 'F', stringConverter=self.correct_tags ,isSpaceFilling=True)])
        
        self.lopt_lev_ojlv.SetEmptyListMsg("Run LOPT first")
        self.lopt_lev_ojlv.SetShowItemCounts(False)
        self.lopt_lev_ojlv.SetAlwaysGroupByColumn(1)
        self.lopt_line_listctrl.EnableCheckBoxes(True)


    def loptGroupKeyConverter(self, energy):
        selected_line_index = self.lopt_levs.loc[self.lopt_levs['Energy'] == energy].index.values[0]
        return self.lopt_levs.at[selected_line_index, 'Designation']
    
    def log_ew_converter(self, log_ew):
        if str(log_ew) == 'nan':
            return '-'
        return f'{log_ew:.2f}'
    
    def wn_converter(self, wn):
        if str(wn) == 'nan':
            return '-'
        return f'{wn:.4f}'
    
    def snr_converter(self, snr):
        if str(snr) == 'nan':
            return '-'
        return f'{snr:.0f}'
    
    def fit_converter(self, fit):
        if str(fit) == 'nan':
            return '-'
        return fit
    
    def star_to_asterix(self, star):
        if star:
            return '*'
        return ''
        
    def correct_tags(self, tag):
        if str(tag) == 'nan':
            return ''
        return tag.upper()
    
    def neg_num_str(self, diff):
        if diff < 0.0:
            return f'{diff:.4f}'
        else:
            return f' {diff:.4f}'
        
    def load_plot_df(self):
        """Loads or creates the matplotlib data from Xgremlin ascii linelist files - maybe better as part of the wizard?"""
        file = f'{self.project_title}_plot.pkl'
        path='/home/christian/Desktop/TASS'

        if not os.path.isfile(file):  # if no existing DataFrame is present
            self.plot_df = pd.concat([pd.read_csv(f, skiprows=4, delim_whitespace=True, names=['wavenumber', f'{f.split("/")[-1].split(".")[0]}']) for f in glob.glob(path + "/*.asc")], ignore_index=True)
            self.plot_df.to_pickle(file) 
        else:
            self.plot_df = pd.read_pickle(file) 
        
        # self.plot_df.sort_values(by=['wavenumber'], ascending=True, inplace=True)
        # self.plot_df.set_index('wavenumber', inplace=True, drop=True)

         
    def display_strans_levs(self):
        """Writes values from a list to the strans_lev_ojlv ObjectListView"""
        self.strans_lev_ojlv.SetObjects(self.strans_levs)

            
    def display_strans_lines(self):
        """Writes lines with designations from self.df to the strans_lines_ojlv ObjectListView"""   
        strans_lines = list(self.df.loc[self.df.main_desig.str.len() > 0 ].transpose().to_dict().values())
        
        for line in strans_lines:
            main_level_string = ''
            other_level_string = ''
            for lev_pair in line['main_desig']:
                if not main_level_string:
                    sep = ''
                else:
                    sep = ';  \t'
                main_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['upper_level'] + ' - ' + lev_pair['lower_level']
                
            for lev_pair in line['other_desig']:
                if not other_level_string:
                    sep = ''
                else:
                    sep = ',     '
                other_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['upper_level'] + ' - ' + lev_pair['lower_level']
                
            line['eq width'] = float(np.log(line['eq width']))
            line['main_desig'] = main_level_string
            line['other_desig'] = other_level_string
            
        self.strans_lines_ojlv.SetObjects(strans_lines)
                    
    
    def create_df(self, lines_file):
        """Creates a new pandas DataFrame from a list of lines in 'lines_file' and saves to a pickle file."""
        self.df = pd.read_csv(lines_file, float_precision='high')  # create new dataframe from the input lines file 
        self.df['main_desig'] = np.empty((len(self.df), 0)).tolist()  # append column of empty lists.
        self.df['other_desig'] = np.empty((len(self.df), 0)).tolist()  # append column of empty lists.
        self.df['user_desig'] = ''
        self.df['line_tags'] = [{'ringing': False, 'noise': False, 'blend': False, 'user_unc': False, 'multiple_lines':False} for x in range(self.df.shape[0])]
        self.df['comments'] = ''
        self.save_df()
        
    
    def save_df(self):
        """Saves the main pandas DataFrame self.df to the pickle file."""
        self.df.to_pickle(self.df_file)    
    
    
    def main_strans(self, strans_levs):
        """Runs strans for the main element under study"""
        
        if self.blank_strans_lev in strans_levs:
            wx.MessageBox('STRANS input contains blank levels. \n\nPlease edit or delete these before running STRANS.', 'Blank Levels Found', 
                          wx.OK | wx.ICON_EXCLAMATION)
            self.frame_statusbar.SetStatusText('')
            return False
                     
        self.df['main_desig'] = np.empty((len(self.df), 0)).tolist()  # replaces any values in main_desig column with empty lists
        desig_list = self.df[['wavenumber', 'main_desig', 'line_tags']].values.tolist()
          
        matched_lines = self.strans(strans_levs, desig_list, self.main_element_name)          
        self.df.update(matched_lines)  
        self.display_strans_lines()

        return True
        
        
    def other_strans(self, other_lev_list):
        """Runs strans for all other elements that could be present in the linelist"""
                     
        self.df['other_desig'] = np.empty((len(self.df), 0)).tolist()  # replaces any values in other_desig column with empty lists
        desig_list = self.df[['wavenumber', 'other_desig', 'line_tags']].values.tolist()
        
        for other_lev in other_lev_list:
            element_name, level_file = other_lev.split(',')
            strans_levs = list(pd.read_csv(level_file).transpose().to_dict().values())
            matched_lines = self.strans(strans_levs, desig_list, element_name)

        self.df.update(matched_lines)  
        self.display_strans_lines()

  
    def strans(self, strans_levs, desig_list, element_name):
        """Creates list of all possible transitions between levels of opposite parity that obey
        the J selection rule. The list is then compared to all lines in the self.db database and lines with 
        wavenumbers that within self.strans_wn_discrim are assigned the labels of the even and odd level."""
                    
        self.frame_statusbar.SetStatusText(f'Running Strans for {element_name}')
                
        strans_levs_even = [x for x in strans_levs if x['parity']==1]
        strans_levs_odd = [x for x in strans_levs if x['parity']==0]
        
        strans_levs_even = sorted(strans_levs_even, key=lambda x: x['j'])
        strans_levs_odd = sorted(strans_levs_odd, key=lambda x: x['j'])
        
        for even_lev in strans_levs_even:    
            j_even = even_lev['j']
            
            if j_even == 0.0:  # J selection rule J != 0 to 0
                left_j = bisect.bisect_left(KeyList(strans_levs_odd, key=lambda x: x['j']), j_even + 1)
            else:  # the other J selection 
                left_j = bisect.bisect_left(KeyList(strans_levs_odd, key=lambda x: x['j']), j_even - 1)
                
            right_j = bisect.bisect_right(KeyList(strans_levs_odd, key=lambda x: x['j']), j_even + 1)
                            
            for odd_lev in strans_levs_odd[left_j:right_j]:                
                label_even = even_lev['label']
                label_odd = odd_lev['label']
                energy_even = even_lev['energy']
                energy_odd = odd_lev['energy']                   
                match_wn = abs(energy_even - energy_odd)
                    
                left = bisect.bisect_left(KeyList(desig_list, key=lambda x: x[0]), match_wn - self.strans_wn_discrim)
                right = bisect.bisect_left(KeyList(desig_list, key=lambda x: x[0]), match_wn + self.strans_wn_discrim)
                    
                for matched_line in desig_list[left:right]:
                    
                    if len(desig_list[left:right]) > 1:  # multiple lines match this transtion
                        matched_line[2]['multiple_lines'] = True
                    
                    if energy_even > energy_odd:  # assign upper and lower levels correctly (LOPT needs them in lower-upper format)
                        upper_lev = label_even
                        lower_lev = label_odd
                    else:
                        upper_lev = label_odd
                        lower_lev = label_even
                        
                    matched_line[1].append({'upper_level':upper_lev, 'lower_level':lower_lev, 'element_name': element_name})
        
        return desig_list
        

    def load_main_config(self):
        """Reads the main TAME config file and sets variables accordingly"""
        self.main_config_file = 'tame.ini'
        self.main_config = configparser.ConfigParser()
        self.main_config.read(self.main_config_file)
        
        self.strans_wn_discrim = self.main_config.getfloat('strans', 'wn_discrim')
        self.project_config_file = self.main_config.get('project', 'project_config')
        self.lopt_default_unc = self.main_config.getfloat('lopt', 'default_unc')
        
        
    def load_project_config(self):
        """Reads the project config file and sets variables accordingly"""
        self.project_config = configparser.ConfigParser()
        self.project_config.read(self.project_config_file)
        
        self.strans_lev_file = self.project_config.get('files', 'strans_lev_file')
        self.strans_lin_file = self.project_config.get('files', 'strans_lin_file')
        self.df_file = self.project_config.get('files', 'df_file')
        self.other_lev_list = self.project_config.get('files', 'other_lev_files').split('\n')

        self.lopt_fixed_levels = self.project_config.get('lopt', 'fixed_levels').split(',')
      
        self.main_element_name = self.project_config.get('tame', 'main_element_name').strip("'")        
        self.project_title = self.project_config.get('tame', 'project_title').strip("'")
    
    
    def load_df(self):
        """Loads the main lines df from file. Creates a new df if this is a new project."""
        if not os.path.isfile(self.df_file):  # if no existing DataFrame is present
            self.create_df(self.strans_lin_file)
            
        self.df = pd.read_pickle(self.df_file)  
    
    def load_project(self):
        """Set all filenames and variables and load/reload all listctrls."""
        self.load_project_config()
        
        self.lopt_inp_file = f'LOPT/{self.main_element_name}_lopt.inp'
        self.lopt_par_file = f'LOPT/{self.main_element_name}_lopt.par'
        self.lopt_fixed_file = f'LOPT/{self.main_element_name}_lopt.fixed'
        self.lopt_lev_file = f'LOPT/{self.main_element_name}_lopt.lev'
        self.lopt_lin_file = f'LOPT/{self.main_element_name}_lopt.lin'
        self.star_discrim = 1.5
          
        self.load_df()  
        self.load_plot_df()           
        self.strans_levs = list(pd.read_csv(self.strans_lev_file, dtype={'parity':float}).transpose().to_dict().values())                 
        self.display_strans_levs() 
        self.SetTitle(f"Term Analysis Made Easy (TAME) - {self.project_title}")
        
        
    def create_project_config(self):
        """Creates config file for a new project"""
        self.project_config = configparser.ConfigParser()
    
    def save_project(self):
        self.save_df()
        self.save_main_config()
        self.save_project_config()
        
    def save_project_config(self):
        with open(self.project_config_file, 'w') as configfile:
            self.project_config.write(configfile)
        
        
    def save_main_config(self):        
        with open(self.main_config_file, 'w') as configfile:
            self.main_config.write(configfile)
            
    def write_lopt_inp(self):
        with open(self.lopt_inp_file, 'w') as inp_file:
            lines = self.df.loc[self.df.main_desig.str.len() > 0 ].values.tolist()  # create list of all lines with a main designation

            if lines == []:  # no lines have a main_designation in self.df ie strans has not been run
                wx.MessageBox('No lines found for LOPT input. Please run STRANS first', 'No Matched Lines', 
                      wx.OK | wx.ICON_EXCLAMATION)
                return False
                
            else:
                for line in lines:
                    snr = f'{line[1]:9.0f}'
                    wn = f'{line[0]:15.4f}'
                    tag = '       '
                    main_desigs = line[6]
                    other_desigs = line[7] 
                    user_desig = line[8]
                    tags = line[9] 
                         
                    if user_desig != '':  # user label for line
                        upper_level = f'{user_desig["upper_level"]:>11}'
                        lower_level = f'{user_desig["lower_level"]:>11}'
                        
                        if all(value == False for value in tags.values()): # no user defined tags for the line
                            unc = f'{line[5]:.4f}'
                        elif tags['user_unc'] != False:
                            unc = f"{tags['user_unc']:.4f}"
                        else:
                            unc = f'{self.lopt_default_unc:.4f}'
                            tag = '       B'
                        
                        lopt_str = f'{snr}{wn} cm-1 {unc}{upper_level}{lower_level}{tag}\n'
                        inp_file.writelines(lopt_str)
                        
                    else:  # no user label for line
                        if len(main_desigs) != 1 or len(other_desigs) !=0:  # multiple identifications for line
                            unc = f'{self.lopt_default_unc:.4f}'
                        elif all(value == False for value in tags.values()): # no user defined tags for the line
                            unc = f'{line[5]:.4f}'
                        elif tags['user_unc'] != False:
                            unc = f"{tags['user_unc']:.4f}"
                        else:
                            unc = f'{self.lopt_default_unc:.4f}'
                            tag = '       B'
                            
                        for desig in main_desigs:
                            upper_level = f'{desig["upper_level"]:>11}'
                            lower_level = f'{desig["lower_level"]:>11}'
                        
                            lopt_str = f'{snr}{wn} cm-1 {unc}{lower_level}{upper_level}{tag}\n'
                            inp_file.writelines(lopt_str)
        return True
        
        
    def write_lopt_par(self):      
        with open('LOPT/lopt_template.par', 'r') as temp_par_file:
            par_lines = temp_par_file.readlines()
            par_lines[0] = f'{self.main_element_name}_lopt.inp{par_lines[0]}'
            par_lines[1] = f'{self.main_element_name}_lopt.fixed{par_lines[1]}'
            par_lines[2] = f'{self.main_element_name}_lopt.lev{par_lines[2]}'
            par_lines[3] = f'{self.main_element_name}_lopt.lin{par_lines[3]}'
            
        with open(self.lopt_par_file, 'w') as par_file:
            par_file.writelines(par_lines)
            
    
    def write_lopt_fixed(self):
        with open(self.lopt_fixed_file, 'w') as fixed_file:
            fixed_strings = []
            
            for level in self.lopt_fixed_levels:
                strans_lev = next((item for item in self.strans_levs if item['label']==level))
                lev_energy = strans_lev['energy'] 
                if lev_energy == 0.0:
                    lev_unc = f'{0.0:.4f}'
                else:
                    lev_unc = f'{2.0:.4f}'
                             
                fixed_strings.append(f'{level:>9}{lev_energy:>9.2f}{lev_unc:>13}\n')
            
            fixed_file.writelines(fixed_strings)    
            
    def get_lopt_output(self):
        self.lopt_levs = pd.read_csv(self.lopt_lev_file, delimiter='\t')
        # print(self.lopt_levs.head)
        #lopt_levs = list(lopt_levs.transpose().to_dict().values())
        lopt_lines_df = pd.read_csv(self.lopt_lin_file, delimiter='\t')
        # print(lopt_lines_df.columns)
        merged_lines = pd.merge_asof(lopt_lines_df[['W_obs', 'S', 'Wn_c', 'E1', 'E2', 'L1', 'L2', 'F', 'uncW_o']].sort_values('W_obs'), 
                                     self.df[['wavenumber', 'peak', 'eq width', 'tags']].sort_values('wavenumber'), 
                                     left_on='W_obs', 
                                     right_on='wavenumber',
                                     tolerance=0.005,
                                     direction='nearest') # match lopt lines to main df file based on nearest wavenumber
        
        merged_lines['dWO-C'] = merged_lines['W_obs'] - merged_lines['Wn_c']
        merged_lines['star'] = np.where(((merged_lines['dWO-C'].abs()*self.star_discrim) > merged_lines['uncW_o']), True, False)
        merged_lines['log_ew'] = np.log(merged_lines['eq width'])
        merged_lines['main_level'] = ''
        merged_lines['other_level'] = ''
        merged_lines['F'] = np.where(merged_lines['F'] == None, ' ', merged_lines['F'])
                
        duplicated_lines = pd.DataFrame(np.repeat(merged_lines.copy().values,2,axis=0))
        duplicated_lines.columns = merged_lines.columns  
                
        t1 = duplicated_lines.iloc[0::2].copy()
        t1['main_level'] = t1['E1']
        t1['other_level'] = t1['L2'] 
        t2 = duplicated_lines.iloc[1::2].copy()
        t2['main_level'] = t2['E2'] 
        t2['other_level'] = t2['L1'] 
        
        # t1 = duplicated_lines.iloc[0::2].copy()
        # t1['main_level'] = t1['L1']
        # t1['other_level'] = t1['L2'] 
        # t2 = duplicated_lines.iloc[1::2].copy()
        # t2['main_level'] = t2['L2'] 
        # t2['other_level'] = t2['L1'] 
        
        duplicated_lines.update(t1)
        duplicated_lines.update(t2)
           
        # duplicated_lines.iloc[0::2]['transition'] = duplicated_lines['L1']
        # duplicated_lines.iloc[1::2]['transition'] = duplicated_lines['L2']
             
        duplicated_lines = list(duplicated_lines.transpose().to_dict().values())
        # print(duplicated_lines)
        self.lopt_lev_ojlv.SetObjects(duplicated_lines) 
        
        self.load_lopt_lev_comments()  # needs to be done properly
        
        
    def load_lopt_lev_comments(self):
        self.lopt_lev_comments_file = 'lopt_lev.pkl'
        if not os.path.isfile(self.lopt_lev_comments_file):  # if no existing DataFrame is present
            self.lopt_lev_comments = self.lopt_levs[['Designation']].copy()
            self.lopt_lev_comments['Comments'] = ''        
            self.lopt_lev_comments.to_pickle(self.lopt_lev_comments_file) 
        else:
            lopt_lev_com_tmp = self.lopt_levs[['Designation']].copy()
            lopt_lev_com_tmp['Comments'] = ''  
            self.lopt_lev_comments = pd.read_pickle(self.lopt_lev_comments_file) 
            self.lopt_lev_comments.combine_first(lopt_lev_com_tmp)  # should combine the two dataframes and keep the values from lopt_lev_comments

        
    def display_lopt_line(self, line):  
        line_dict = list(line.transpose().to_dict().values()).pop()   
        self.user_unc_txtctrl.SetValue('')  # sets back to empty when new line selected        
        self.lopt_line_panel_header.SetLabel(f"Line: {line_dict['wavenumber']:.4f} {self.cm_1}")
        
        ### Get all designations and update the desig listctrl ###
        self.lopt_line_listctrl.DeleteAllItems()
        
        for i, desig in enumerate(line_dict['main_desig'] + line_dict['other_desig']):
            list_ctrl_list = [desig['element_name'], 
                              desig['upper_level'], 
                              desig['lower_level'],
                              'wn_diff']
            
            self.lopt_line_listctrl.Append(list_ctrl_list)
            
            if desig == line_dict['user_desig']:
                self.lopt_line_listctrl.CheckItem(i, True)
        
        ### Update the checkboxes ###   
        line_tags = line_dict['line_tags']        
        self.ringing_chkbox.SetValue(line_tags['ringing'])
        self.noise_chkbox.SetValue(line_tags['noise'])
        self.blend_chkbox.SetValue(line_tags['blend'])
        
        if line_tags['user_unc'] != False:  # if not False
            self.user_unc_chkbox.SetValue(True)
            self.user_unc_txtctrl.SetValue(str(line_tags['user_unc']))
        else:
            self.user_unc_chkbox.SetValue(False)
            self.user_unc_txtctrl.SetValue('')
            
        with wx.EventBlocker(self):  # prevents the textctrl event from firing and writing this value to the df
            self.lopt_line_comments_txtctrl.SetValue(line_dict['comments'])  
        
        ### Plot spectra ###        
        self.plot_lopt_line(line_dict['wavenumber'], line_dict['peak'], line_dict['width'])
    
        #self.lopt_line_panel.SetPosition((0,0))
        self.lopt_level_panel.Hide()
        self.lopt_line_panel.Show()
        self.sizer_8.Layout()
     
    def plot_lopt_line(self, wavenumber, peak, width):   
        """Plots the spectra from self.plot_df to the matplotlib plot. There are scaling factors that can be changed to show
        more/less of the plot area."""
        plot_width = 1.
        # plot_y_scale = 1.1
        plot_x_scale = 5.
        
        self.matplotlib_canvas.clear()  
        ax = self.matplotlib_canvas.gca()   
        # df_part = self.plot_df.loc[(self.plot_df['wavenumber'] < wavenumber+plot_width) & (self.plot_df['wavenumber'] > wavenumber-plot_width)]
        #Test to see if below or above is faster
        df_part_1 = self.plot_df.loc[(self.plot_df['wavenumber'] < wavenumber+plot_width)]
        df_part = df_part_1.loc[(df_part_1['wavenumber'] > wavenumber-plot_width)]
        
        
        spectras = df_part.columns[df_part.notna().any()].tolist()[1:]  # gives columns that do not contain only NaN.
        
        for spectra in spectras:
            df_part.plot(kind='line', x='wavenumber', y=spectra, ax=ax)
            
        self.matplotlib_canvas.axes.set_xlabel('Wavenumber (cm-1)')
        self.matplotlib_canvas.axes.set_ylabel('SNR')
        ax.ticklabel_format(useOffset=False)
        # ax.set_ylim([df_part.min().min()*plot_y_scale, peak*plot_y_scale])
        ax.set_xlim([wavenumber - (width/1000)*plot_x_scale, wavenumber + (width/1000)*plot_x_scale])
        # ax.axvline(x=wavenumber)  # due to calibration this is not in the right place and looks odd

        self.matplotlib_canvas.draw()     
    
    def display_lopt_lev(self, level):
        lev_dict = list(level.transpose().to_dict().values()).pop()         
        self.lopt_lev_panel_header.SetLabel(f"Level: {lev_dict['Designation']}")
        self.lopt_level_listctrl.DeleteAllItems()
        
        if lev_dict['Comments'] is np.nan:
            lopt_comments = ''
        else:
            lopt_comments = lev_dict['Comments']
        
        lopt_level_list = [f"   {lev_dict['Energy']:.4f}",
                           lev_dict['D1'],
                           lev_dict['D2'],
                           lev_dict['D3'],
                           lev_dict['N_lines'],
                           lopt_comments]        
        self.lopt_level_listctrl.Append(lopt_level_list)   
         
        next_row = self.lopt_lev_ojlv.GetObjectAt(self.lopt_lev_ojlv.GetFocusedRow()+1)  # gives the first line row in the group      
        selected_lev = self.group_column.GetGroupKeyAsString(self.group_column.GetGroupKey(next_row))  # gives the group title of the selected level, i.e. the level designation
        self.selected_lev_index = self.lopt_lev_comments.loc[self.lopt_lev_comments['Designation'] == selected_lev].index.values[0]  # get index of the level     
        
        with wx.EventBlocker(self):  # prevents the textctrl event from firing and writing this value to the df
            self.lopt_level_comments.SetValue(self.lopt_lev_comments.at[self.selected_lev_index, 'Comments']) 

        
        self.lopt_level_panel.Show()
        self.lopt_line_panel.Hide()
        self.sizer_8.Layout()
        
    def update_df_cell(self, wavenumber, column, value):
        """Updates a cell of the main self.df dataframe, specified by wavenumber and column, with value"""   
        selected_line_index = self.df.loc[self.df['wavenumber'] == wavenumber].index.values[0]
        self.df.at[selected_line_index, column] = value  # just updates a single value
        print(self.df.iloc[selected_line_index])
        
    def get_df_cell(self, wavenumber, column):
        selected_line_index = self.df.loc[self.df['wavenumber'] == wavenumber].index.values[0]
        return self.df.at[selected_line_index, column]
    

    def search_listview(self, event, listview):
        """Searches the primary column of the listview with the string typed into the search ctrl. If the entered
        text is not found, the serachctrl background turns red."""
        search_bar = event.GetEventObject()
        search_str = search_bar.GetValue()         

        if not listview._FindByTyping(listview.GetPrimaryColumn(), search_str):
            search_bar.SetBackgroundColour(wx.RED)
        else:
            search_bar.SetBackgroundColour(wx.WHITE)
            
    def is_float(self, value):
        """Checks if a given value is a number or not"""
        try:
            float(value)
            return True
        except ValueError:
            return False
        
        

### Event-driven functions ###  

    def on_lopt_lev_comments(self, event):
        text = self.lopt_level_comments.GetValue()   
        self.lopt_lev_comments.at[self.selected_lev_index, 'Comments'] = text  # updates level with the comments

    def on_lopt_line_comments(self, event):
        text = self.lopt_line_comments_txtctrl.GetValue()
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']
        self.update_df_cell(selected_wn, 'comments', text)

    def on_ringing_tag(self, event):      
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']         
        line_tags = self.get_df_cell(selected_wn, 'line_tags')
        line_tags['ringing'] = self.ringing_chkbox.GetValue()
        self.update_df_cell(selected_wn, 'line_tags', line_tags)

    def on_noise_tag(self, event):  
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']         
        line_tags = self.get_df_cell(selected_wn, 'line_tags')
        line_tags['noise'] = self.noise_chkbox.GetValue()
        self.update_df_cell(selected_wn, 'line_tags', line_tags)

    def on_blend_tag(self, event):  
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']         
        line_tags = self.get_df_cell(selected_wn, 'line_tags')
        line_tags['blend'] = self.blend_chkbox.GetValue()
        self.update_df_cell(selected_wn, 'line_tags', line_tags)

    def on_user_unc_tag(self, event): 
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']         
        line_tags = self.get_df_cell(selected_wn, 'line_tags')
        
        if self.user_unc_chkbox.GetValue():
            user_unc = self.user_unc_txtctrl.GetLineText(0)
            
            if self.is_float(user_unc) and float(user_unc) > 0.0 and float(user_unc) < 10.0:
                line_tags['user_unc'] = float(user_unc)
            else:
                wx.MessageBox('User uncertainty must be a number between 0.0 and 9.9999', 'Incorrect Uncertainty', 
                      wx.OK | wx.ICON_EXCLAMATION)
                self.user_unc_chkbox.SetValue(False)
                self.user_unc_txtctrl.SetValue('')
                return            
        else:
            line_tags['user_unc'] = False
        
        self.update_df_cell(selected_wn, 'line_tags', line_tags)


    def on_lopt_trans_checked(self, event): 
        user_desig = {}        
        line_index = event.GetIndex()
        
        user_desig['upper_level'] = self.lopt_line_listctrl.GetItem(line_index, 1).GetText()
        user_desig['lower_level'] = self.lopt_line_listctrl.GetItem(line_index, 2).GetText()
        user_desig['element_name'] = self.lopt_line_listctrl.GetItem(line_index, 0).GetText()        
        
        with wx.EventBlocker(self):  # stops the unchecked event from firing  
            for i in range(self.lopt_line_listctrl.GetItemCount()):
                if i != line_index:
                    self.lopt_line_listctrl.CheckItem(i, False)
                    
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']       
        self.update_df_cell(selected_wn, 'user_desig', user_desig)

    def on_lopt_trans_unchecked(self, event):
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']       
        self.update_df_cell(selected_wn, 'user_desig', '')

    def on_partial_strans(self, event):  # Run >> Strans(partial)
        if self.main_strans(self.strans_levs):  # if strans ran succesfully
            self.main_panel.ChangeSelection(0)  # changes the notebook tab to LOPT
            self.frame_statusbar.SetStatusText('Strans Complete')
        
    def on_full_strans(self, event):  
        if self.main_strans(self.strans_levs): # if strans ran succesfully
            self.other_strans(self.other_lev_list)
            self.main_panel.ChangeSelection(0)  # changes the notebook tab to LOPT
            self.frame_statusbar.SetStatusText('Strans Complete')
           
    def on_strans_del(self, event):  
        """Delete levels from self.strans_lev_file. This will be a permanent change."""
        selected_levs = self.strans_lev_ojlv.GetSelectedObjects()
        
        if selected_levs:  # if selection not empty
            if len(selected_levs) == 1:
                message = 'Are you sure you want to delete this level?'
                title = 'Delete Level?'
            else:
                message = f'Are you sure you want to delete these {len(selected_levs)} levels?'
                title = 'Delete Levels?'
        
            if wx.MessageBox(message, title, 
                          wx.YES_NO | wx.NO_DEFAULT | wx.ICON_EXCLAMATION) == wx.YES:
                
    
                self.strans_levs = [x for x in self.strans_levs if x not in selected_levs]            
                self.display_strans_levs()
                
                # XXX will need to add a corresponding part in save to write these values back out to the strans file.               
            else:
                return
            
    def on_strans_add(self, event):  
        next_row = 0

        self.strans_levs.insert(next_row, self.blank_strans_lev)
        self.display_strans_levs()
        
    def on_strans_save(self, event):  
        print("Event handler 'on_strans_save' not implemented!")
        event.Skip()
        
    def on_strans_lev_edit(self, event):  
        self.strans_levs = self.strans_lev_ojlv.GetObjects()  # updates the edited cells to strans_levs

            
    def on_Export_Linelist(self, event):  
        with wx.FileDialog(self, "Export Matched Linelist", wildcard="Linelist Files (*.lin)|*.lin",
                       style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            filename = fileDialog.GetPath()
                
            # try:
            with open(filename, 'w') as file:
                    lines = self.df.loc[self.df.main_desig.str.len() > 0 ].values.tolist()
                    file.writelines('wavenumber,snr,fwhm,eq_width,fit,unc,main_element,other_elements\n') 
                    
                    for line in lines:
                        main_level_string = ''
                        other_level_string = ''
                        for lev_pair in line[6]:
                            if not main_level_string:
                                sep = ''
                            else:
                                sep = ',     '
                            main_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['upper_level'] + ' - ' + lev_pair['lower_level']
                            
                        for lev_pair in line[7]:
                            if not other_level_string:
                                sep = ''
                            else:
                                sep = ',     '
                            other_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['upper_level'] + ' - ' + lev_pair['lower_level']
                            
                        
                        line[0] = f'{line[0]:.4f}'
                        line[1] = f'{line[1]:.0f}' 
                        line[2] = f'{line[2]:.0f}' 
                        line[3] = f'{line[3]:.0f}' 
                        line[5] = f'{line[5]:.4f}' 
                        line[6] = main_level_string
                        line[7] = other_level_string
                        
                        file.writelines(','.join(line[:7]) + '\n')   
                        self.frame_statusbar.SetStatusText(f'Matched linelist exported to {filename}')       
                
            # except:
            #     wx.MessageBox(f'Unsupported .ini file. Please select a TAME project file.', 'Unsupported File', 
            #           wx.OK | wx.ICON_EXCLAMATION)
                
    def on_Open(self, event):  # File >> Open
        save_dlg = wx.MessageBox(f'Do you want to save changes to {self.project_title}', 'Save Changes?', 
                                 wx.YES_NO | wx.CANCEL | wx.CANCEL_DEFAULT | wx.ICON_INFORMATION)
        if save_dlg == wx.ID_YES:
            self.save_project()
        elif save_dlg == wx.ID_CANCEL:
            return

        
        with wx.FileDialog(self, "Open TAME project file", wildcard="project files (*.ini)|*.ini",
                       style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            self.project_config_file = fileDialog.GetPath()
            try:
                self.load_project()
                self.main_config.set('project', 'project_config', self.project_config_file)
                self.save_main_config()
                self.frame_statusbar.SetStatusText('Project loaded successfully')
                
            except:
                wx.MessageBox('Unsupported or out of date .ini file. Please select a valid TAME project file.', 'Unsupported File', 
                      wx.OK | wx.ICON_EXCLAMATION)
     
                
                
    def set_fixed_levels(self):
        fixed_level_dialog = fixedLevels(self)
        
        if fixed_level_dialog.ShowModal() == wx.ID_CANCEL:
            return False
        
        self.lopt_fixed_levels = fixed_level_dialog.fixed_levels     
        
        level_config_string = ','.join(self.lopt_fixed_levels)
        print(level_config_string)
        self.project_config.set('lopt', 'fixed_levels', level_config_string)
        
        return True
   
    
    def on_lopt(self, event):
        """Create/update all neccessary files for LOPT input and then call LOPT"""
        
        if self.lopt_fixed_levels == ['']:  # no fixed levels have been selected
            if not self.set_fixed_levels():  # if cancel button pressed on fixed level dialog
                return
        
        self.frame_statusbar.SetStatusText('Writing LOPT input files...')
        if self.write_lopt_inp():  # lines in the strans linelist
            self.write_lopt_par()
            self.write_lopt_fixed()
            
            self.frame_statusbar.SetStatusText('Running LOPT...')    
            try:
                p = subprocess.run(['java', '-jar', 'Lopt.jar', self.lopt_par_file.split('/')[-1]], cwd='LOPT/', capture_output=True, text=True).stdout.split('\n')  # run LOPT and get output as a list of lines
                
                rss = [x for x in p if 'RSS' in x]  # gives the RSS\degrees_of_freedom line
                tot_time = [x for x in p if 'Total time' in x]  # gives the total time line
                self.frame_statusbar.SetStatusText(f'LOPT ran successfully:  {rss[0]}. {tot_time[0]}.')  
                
                self.get_lopt_output()
                self.main_panel.ChangeSelection(1)  # changes the notebook tab to LOPT
            
            except FileNotFoundError as fnf:
                print(fnf)
                if 'java' in str(fnf):
                    self.frame_statusbar.SetStatusText('LOPT error') 
                    wx.MessageBox('Java Runtime Environment (JRE) is not installed on this machine. \n\nPlease install and restart TAME.', 'Missing Java runtime', 
                          wx.OK | wx.ICON_EXCLAMATION)
            # except IndexError:
            #     print('v')
            #     self.frame_statusbar.SetStatusText('Run STRANS first') 
        else:
            self.frame_statusbar.SetStatusText('Run STRANS first') 
                        
    def on_Save(self, event):  
        self.save_project()
        self.frame_statusbar.SetStatusText('Project saved')                   
        
    def on_Exit(self, event):  # File >> Exit
        save_dlg = wx.MessageBox(f'Do you want to save changes to {self.project_title}', 'Save Changes?', 
                                 wx.YES_NO | wx.CANCEL | wx.CANCEL_DEFAULT | wx.ICON_INFORMATION)
        if save_dlg == wx.YES:
            self.save_project()
            self.Destroy()
        elif save_dlg == wx.NO:
            self.Destroy() 
        else:
            return       
        
    def on_click_lopt_levs(self, event):  
        try:
            selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']
            selected_line = self.df.loc[self.df['wavenumber'] == selected_wn]  
            
            self.display_lopt_line(selected_line)

        except TypeError: # group header or blank row selected
            next_row = self.lopt_lev_ojlv.GetObjectAt(self.lopt_lev_ojlv.GetFocusedRow()+1)
            
            if not next_row == None:  # group header selected and not the blank line above it         
                next_row_lev = next_row['main_level']
                # print(next_row_lev)
                selected_lev = self.lopt_levs.loc[self.lopt_levs['Energy'] == next_row_lev]
                
                self.display_lopt_lev(selected_lev)    
            else:
                self.lopt_level_panel.Hide()
                self.lopt_line_panel.Hide()
                self.sizer_8.Layout()
                
        
            event.Skip()    
        
    def on_col_clicked(self, event):  # wxGlade: mainWindow.<event_handler>
        #self.lopt_lev_ojlv.SortBy(1)
        # self.lopt_lev_ojlv.RebuildGroups()
        event.Skip()
        
    def on_strans_lev_search(self, event):
        self.search_listview(event, self.strans_lev_ojlv)
        
    def on_strans_line_search(self, event):  
        search_bar = event.GetEventObject()
        try:
            self.search_listview(event, self.strans_lines_ojlv)
            search_bar.SetBackgroundColour(wx.WHITE)
        except:
            search_bar.Clear()
            search_bar.SetBackgroundColour(wx.RED)
            
    def on_New(self, event):  
        self.new_proj = newProject(self)
        self.new_proj.ShowModal()  
        
        
        try:
            self.template_config_file = 'template.ini'
            self.new_config = configparser.ConfigParser()
            self.new_config.read(self.template_config_file)
            
            self.new_config.set('files', 'strans_lev_file', self.new_proj.main_element_lev_file)
            self.new_config.set('files', 'strans_lin_file', self.new_proj.linelist_file)
            self.new_config.set('files', 'df_file', self.new_proj.df_file)
            
            other_lev_files = ''
            for file in self.new_proj.other_element_lev_files:
                other_lev_files += f'{file["shortname"]},{file["filename"]}\n'
                       
            self.new_config.set('files', 'other_lev_files', other_lev_files)
            
            self.new_config.set('tame', 'project_title', self.new_proj.project_name)
            self.new_config.set('tame', 'main_element_name', self.new_proj.main_element_name)
            
            self.new_config_file = self.new_proj.project_name + '.ini'
            
            with open(self.new_config_file, 'w') as configfile:
                self.new_config.write(configfile)
                
        except AttributeError:  # if the user cancelled the new project process.
            event.Skip()
  
        
  
    
  
    
  
class fixedLevels(fixedLevelsDialog):
    def __init__(self, *args, **kwds):
        fixedLevelsDialog.__init__(self, *args, **kwds)  
        self.fixed_levels = []
        self.fixed_level_lc.EnableCheckBoxes(True)
        self.fixed_level_lc.DeleteAllItems()
                
        for level in self.GetParent().strans_levs:
            self.fixed_level_lc.Append([level['label'], f"{level['energy']:.4f}"])
            
    def on_fixed_lev_checked(self, event):
        level = self.fixed_level_lc.GetItem(event.GetIndex(), 0).GetText()
        self.fixed_levels.append(level)
               
    def on_fixed_lev_unchecked(self, event):
        level = self.fixed_level_lc.GetItem(event.GetIndex(), 0).GetText()
        self.fixed_levels = [x for x in self.fixed_levels if x != level]
    
              
class newProject(newProjectDialog):
    def __init__(self, *args, **kwds):
        newProjectDialog.__init__(self, *args, **kwds)
        self.button_9.Disable()  # greys out back button on first window
        self.other_element_lev_files = []
        
        self.config_file = 'TEST'
        
        
        self.other_element_levs_ojlv.SetColumns([
            ColumnDefn("Filename", "left", 500, 'filename'),          
            ColumnDefn("Element Short Name", "left", 50, 'shortname', isSpaceFilling=True, isEditable=True)])
        
        self.other_element_levs_ojlv.SetEmptyListMsg('')
        self.other_element_levs_ojlv.cellEditMode = self.other_element_levs_ojlv.CELLEDIT_SINGLECLICK
        
        
    def on_main_lev_file_btn(self, event):
        with wx.FileDialog(self, "Select level file", wildcard="level files (*.lev)|*.lev",
                       style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            self.main_element_lev_file = fileDialog.GetPath()
            self.main_lev_file_tc.SetValue(str(self.main_element_lev_file))

    def on_linelist_file_btn(self, event):  # wxGlade: NewProjectFrame.<event_handler>
        with wx.FileDialog(self, "Select linelist file", wildcard="linelist files (*.lin)|*.lin",
                       style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            self.linelist_file = fileDialog.GetPath()
            self.linelist_file_tc.SetValue(str(self.linelist_file))

    def on_plot_files_btn(self, event):  # wxGlade: NewProjectFrame.<event_handler>
        with wx.FileDialog(self, "Select linelist file", wildcard="plot files (*.asc)|*.asc",
                       style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            self.plot_files = fileDialog.GetPaths()
            for file in self.plot_files:
                self.plot_file_tc.write(str(file) + '\n')

    def on_cancel(self, event):
        self.Destroy()

    def on_next(self, event):
        self.project_name = self.project_name_tc.GetValue()
        self.df_file = self.project_name.split('.')[0] + '.pkl'
        self.main_element_name = self.element_name_tc.GetValue()
        
        if not self.project_name or not self.main_element_name or not self.linelist_file or not self.plot_files:
            wx.MessageBox('Please fill in all fields before continuing', 'Missing Project Parameters', 
                      wx.OK | wx.ICON_EXCLAMATION)
        else:
            self.new_project_nb.ChangeSelection(1)

    def on_other_lev_files_btn(self, event):
        with wx.FileDialog(self, "Select linelist file", wildcard="level files (*.lev)|*.lev",
                       style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            lev_files = fileDialog.GetPaths()
            
            for file in lev_files:
                self.other_lev_files_tc.write(str(file) + '\n')                
                self.other_element_lev_files.append({'filename':file, 'shortname':''})
                
            self.other_element_levs_ojlv.AddObjects(self.other_element_lev_files)

    def on_back(self, event):
        self.new_project_nb.ChangeSelection(0)

    def on_ok(self, event):
        
        if not self.other_element_lev_files:
            wx.MessageBox('Please fill in all fields before continuing', 'Missing Project Parameters', 
                      wx.OK | wx.ICON_EXCLAMATION)     
        else:
            
            if '' in [x['shortname'] for x in self.other_element_lev_files]:
                wx.MessageBox('Please fill in short names for all elements before continuing', 'Missing Project Parameters', 
                      wx.OK | wx.ICON_EXCLAMATION)
            else:              
                # SAVE all parameters to a new config file.
                self.Destroy()
        
        
    

            
        
class KeyList(object):
    """Key function for bisect search. This allows bisect to work on specific elements within lists of iterable objects"""
    def __init__(self, l, key):
        self.l = l
        self.key = key
    def __len__(self):
        return len(self.l)
    def __getitem__(self, index):
        return self.key(self.l[index])


class MyApp(wx.App):
    def OnInit(self):
        self.frame = MyFrame(None, wx.ID_ANY, "")
        # self.SetTopWindow(self.frame)
        self.frame.Show()
        return True


if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()