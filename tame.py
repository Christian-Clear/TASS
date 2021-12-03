#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import wx
import numpy as np
import pandas as pd
import os
from lib.tame_gui import mainWindow, newProjectDialog, fixedLevelsDialog, propertiesDialog, lostLinesDialog
import os.path
import bisect
import configparser
import subprocess
from lib.ObjectListView import ColumnDefn, OLVEvent
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from shutil import copy

import warnings  # only here to stop deprecation warning of objectlistview from clogging up terminal
warnings.filterwarnings("ignore", category=DeprecationWarning)

TAME_VERSION_STRING = 'Version: 0.0.1'

class MyFrame(mainWindow):
    """The main TAME window. Everything is run within this class."""
    def __init__(self, *args, **kwds):
        """Initialise mainWindow, set constants etc. and then load the configs for tame."""
        mainWindow.__init__(self, *args, **kwds)
        
        self.set_constants()
        self.configure_layout()
        self.configure_listviews()
        
        try:
            self.load_main_config()
            self.load_project()
            self.project_loaded = True          
        except:
            self.project_loaded = False
        
    def set_constants(self):
        """Set the constants that are used throughout TAME."""
        self.cm_1 = 'cm\u207B\u00B9'  # unicode for inverse centimetres
        self.blank_strans_lev = {'label': '', 'j':0.0 , 'energy':0.0 , 'parity':0}
        self.levhams_selected_levs = {}
        self.groupHeaderColour = wx.Colour(159, 185, 250, 249)  # BLUE
        self.evenRowsBackColour = wx.Colour(240, 248, 255)  # ALICE BLUE
        self.oddRowsBackColour = wx.Colour(255, 250, 205)  # LEMON CHIFFON
    
   
    def configure_layout(self):
        """ Configure the main TAME screen layout and positioning."""
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
        """Set column definitions and settings for all object and group listviews and basic wx.listctrls."""
        self.strans_lev_ojlv.SetColumns([
            ColumnDefn("Level", "left", 100, 'label', isEditable=True),
            ColumnDefn("J", "left", 50, 'j', stringConverter="%.1f", isEditable=True),
            ColumnDefn(f"Energy ({self.cm_1})", "left", 120, 'energy', stringConverter="%.4f", isEditable=True),           
            ColumnDefn("Parity", "left", 50, 'parity', stringConverter="%d", isSpaceFilling=True, isEditable=True)])
    
        self.strans_lev_ojlv.cellEditMode = self.strans_lev_ojlv.CELLEDIT_DOUBLECLICK
        self.strans_lev_ojlv.Bind(OLVEvent.EVT_CELL_EDIT_FINISHED, self.on_finish_strans_lev_edit)
        self.strans_lev_ojlv.Bind(OLVEvent.EVT_CELL_EDIT_STARTING, self.on_start_strans_lev_edit)  
        self.strans_lev_ojlv.Bind(OLVEvent.EVT_SORT, self.on_strans_lev_sort)
        self.strans_lev_ojlv.SetEmptyListMsg("")
                
        self.strans_lines_ojlv.SetColumns([
            ColumnDefn(f'Wavenumber ({self.cm_1})', 'left', 150, 'wavenumber', stringConverter="%.4f"),
            ColumnDefn('SNR', 'left', 40, 'peak', stringConverter="%d"),
            ColumnDefn('FHWM (mK)', 'left', 90, 'width', stringConverter="%d"),
            ColumnDefn('log(Eq. Width)', 'left', 110, 'eq width', stringConverter="%.2f"),
            ColumnDefn('Fit', 'left', 30, 'tags'),
            ColumnDefn(f'Unc. ({self.cm_1})', 'left', 90, 'unc', stringConverter="%.4f"),
            ColumnDefn('Main Element Transitions', 'left', 500, 'main_desig'),
            ColumnDefn('Other Element Transitions', 'left', 500, 'other_desig', isSpaceFilling=True)])
        
        self.strans_lines_ojlv.SetEmptyListMsg("Run Line Matching")
        
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
        
        self.lopt_lev_ojlv.SetEmptyListMsg("Run Level Optimisation first")
        self.lopt_lev_ojlv.SetShowItemCounts(False)
        self.lopt_lev_ojlv.SetAlwaysGroupByColumn(1)  # Level energy column  
        self.lopt_lev_ojlv.Bind(OLVEvent.EVT_SORT, self.on_lopt_lev_sort)
        
        
        self.levhams_output_ojlv.SetColumns([
            ColumnDefn(f'Level Energy ({self.cm_1})', 'left', 150, 'energy', stringConverter=self.levhams_float_converter_4),
            ColumnDefn(f'Line Wn ({self.cm_1})', 'left', 130, 'wavenumber', stringConverter=self.levhams_float_converter_4),
            ColumnDefn('SNR', 'left', 50, 'peak', stringConverter=self.levhams_float_converter_0),
            ColumnDefn('Intensity', 'left', 70, 'eq width', stringConverter=self.levhams_log_ew_converter),
            ColumnDefn('Level', 'left', 100, 'level'),
            ColumnDefn('J', 'left', 40, 'j'),
            ColumnDefn(f'Predicted Energy ({self.cm_1})', 'left', 170, 'pred_energy', stringConverter=self.levhams_float_converter_4),
            ColumnDefn(f'Average Energy ({self.cm_1})', 'left', 170, 'avg_energy', stringConverter=self.levhams_float_converter_4),
            ColumnDefn(f'Separation ({self.cm_1})', 'left', 170, 'sep', stringConverter=self.levhams_float_converter_4)])        
                
        self.lopt_line_listctrl.EnableCheckBoxes(True)
        self.levhams_level_listctrl.EnableCheckBoxes(True)
        
        
        

    ### String/Group Converters for Object and Group Listviews #################
    def levhams_float_converter_4(self, num):
        try:
            return f'{float(num):.4f}'
        except:
            return ''
        
    def levhams_float_converter_0(self, num):
        try:
            return f'{float(num):.0f}'
        except:
            return ''
    
    def loptGroupKeyConverter(self, energy):
        """Convert energy of group to the level designation."""
        selected_line_index = self.lopt_levs.loc[self.lopt_levs['Energy'] == energy].index.values[0]
        return self.lopt_levs.at[selected_line_index, 'Designation']
    
    def log_ew_converter(self, log_ew):
        """Convert equivalent width (from .lin file) to log(ew)."""  # XXX Check that this is right? log_ew to log(ew)?
        if str(log_ew) == 'nan':  # Correct formatting if this is a virtual line from LOPT.
            return '-'
        return f'{log_ew:.2f}'
    
    def levhams_log_ew_converter(self, log_ew):
        """Convert equivalent width (from .lin file) to log(ew)."""  # XXX Check that this is right? log_ew to log(ew)?
        try:
            return f'{np.log(log_ew):.2f}'
        except:
            return ''
    
    def wn_converter(self, wn):
        """Convert wn to 4 decimal places."""
        if str(wn) == 'nan':  # Correct formatting if this is a virtual line from LOPT.
            return '-'
        return f'{wn:.4f}'
    
    def snr_converter(self, snr):
        """Convert SNR to integer."""
        if str(snr) == 'nan':  # Correct formatting if this is a virtual line from LOPT.
            return '-'
        return f'{snr:.0f}'
    
    def fit_converter(self, fit):
        """Correct formatting for virtual lines."""
        if str(fit) == 'nan':  # Correct formatting if this is a virtual line from LOPT.
            return '-'
        return fit
    
    def star_to_asterix(self, star):
        """Convert boolean of star to asterix or empty string."""
        if star:
            return '*'
        return ''
        
    def correct_tags(self, tag):
        """Correct formatting for virtual lines and line tags"""
        if str(tag) == 'nan':  # Correct formatting if this is a virtual line from LOPT.
            return ''
        elif tag == 'b':
            return 'U'
        elif tag == 'Q':
            return 'M'
        else:
            return tag.upper()
    
    def neg_num_str(self, diff):
        """Convert diff to 4 d.p. float and make it right alinged. Means that negative numbers aren't shifted to the right."""
        return f'{diff:>7.4f}'
    
    ############################################################################
    
    def load_main_config(self):
        """Reads the main TAME config file and sets variables accordingly"""
        self.main_config_file = 'config/tame.ini'
        self.main_config = configparser.ConfigParser()
        self.main_config.read(self.main_config_file)
        
        self.project_config_file = self.main_config.get('project', 'project_config')
        self.lopt_default_unc = self.main_config.getfloat('lopt', 'default_unc')
         
    def load_project_config(self):
        """Reads the project config file and sets variables accordingly"""
        self.project_config = configparser.ConfigParser()
        self.project_config.read(self.project_config_file)

        self.strans_lev_file = self.project_config.get('files', 'strans_lev_file')
        self.strans_lin_file = self.project_config.get('files', 'strans_lin_file')
        self.df_file = self.project_config.get('files', 'df_file')
        self.plot_df_file = self.project_config.get('files', 'plot_file')
        self.lopt_lev_comments_file = self.project_config.get('files', 'lopt_lev_comments_file')
        self.other_lev_list = self.project_config.get('files', 'other_lev_files').split('\n')

        self.lopt_fixed_levels = self.project_config.get('lopt', 'fixed_levels').split(',')
        self.star_discrim = self.project_config.getfloat('lopt', 'star_discrim')
        self.lopt_plot_width = self.project_config.getfloat('lopt', 'plot_width')
        
        self.strans_wn_discrim = self.project_config.getfloat('strans', 'wn_discrim')
      
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
        
        self.lopt_inp_file = f'lopt/{self.main_element_name}_lopt.inp'
        self.lopt_par_file = f'lopt/{self.main_element_name}_lopt.par'
        self.lopt_fixed_file = f'lopt/{self.main_element_name}_lopt.fixed'
        self.lopt_lev_file = f'lopt/{self.main_element_name}_lopt.lev'
        self.lopt_lin_file = f'lopt/{self.main_element_name}_lopt.lin'
          
        self.load_df() 
        self.load_plot_df()         
        self.strans_levs = list(pd.read_csv(self.strans_lev_file, dtype={'parity':float}).transpose().to_dict().values())                 
        self.display_strans_levs() 
        self.load_lopt_lev_comments()
        self.SetTitle(f"Term Analysis Made Easy (TAME) - {self.project_title}")
    
        
    def load_plot_df(self):
        """Loads plot_df for the matplotlib plot from the user-selected Xgremlin ascii linelist files.
        The plot_df file is created as part of the new project process."""
        # path=os.path.dirname(self.plot_df_file)

        # if not os.path.isfile(self.plot_df_file):  # if no existing DataFrame is present
        #     self.plot_df = pd.concat([pd.read_csv(f, skiprows=4, delim_whitespace=True, names=['wavenumber', f'{f.split("/")[-1].split(".")[0]}']) for f in glob.glob(path + "/*.asc")], ignore_index=True)
        #     self.plot_df.to_pickle(self.plot_df_file) 
        # else:
        #     self.plot_df = pd.read_pickle(self.plot_df_file) 
        self.plot_df = pd.read_pickle(self.plot_df_file) 
        # self.plot_df.sort_values(by=['wavenumber'], ascending=True, inplace=True)
        # self.plot_df.set_index('wavenumber', inplace=True, drop=True)

          
    
    def create_df(self, lines_file):
        """Creates a new pandas DataFrame from a list of lines in 'lines_file' and saves to a pickle file."""
        self.df = pd.read_csv(lines_file, float_precision='high')  # create new dataframe from the input lines file 
        self.df['main_desig'] = np.empty((len(self.df), 0)).tolist()  # append column of empty lists.
        self.df['other_desig'] = np.empty((len(self.df), 0)).tolist()  # append column of empty lists.
        self.df['user_desig'] = ''
        self.df['line_tags'] = [{'ringing': False, 'incorr_assign': False, 'noise': False, 'blend': False, 'user_unc': False, 'multiple_lines':False} for x in range(self.df.shape[0])]
        self.df['comments'] = ''
        self.save_df()
        
    def save_df(self):
        """Saves the main pandas DataFrame self.df to the pickle file."""
        self.df = self.df.sort_values(by=['wavenumber'])
        self.df.to_pickle(self.df_file)    
    
    def main_strans(self, strans_levs):
        """Runs strans for the main element under study"""
        
        self.strans_levs = self.strans_lev_ojlv.GetObjects()
        
        if self.blank_strans_lev in strans_levs:
            wx.MessageBox('STRANS input contains blank levels. \n\nPlease edit or delete these before running STRANS.', 'Blank Levels Found', 
                          wx.OK | wx.ICON_EXCLAMATION)
            self.frame_statusbar.SetStatusText('')
            return False
        
        if not all(len(label) <= 10 for label in [lev['label'] for lev in strans_levs]):
            wx.MessageBox('One or more level labels are too long. The maximum label length is 10 characters. \n\nPlease edit or delete these before running STRANS.', 'Level Label(s) Too Long', 
                          wx.OK | wx.ICON_EXCLAMATION)
            self.frame_statusbar.SetStatusText('')
            return False
        
        self.df['main_desig'] = np.empty((len(self.df), 0)).tolist()  # replaces any values in main_desig column with empty lists
        tag_sep_linelist = self.get_tag_sep_linelist()
        desig_list = self.df[['wavenumber', 'main_desig', 'line_tags']].values.tolist()  # list of all lines in the linelist
          
        matched_lines = self.strans(strans_levs, desig_list, self.main_element_name, tag_sep_linelist)          
        self.df.update(matched_lines)  # update the main df with designations from strans
        self.display_strans_lines()

        return True    
    
    def get_tag_sep_linelist(self):
        """returns a dictionary of lists of lines separated by tag type"""        
        desig_list = {}
        tags = self.df.tags.unique()  # get list of all tags for lines (L, G, P etc.)
        
        for tag in tags:
            desig_list[tag] = self.df[['wavenumber', 'main_desig', 'line_tags']].loc[self.df['tags'] == tag].values.tolist()    
            
        return desig_list
    
    def other_strans(self, other_lev_list):
        """Runs strans for all other elements that could be present in the linelist"""                     
        self.df['other_desig'] = np.empty((len(self.df), 0)).tolist()  # replaces any values in other_desig column with empty lists
        desig_list = self.df[['wavenumber', 'other_desig', 'line_tags']].values.tolist()
        tag_sep_linelist = self.get_tag_sep_linelist()
        
        for other_lev in other_lev_list:
            element_name, level_file = other_lev.split(',')
            strans_levs = list(pd.read_csv(level_file).transpose().to_dict().values())
            matched_lines = self.strans(strans_levs, desig_list, element_name, tag_sep_linelist)

        self.df.update(matched_lines)  # update the main df with designations from strans  
        self.display_strans_lines()

    def strans(self, strans_levs, desig_list, element_name, tag_sep_linelist):
        """Creates list of all possible transitions between levels of opposite parity that obey
        the J selection rule. The list is then compared to all lines in the self.df database and lines with 
        wavenumbers that match within self.strans_wn_discrim are assigned the labels of the even and odd level.
        Inputs:
            strans_levs: list of levels to be used in strans
            desig_list: dict of lists of line dicts separated into tags to be matched by strans
            element_name: name of level's element'
        """     
        self.frame_statusbar.SetStatusText(f'Running Line Matching for {element_name}')
                
        strans_levs_even = [x for x in strans_levs if x['parity']==1]
        strans_levs_odd = [x for x in strans_levs if x['parity']==0]
        
        strans_levs_even = sorted(strans_levs_even, key=lambda x: x['j'])  # sorting by j value.
        strans_levs_odd = sorted(strans_levs_odd, key=lambda x: x['j'])
        
        self.strans_wn_discrim = {'P': 0.1, 'L': 0.02, 'G': 0.02, 'I': 0.02, 'F': 0.02}
        
        for even_lev in strans_levs_even:    
            j_even = even_lev['j']
            
            if j_even == 0.0:  # J selection rule J != 0 to 0
                left_j = bisect.bisect_left(KeyList(strans_levs_odd, key=lambda x: x['j']), j_even + 1)  # returns index of leftmost match
            else:  # the other J selection 
                left_j = bisect.bisect_left(KeyList(strans_levs_odd, key=lambda x: x['j']), j_even - 1)
                
            right_j = bisect.bisect_right(KeyList(strans_levs_odd, key=lambda x: x['j']), j_even + 1)  # returns index of leftmost match
               
            for odd_lev in strans_levs_odd[left_j:right_j]:  # strans_levs_odd now reduced to levels with J values in line with selection rules.             
                label_even = even_lev['label']
                label_odd = odd_lev['label']
                energy_even = even_lev['energy']
                energy_odd = odd_lev['energy']                   
                match_wn = abs(energy_even - energy_odd)
                
                matched_lines = []
                
                for tag in tag_sep_linelist:  # find all matched lines, with wn matching tolerance set by tag type
                    wn_discrim = self.strans_wn_discrim[tag]  # get wn_discrim for the type of tag
                                        
                    left = bisect.bisect_left(KeyList(tag_sep_linelist[tag], key=lambda x: x[0]), match_wn - wn_discrim)  # x[0] is line wavenumber
                    right = bisect.bisect_left(KeyList(tag_sep_linelist[tag], key=lambda x: x[0]), match_wn + wn_discrim)
                
                    matched_lines += tag_sep_linelist[tag][left:right]  # add the list of matched lines to the main matched_lines list

                for matched_line in matched_lines:
                    #matched line is a list of:[line wavenumber, [level assignment dicts], {line tags}]
                    
                    desig_index = next(i for i,v in enumerate(desig_list) if matched_line[0] in v)  # this gives the index of matched_line in the main desig_list
                    
                    if len(matched_lines) > 1:  # multiple lines match this transtion
                        desig_list[desig_index][2]['multiple_lines'] = True
                    
                    if energy_even > energy_odd:  # assign upper and lower levels correctly (LOPT needs them in lower-upper format)
                        upper_lev = label_even
                        lower_lev = label_odd
                    else:
                        upper_lev = label_odd
                        lower_lev = label_even    
                    
                    # because we found the desig_list index, we can modify the desig_list item directly                        
                    desig_list[desig_index][1].append({'upper_level':upper_lev, 'lower_level':lower_lev, 'element_name': element_name})  # this is being added to the lines in desig_list that were matched.
            
        return desig_list
    
             
    def display_strans_levs(self):
        """Writes values from self.strans_levs list to the strans_lev_ojlv ObjectListView"""
        self.strans_lev_ojlv.SetObjects(self.strans_levs)
       
    def display_strans_lines(self):
        """Writes lines with designations from self.df to the strans_lines_ojlv ObjectListView"""   
        strans_lines = list(self.df.loc[self.df.main_desig.str.len() > 0 ].transpose().to_dict().values())  # convert to list of dicts
        
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
    
    def save_project(self):
        """Saves the project."""
        self.save_project_config()
        self.save_main_config()
        self.save_df()
        self.save_strans_levs()
        self.save_lev_comments_df()
    
    def save_lev_comments_df(self):
        """Saves the lopt level comments to the correct pickle file."""
        self.lopt_lev_comments.to_pickle(self.lopt_lev_comments_file) 
   
    def save_strans_levs(self):
        """Writes the levels in TAME to the .lev file. This is needed for user changes that have been made within TAME."""
        with open(self.strans_lev_file, 'w') as lev_file:
            lev_file.write('label,j,energy,parity\n')
            
            for lev in self.strans_levs:
                lev_file.write(f"{lev['label']},{lev['j']},{lev['energy']},{lev['parity']}\n")
               
    def save_project_config(self):
        """Save the project config."""  
        print(self.project_config_file)
        with open(self.project_config_file, 'w') as configfile:
            self.project_config.write(configfile)
                                      
        print('project config saved')
        
    def save_main_config(self):  
        """Saves the main TAME config."""
        with open(self.main_config_file, 'w') as configfile:
            self.main_config.write(configfile)
            
    def write_lopt_inp(self):
        """Writes the LOPT input file. Taking into account user selected tags, uncertainties and multiply identified lines."""
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
                         
                    if user_desig != '':  # there is a user selected level for the line
                        if user_desig['element_name'] == self.main_element_name:  # only if the user selected transition is of the main element
                            upper_level = f'{user_desig["upper_level"]:>12}'
                            lower_level = f'{user_desig["lower_level"]:>12}'
                            
                            if all(value == False for value in tags.values()): # no user defined tags for the line
                                unc = f'{line[5]:.4f}'
                            elif tags['user_unc'] != False:
                                unc = f"{tags['user_unc']:.4f}"
                                tag = '       B'
                            elif tags['multiple_lines'] == True:  # multiple lines could have been a transition, but the user has selected one. 
                                unc = f'{line[5]:.4f}'
                            else:
                                unc = f'{self.lopt_default_unc:.4f}'
                                tag = '       B'
                            
                            lopt_str = f'{snr}{wn} cm-1 {unc}{lower_level}{upper_level}{tag}\n'
                            inp_file.writelines(lopt_str)
                        
                    else:  # no user label for line
                        if len(main_desigs) != 1 or len(other_desigs) !=0:  # multiple identifications for line
                            unc = f'{self.lopt_default_unc:.4f}'
                            tag = '       Q'
                        elif all(value == False for value in tags.values()): # no user defined tags for the line
                            unc = f'{line[5]:.4f}'
                        elif tags['user_unc'] != False:
                            unc = f"{tags['user_unc']:.4f}"
                            tag = '       B'
                        else:
                            unc = f'{self.lopt_default_unc:.4f}'
                            tag = '       B'
                            
                        for desig in main_desigs:
                            upper_level = f'{desig["upper_level"]:>12}'
                            lower_level = f'{desig["lower_level"]:>12}'
                        
                            lopt_str = f'{snr}{wn} cm-1 {unc}{lower_level}{upper_level}{tag}\n'
                            inp_file.writelines(lopt_str)
                            
            # for level in self.lopt_fixed_levels[1:]:  # skip ground
            #     strans_lev = next((item for item in self.strans_levs if item['label']==level))
            #     lev_energy = strans_lev['energy'] 
                
            #     lopt_str = f'     9999{lev_energy:>15} cm-1 0.0000{self.lopt_fixed_levels[0]:>12}{level:>12}\n'
            #     inp_file.writelines(lopt_str)
                
        return True
               
    def write_lopt_par(self):
        """Gets text from the LOPT .par template file and writes .par file for the project."""
        with open('lopt/lopt_template.par', 'r') as temp_par_file:
            par_lines = temp_par_file.readlines()
            par_lines[0] = f'{self.main_element_name}_lopt.inp{par_lines[0]}'
            par_lines[1] = f'{self.main_element_name}_lopt.fixed{par_lines[1]}'
            par_lines[2] = f'{self.main_element_name}_lopt.lev{par_lines[2]}'
            par_lines[3] = f'{self.main_element_name}_lopt.lin{par_lines[3]}'
            
        with open(self.lopt_par_file, 'w') as par_file:
            par_file.writelines(par_lines)
                
    def write_lopt_fixed(self):
        """Writes the fixed levels for LOPT. If the ground level is selected, then unc = 0, otherwise = 2.0."""
        with open(self.lopt_fixed_file, 'w') as fixed_file:
            fixed_strings = []
            
            for level in self.lopt_fixed_levels:
                strans_lev = next((item for item in self.strans_levs if item['label']==level))
                lev_energy = strans_lev['energy'] 
                if lev_energy == 0.0:
                    lev_unc = f'{0.0:.4f}'
                else:
                    lev_unc = f'{2.0:.4f}'
                             
                fixed_strings.append(f'{level:>9}{lev_energy:>13.4f}{lev_unc:>13}\n')
            
            fixed_file.writelines(fixed_strings)    
            
    def get_lopt_output(self):
        """Gets output from the LOPT output files. Puts these into dataframes. Duplicates lines so that the line appears
        in the output GroupListView twice - once for each level in the transition."""
        self.lopt_lev_pos = self.lopt_lev_ojlv.GetTopItem()
        self.lopt_lev_select_row = self.lopt_lev_ojlv.GetFocusedRow()          
        self.lopt_lev_groups_expanded = []
        
        if self.lopt_lev_ojlv.groups:
            for i, group in enumerate(self.lopt_lev_ojlv.groups):
                self.lopt_lev_groups_expanded.append((i, group.isExpanded))

        self.lopt_levs = pd.read_csv(self.lopt_lev_file, delimiter='\t')
        lopt_lines_df = pd.read_csv(self.lopt_lin_file, delimiter='\t')
        merged_lines = pd.merge_asof(lopt_lines_df[['W_obs', 'S', 'Wn_c', 'E1', 'E2', 'L1', 'L2', 'F', 'uncW_o']].sort_values('W_obs'), 
                                     self.df[['wavenumber', 'peak', 'eq width', 'tags']].sort_values('wavenumber'), 
                                     left_on='W_obs', 
                                     right_on='wavenumber',
                                     tolerance=0.005,
                                     direction='nearest') # match lopt lines to main df file based on nearest wavenumber
        
        merged_lines['dWO-C'] = merged_lines['W_obs'] - merged_lines['Wn_c']
        merged_lines['star'] = np.where((merged_lines['dWO-C'].abs() > (merged_lines['uncW_o']*self.star_discrim)), True, False)
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
        
        duplicated_lines.update(t1)
        duplicated_lines.update(t2)
             
        duplicated_lines = list(duplicated_lines.transpose().to_dict().values())
        self.lopt_lev_ojlv.SetObjects(duplicated_lines)         
        self.lopt_output_lines = duplicated_lines        
        # self.load_lopt_lev_comments()    
        self.lopt_lev_ojlv.EnsureVisible(min(self.lopt_lev_pos + self.lopt_lev_ojlv.GetCountPerPage() - 1, self.lopt_lev_ojlv.GetItemCount() - 1))  # ensures that the levels stay in the same scrolled position after lopt has run.
           
        for group in self.lopt_lev_groups_expanded:
            if not group[1]:
                self.lopt_lev_ojlv.Collapse(self.lopt_lev_ojlv.groups[group[0]])
                
        if self.lopt_lev_select_row != -1:  # -1 means no row was focussed.
            self.lopt_lev_ojlv.Select(self.lopt_lev_select_row)  # puts the focus back on the line that was focussed before lopt ran.
        
        
        
    def load_lopt_lev_comments(self):
        """Load the comments for each level. If the pkl file is missing - create it from the strans levels"""        
        if not os.path.isfile(self.lopt_lev_comments_file):  # if no existing DataFrame is present
            strans_desigs = [lev['label'] for lev in self.strans_levs]
            self.lopt_lev_comments = pd.DataFrame(strans_desigs, columns=['Designation'])
            self.lopt_lev_comments['Comments'] = ''
            self.lopt_lev_comments.to_pickle(self.lopt_lev_comments_file) 
        else:
            self.lopt_lev_comments = pd.read_pickle(self.lopt_lev_comments_file)
            # old_lopt_levs = self.lopt_lev_comments['Designation'].values.tolist()        
            # current_lopt_levs= self.lopt_levs['Designation'].values.tolist()            
            # new_levs = [{'Designation':x, 'Comments':''} for x in current_lopt_levs if x not in old_lopt_levs]

            # for lev in new_levs:
            #     self.lopt_lev_comments = self.lopt_lev_comments.append(lev, ignore_index=True)
      
    def display_lopt_line(self, line):  
        """Sets values for the various controls and the line plot in the LOPT line panel.""" 
        
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
        self.incorr_assign_chkbox.SetValue(line_tags['incorr_assign'])
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
        
        self.lopt_level_panel.Hide()
        self.lopt_line_panel.Show()
        self.sizer_8.Layout()
        
     
    def plot_lopt_line(self, wavenumber, peak, width):   
        """Plots the spectra from self.plot_df to the matplotlib plot. There are scaling factors that can be changed to show
        more/less of the plot area."""

        if not self.plot_df.empty:  # plot_df will be empty if the user did not select any plot files at project creation
            plot_width = self.lopt_plot_width / 2  # as this is total width
            # plot_y_scale = 1.1
            # plot_x_scale = 5.
            self.matplotlib_canvas.clear()  
            ax = self.matplotlib_canvas.gca()   
            
            low = abs(self.plot_df['wavenumber'] - (wavenumber-plot_width)).idxmin() + 1
            high = abs(self.plot_df['wavenumber'] - (wavenumber+plot_width)).idxmin() + 1
            df_part = self.plot_df.iloc[low:high]
               
            spectras = df_part.columns[df_part.notna().any()].tolist()[1:]  # gives columns that do not contain only NaN.
    
            for spectra in spectras:
                df_part.plot(kind='line', x='wavenumber', y=spectra, ax=ax)
    
            self.matplotlib_canvas.axes.set_xlabel('Wavenumber (cm-1)')
            self.matplotlib_canvas.axes.set_ylabel('SNR')
            ax.ticklabel_format(useOffset=False)
            # ax.set_ylim([df_part.min().min()*plot_y_scale, peak*plot_y_scale])
            #ax.set_xlim([wavenumber - (width/1000)*plot_x_scale, wavenumber + (width/1000)*plot_x_scale])
            # ax.axvline(x=wavenumber)  # due to calibration this is not in the right place and looks odd
    
            self.matplotlib_canvas.draw()     
    
    def display_lopt_lev(self, level):
        """Display info about the LOPT level. Also the user comments."""
        
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
        """Updates a cell of the main self.df dataframe, specified by wavenumber and column, with value."""   
        selected_line_index = self.df.loc[self.df['wavenumber'] == wavenumber].index.values[0]
        self.df.at[selected_line_index, column] = value  # just updates a single value
        print(self.df.iloc[selected_line_index])
        
    def get_df_cell(self, wavenumber, column):
        """Gets the value of a cell in the main self.df dataframe for a given wavenumber and column."""
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
        
    def set_fixed_levels(self):
        """User dialog window for setting the fixed levels for LOPT."""
        fixed_level_dialog = fixedLevels(self)
        
        if fixed_level_dialog.ShowModal() == wx.ID_CANCEL:
            return False
                
        self.lopt_fixed_levels = fixed_level_dialog.fixed_levels
        level_config_string = ','.join(self.lopt_fixed_levels)
        self.project_config.set('lopt', 'fixed_levels', level_config_string)
        
        return True
        
       
    def export_linelist(self, full_list):
        """Exports the STRANS output linelist with matched transitions. full_list determines whether only lines with
        transitions involving the main element are outputted, or the entire linelist is regardless of whether a 
        transition has been found by STRANS or not."""
        with wx.FileDialog(self, "Export Matched Linelist", wildcard="Linelist Files (*.lin)|*.lin",
                       style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            filename = fileDialog.GetPath()

            with open(filename, 'w') as file:
                
                if full_list:  # user has selected that all lines should be in the outputted linelist
                    lines = self.df.values.tolist()
                    linelist_type = 'Complete'
                else:
                    lines = self.df.loc[self.df.main_desig.str.len() > 0 ].values.tolist()
                    linelist_type = 'Matched'
                    
                file.writelines('wavenumber,snr,fwhm,eq_width,fit,unc,main_element,other_elements\n')  # header
                
                for line in lines:
                    main_level_string = ''
                    other_level_string = ''
                    
                    for lev_pair in line[6]:
                        if not main_level_string:
                            sep = ''
                        else:
                            sep = '\t'                                
                        main_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['upper_level'] + ' - ' + lev_pair['lower_level']
                        
                    for lev_pair in line[7]:
                        if not other_level_string:
                            sep = ''
                        else:
                            sep = '\t'                                
                        other_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['upper_level'] + ' - ' + lev_pair['lower_level']
                    
                    line[0] = f'{line[0]:.4f}'
                    line[1] = f'{line[1]:.0f}' 
                    line[2] = f'{line[2]:.0f}' 
                    line[3] = f'{line[3]:.0f}' 
                    line[5] = f'{line[5]:.4f}' 
                    line[6] = main_level_string
                    line[7] = other_level_string
                    
                    file.writelines(','.join(line[:8]) + '\n') 
                    
        self.frame_statusbar.SetStatusText(f'{linelist_type} linelist exported to {filename}')   
      
        
      
    def levhams_match_tol(self, line_list, tol):
        """Generator function to split predicted lines into groups with each element separated from its neighbour 
        by a maximium given tolerance.
        """
        matches = []
        last = line_list[0]['pred_energy']
        
        for element in line_list:
            if element['pred_energy'] - last > tol:
                yield matches
                matches = []
                
            matches.append(element)
            last = element['pred_energy']
        yield matches    
        
    
    def display_levhams_levs(self, levels):
        """Add the list of grouped level lines from levhams to the listctrl if there are >= the number of selected lines
        in each group. Also highlights the summary row."""
                     
        num_levels = 0
        
        self.levhams_output_ojlv.DeleteAllItems()  # clear the ojlv
        
        for level in levels:            
            if len(level) >= self.levhams_min_matches:
                sum_lines = 0.
                num_levels += 1
                
                for line in level:
                    self.levhams_output_ojlv.AddObject(line)
                    sum_lines += line['pred_energy']
                
                sep = level[-1]['pred_energy'] - level[0]['pred_energy']  # the separation in the level (highest predicted energy minus the lowest) 
                self.levhams_output_ojlv.AddObject({'avg_energy': sum_lines/len(level), 'sep': sep})
                
        self.frame_statusbar.SetStatusText(f'{num_levels} predicted levels found.')
        
        if num_levels == 0:
            wx.MessageBox('No predicted levels found. Please add more levels or change search parameters.', 'No Predicted Levels Found', 
                              wx.OK | wx.ICON_EXCLAMATION)
        
        
      
        
### Event-driven functions ###  

    def on_lopt_lev_sort(self, event):
        """Handles clicks on the columns 0 and 10 which otherwise raise exceptions due to TypeErrors."""
        if event.GetSortColumn() == 0:  # group expand/contract column
            event.Veto()
            num_expanded = 0
            
            for group in self.lopt_lev_ojlv.groups:  # works out if more than 10 groups are expanded
                if group.isExpanded == True:
                    num_expanded +=1
                    
            if num_expanded <= 10:
                self.lopt_lev_ojlv.ExpandAll()
            else:
                self.lopt_lev_ojlv.CollapseAll()             
            
        elif event.GetSortColumn() == 10:  # tags column
            event.Veto()
            

    def on_strans_lev_sort(self, event):
        print(event.GetSortColumn())
        # XXX may want tp put a custom sort in here for the strans levels if the user selects parity to sort by

    def on_levhams_right_click(self, event):  
        self.PopupMenu(LevhamsPopupMenu(self))
        

    def on_notebook_page_change(self, event): 
        """Detects a change in notebook page. This is needed to ensure that the LEVHAMS page always uses the most
        up to date level values.
        """
        if self.main_panel.GetSelection() == 2:  # Update LEVHAMS level listctrl with the latest strans values            
            with wx.EventBlocker(self):  # prevents the event from firing            
                self.levhams_level_listctrl.DeleteAllItems()  # otherwise just appends to prev. values
                
                if not self.levhams_selected_levs:  # populate dict if doesn't exist already
                    for level in self.strans_levs:
                        self.levhams_selected_levs[level['label']] = False
                
                for i, level in enumerate(self.strans_levs):                    
                    levhams_level_list = [
                        level['label'],
                        level['j'],
                        level['energy'],
                        level['parity']]          
                    self.levhams_level_listctrl.Append(levhams_level_list)  
                    
                    if level['label'] in self.levhams_selected_levs:  # set checkboxes
                        self.levhams_level_listctrl.CheckItem(i, self.levhams_selected_levs[level['label']])
                    else:
                        self.levhams_selected_levs[level['label']] = False
                        
    def on_levhams_lev_checked(self, event):  
        """Updates levhams_selected_lines based on user checkbox selection"""           
        level_label = self.levhams_level_listctrl.GetItem(event.GetIndex(), 0).GetText()
        self.levhams_selected_levs[level_label] = True
        

    def on_levhams_lev_unchecked(self, event): 
        """Updates levhams_selected_lines based on user checkbox selection"""        
        level_label = self.levhams_level_listctrl.GetItem(event.GetIndex(), 0).GetText()
        self.levhams_selected_levs[level_label] = False
        
        
    def on_levhams(self, event): 
        """Runs the main levhams code when the user selects run."""
        
        if self.main_panel.GetSelection() != 2:
            self.main_panel.ChangeSelection(2)
            self.on_notebook_page_change(None)
            
        self.frame_statusbar.SetStatusText('Predicting Levels ...') 
        
        self.levhams_tol = self.levhams_tol_spinctrl.GetValue() # XXX add this to the main params and load it at startup.
        self.levhams_min_matches = self.levhams_num_spinctrl.GetValue()  # XXX add this to main params and load at startup   
        self.levhams_use_all_lines = self.levhams_all_rbutton.GetValue()  # if the all lines radio button is selected
        self.levhams_wn_max = self.levhams_wn_max_spinctrl.GetValue()
        self.levhams_wn_min = self.levhams_wn_min_spinctrl.GetValue()
                
        selected_levs = [lev for lev in self.strans_levs if self.levhams_selected_levs[lev['label']]]
        
        if selected_levs:  
            
            pred_lines = []
            all_lines = list(self.df[['wavenumber','peak','eq width','unc', 'main_desig']].transpose().to_dict().values()) # xxx here is where you would run a query for all lines without a desig if you wanted to implement that
            
            if not self.levhams_use_all_lines:
                all_lines = [x for x in all_lines if not x['main_desig']]
            
            for level in selected_levs:  
                for line in all_lines:
                    pred_lines.append({'pred_energy':level['energy'] - line['wavenumber'], 
                                       'level':level['label'], 
                                       'j':level['j'],
                                       'energy':level['energy'],
                                       'wavenumber':line['wavenumber'],
                                       'peak':line['peak'],
                                       'eq width':line['eq width'],
                                       'unc':line['unc']})
                    pred_lines.append({'pred_energy':level['energy'] + line['wavenumber'], 
                                       'level':level['label'], 
                                       'j':level['j'],
                                       'energy':level['energy'],
                                       'wavenumber':line['wavenumber'],
                                       'peak':line['peak'],
                                       'eq width':line['eq width'],
                                       'unc':line['unc']})
            
            pred_lines = sorted(pred_lines, key=lambda k: k['pred_energy'])  # sort by wavenumber            
            pred_lines = [x for x in pred_lines if x['pred_energy'] >= self.levhams_wn_min and x['pred_energy'] <= self.levhams_wn_max]
            
            if pred_lines:  # if there are any predicted lines
                pred_levels = list(self.levhams_match_tol(pred_lines, self.levhams_tol))  # separates predicted lines into groups
                self.display_levhams_levs(pred_levels)
            else:
                wx.MessageBox('No predicted levels found. Please add more levels or change search parameters.', 'No Predicted Levels Found', 
                              wx.OK | wx.ICON_EXCLAMATION)
            
        else:
            wx.MessageBox('No levels selected for Level Prediction. Please tick the levels you wish to use first.', 'No Levels Selected', 
                              wx.OK | wx.ICON_EXCLAMATION)
            
        

    def on_export_lopt_levs(self, event):    
        """Exports a sorted and formatted list of all LOPT levels and their lines."""
        with wx.FileDialog(self, "Export LOPT Sorted Levels", wildcard="Linelist Files (*.llf)|*.llf",
                       style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:
           
            # if 'self.lopt_levs' in globals():  # LOPT has been run            
            try:
                if fileDialog.ShowModal() == wx.ID_CANCEL:
                    return     # the user changed their mind        
                filename = fileDialog.GetPath()
                
                if filename.split('.')[-1] != 'llf':   #if user types new filename in dialog
                    filename += '.llf'                
                lopt_levels = list(self.lopt_levs.transpose().to_dict().values())
    
                with open(filename, 'w') as file:            
                    for level in lopt_levels: 
                        file.write('--------------------------------------------\n')
                        file.write('Config      Energy (cm-1)   D1        D2        D3        No. of Lines  Comments\n')
                        file.write(f"{level['Designation']:<12}{level['Energy']:<16}{level['D1']:<10}{level['D2']:<10}{level['D3']:<10}{level['N_lines']:<14}{level['Comments']}\n\n")
                        file.write('   Fit  log(EW) SNR     Wn_Obs (cm-1)  Unc_Wn_obs  del_Wn(Obs-C)  Level       Tags\n')
                        
                        lopt_lines = [x for x in self.lopt_output_lines if x['other_level'] == level['Designation']]
                        
                        for line in lopt_lines:
                            if line['star']:
                                star = '*'
                            else:
                                star= ' '
                            
                            if line['L1'] == line['other_level']:
                                trans_lev = line['L2']
                            else:
                                trans_lev = line['L1']
                                
                            if line['F'] == 'Q':
                                flag = 'M'
                            else:
                                flag = ' '
                            
                            file.write(f"{star:<3}{line['tags']:<5}{line['log_ew']:<8.2f}{line['peak']:<8.0f}{line['wavenumber']:<15.4f}{line['uncW_o']:<12.4f}{line['dWO-C']:>7.4f}        {trans_lev:<12}{flag}\n")
                        
                        file.write('--------------------------------------------\n')
                        
            except AttributeError:
                wx.MessageBox('No levels found from LOPT output. Please run LOPT first', 'No Levels Found', 
                              wx.OK | wx.ICON_EXCLAMATION)
        
    def on_preferences(self, event):
        """Display TAME preferences dialog and set global variables accordingly."""
        self.prop_dialog = preferenceDialog(self)
        
        if self.prop_dialog.ShowModal() == wx.ID_OK:
            self.strans_wn_discrim = float(self.prop_dialog.strans_wn_discrim.GetValue())
            self.star_discrim = float(self.prop_dialog.lopt_delwn_discrim.GetValue())
            self.lopt_plot_width = float(self.prop_dialog.lopt_plot_width.GetValue()) 
            
            self.project_config.set('lopt', 'star_discrim', str(self.star_discrim))
            self.project_config.set('lopt', 'plot_width', str(self.lopt_plot_width))
            self.project_config.set('lopt', 'wn_discrim', str(self.strans_wn_discrim))
                       
    def on_lopt_lev_comments(self, event):
        """Updates the lopt_lev_comments df with the user entered comments."""
        text = self.lopt_level_comments.GetValue()   
        self.lopt_lev_comments.at[self.selected_lev_index, 'Comments'] = text  # updates level with the comments

    def on_lopt_line_comments(self, event):
        """Updates the main df with the user entered comments."""
        text = self.lopt_line_comments_txtctrl.GetValue()
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']
        self.update_df_cell(selected_wn, 'comments', text)
             
    def on_incorr_assign_tag(self, event): 
        """Updates the main df with the user selected tag."""
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']         
        line_tags = self.get_df_cell(selected_wn, 'line_tags')
        line_tags['incorr_assign'] = self.incorr_assign_chkbox.GetValue()
        self.update_df_cell(selected_wn, 'line_tags', line_tags)

    def on_ringing_tag(self, event):   
        """Updates the main df with the user selected tag."""
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']         
        line_tags = self.get_df_cell(selected_wn, 'line_tags')
        line_tags['ringing'] = self.ringing_chkbox.GetValue()
        self.update_df_cell(selected_wn, 'line_tags', line_tags)

    def on_noise_tag(self, event): 
        """Updates the main df with the user selected tag.""" 
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']         
        line_tags = self.get_df_cell(selected_wn, 'line_tags')
        line_tags['noise'] = self.noise_chkbox.GetValue()
        self.update_df_cell(selected_wn, 'line_tags', line_tags)

    def on_blend_tag(self, event):  
        """Updates the main df with the user selected tag."""
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']         
        line_tags = self.get_df_cell(selected_wn, 'line_tags')
        line_tags['blend'] = self.blend_chkbox.GetValue()
        self.update_df_cell(selected_wn, 'line_tags', line_tags)

    def on_user_unc_tag(self, event):
        """Updates the main df with the user selected tag and its uncertainty."""
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
        """Updates the main df with the user selected line transition."""
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
        """Updates the main df with the user selected tag."""
        selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']       
        self.update_df_cell(selected_wn, 'user_desig', '')

    def on_partial_strans(self, event):
        """Runs STRANS only for the main element. Could save user time for large numbers of impurity files, that do
        not change between STRANS runs."""
        if self.main_strans(self.strans_levs):  # if strans ran succesfully
            self.main_panel.ChangeSelection(0)  # changes the notebook tab to STRANS
            self.frame_statusbar.SetStatusText('Line Matching Complete')
        
    def on_full_strans(self, event):  
        """Runs STRANS for all elements."""
        if self.main_strans(self.strans_levs): # if strans ran succesfully
            self.other_strans(self.other_lev_list)
            self.main_panel.ChangeSelection(0)  # changes the notebook tab to STRANS
            self.frame_statusbar.SetStatusText('Line Matching Complete')
           
    def on_strans_del(self, event):  
        """Delete levels from self.strans_lev_file. This will be a permanent change when saved by user."""
        selected_levs = self.strans_lev_ojlv.GetSelectedObjects()
        
        if selected_levs:  # if selection not empty
            if len(selected_levs) == 1:
                message = 'Are you sure you want to delete this level?'
                title = 'Delete Level?'
            else:
                message = f'Are you sure you want to delete these {len(selected_levs)} levels?'
                title = 'Delete Levels?'
        
            if wx.MessageBox(message, title, wx.YES_NO | wx.NO_DEFAULT | wx.ICON_EXCLAMATION) == wx.YES:                  
                self.strans_levs = [x for x in self.strans_levs if x not in selected_levs]            
                self.display_strans_levs()
            else:
                return
            
    def on_strans_add(self, event):  
        """Add a blank line to the STRANS input levels and display it."""
        self.strans_levs.insert(0, {'label': '', 'j':0.0 , 'energy':0.0 , 'parity':0})  # inserts blank line at head of the table
        self.display_strans_levs()
        
    def on_strans_save(self, event):  
        print("Event handler 'on_strans_save' not implemented!")
        event.Skip()
        
    def on_start_strans_lev_edit(self, event):
        """Triggers when a user double clicks a cell. Keeps previous cell value so dfs can be updated."""
        self.edited_cell_prev_value = event.cellValue    
    
    def on_finish_strans_lev_edit(self, event):  
        """Update the strans_levs list and lopt_lev_comments with user changes. Updates comment df if the user changes
        the designation of a level or creates a new row if a new level is added."""
        if not self.is_float(self.edited_cell_prev_value):  # is the level designation (as all other values are floats)
            new_desig = self.strans_lev_ojlv.GetLastEditedObject()['label'] 
            
            try:  # this will work if the level being edited is in self.lopt_lev_comments
                selected_line_index = self.lopt_lev_comments.loc[self.lopt_lev_comments['Designation'] == self.edited_cell_prev_value].index.values[0]
                self.lopt_lev_comments.at[selected_line_index, 'Designation'] = new_desig # just updates a single value
            except IndexError:  # a newly added level has been edited by the user
                self.lopt_lev_comments = self.lopt_lev_comments.append({'Designation': new_desig, 'Comments': ''}, ignore_index=True)
        
        self.strans_levs = self.strans_lev_ojlv.GetObjects()  # updates the edited cells to strans_levs
            
    def on_export_matched_linelist(self, event):
        """Export matched linelist."""
        self.export_linelist(False)  # export only matched linelist

    def on_export_full_linelist(self, event):
        """Export full linelist."""
        self.export_linelist(True)  # export complete linelist                
                
    def on_Open(self, event):
        """Opens user selected project and loads all variables and dataframes."""
        if self.project_loaded:
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
            except ValueError:  # XXX need to change this back to generic error - or better yet - to just the specific one
                wx.MessageBox('Unsupported or out of date .ini file. Please select a valid TAME project file.', 'Unsupported File', 
                      wx.OK | wx.ICON_EXCLAMATION)    
                
    def on_Save_As(self, event):  
        """Save current project as a separate new project, copying releveant files and updating the project and main
        config files."""
        with wx.FileDialog(self, "Save TAME project as", wildcard="project file (*.ini)|*.ini",
                        style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
                       
            sa_folder, project_title = os.path.split(fileDialog.GetPath())
            sa_folder += '/'
            project_file = project_title.replace(' ', '_')
            
            if project_title.split('.')[-1] != 'ini':   #if user types new filename in dialog
                sa_ini_file = project_file + '.ini'
                
        sa_project_config_file =  sa_folder + sa_ini_file      
        sa_strans_lev_file = sa_folder + project_file + '_input.lev'
        sa_strans_lin_file = sa_folder + project_file + '_input.lin'
        sa_main_df_file = sa_folder + project_file + '.pkl'
        sa_plot_df_file = sa_folder + project_file + '_plot.pkl'
        sa_lev_comments_file = sa_folder + project_file + '_lev_comments.pkl'
        
        copy(self.project_config_file, sa_project_config_file)
        copy(self.strans_lev_file, sa_strans_lev_file)
        copy(self.strans_lin_file, sa_strans_lin_file)
        copy(self.df_file, sa_main_df_file)
        copy(self.plot_df_file, sa_plot_df_file)
        copy(self.lopt_lev_comments_file, sa_lev_comments_file)
        
        self.project_config.set('files', 'strans_lev_file', sa_strans_lev_file)
        self.project_config.set('files', 'strans_lin_file',sa_strans_lin_file)
        self.project_config.set('files', 'df_file',sa_main_df_file)
        self.project_config.set('files', 'plot_file',sa_plot_df_file)
        self.project_config.set('files', 'lopt_lev_comments_file',sa_lev_comments_file)
        self.project_config.set('tame', 'project_title', project_title)
        
        self.main_config.set('project', 'project_config', sa_project_config_file)
        
        self.project_config_file = sa_project_config_file
        self.save_project()
        self.load_project()
                                  
    def on_edit_fixed_levels(self, event):
        """Make changes to user defined fixed levels."""
        self.set_fixed_levels()          
    
    def on_lopt(self, event):
        """Create/update all neccessary files for LOPT input and then call LOPT"""        
        if self.lopt_fixed_levels == [''] or self.lopt_fixed_levels == []:  # no fixed levels have been selected
            if not self.set_fixed_levels():  # if cancel button pressed on fixed level dialog
                return
        
        self.frame_statusbar.SetStatusText('Writing LOPT input files...')
        
        if self.write_lopt_inp():  # if there are lines in the strans linelist
            self.write_lopt_par()
            try:
                self.write_lopt_fixed()   
            except:
                self.set_fixed_levels()
                self.write_lopt_fixed()
                
            self.frame_statusbar.SetStatusText('Running LOPT...')   
            
            try:
                p = subprocess.run(['java', '-jar', 'Lopt.jar', self.lopt_par_file.split('/')[-1]], cwd='lopt/', capture_output=True, text=True).stdout.split('\n')  # run LOPT and get output as a list of lines
                rss = [x for x in p if 'RSS' in x]  # gives the RSS\degrees_of_freedom line
                tot_time = [x for x in p if 'Total time' in x]  # gives the total time line
                self.frame_statusbar.SetStatusText(f'LOPT ran successfully:  {rss[0]}. {tot_time[0]}.')  
                
                self.get_lopt_output()
                self.main_panel.ChangeSelection(1)  # changes the notebook tab to LOPT
            except FileNotFoundError as fnf:                
                if 'java' in str(fnf):
                    self.frame_statusbar.SetStatusText('LOPT error') 
                    wx.MessageBox('Java Runtime Environment (JRE) is not installed on this machine. \n\nPlease install and restart Tame.', 'Missing Java runtime', 
                          wx.OK | wx.ICON_EXCLAMATION)
                    
            # except:  # XXX this needs work - put it in a scrolled dialog
            #     wx.MessageBox('\n'.join(p), 'LOPT Issue', 
            #               wx.OK | wx.ICON_EXCLAMATION)
        else:
            self.frame_statusbar.SetStatusText('Run Line Matching first') 
                        
    def on_Save(self, event):  
        """Saves all user changes for the project."""
        self.save_project()
        self.frame_statusbar.SetStatusText('Project saved')                   
        
    def on_Exit(self, event):  
        """Exits TAME after ensuring changes have been saved if user wanted them to be."""    
        if self.project_loaded:
            save_dlg = wx.MessageBox(f'Do you want to save changes to {self.project_title}', 'Save Changes?', 
                                     wx.YES_NO | wx.CANCEL | wx.CANCEL_DEFAULT | wx.ICON_INFORMATION)
            if save_dlg == wx.YES:
                self.save_project()
                self.Destroy()
            elif save_dlg == wx.NO:
                self.Destroy() 
            else:
                return   
        self.Destroy()
        
    def on_click_lopt_levs(self, event):  
        """Event for user selecting a line,level or blank line in the LOPT GroupListView."""
        try:
            selected_wn = self.lopt_lev_ojlv.GetSelectedObject()['wavenumber']
            selected_line = self.df.loc[self.df['wavenumber'] == selected_wn]              
            self.display_lopt_line(selected_line)

        except TypeError: # group header or blank row selected
            next_row = self.lopt_lev_ojlv.GetObjectAt(self.lopt_lev_ojlv.GetFocusedRow()+1)
            
            if not next_row == None:  # group header selected and not the blank line above it         
                next_row_lev = next_row['main_level']
                selected_lev = self.lopt_levs.loc[self.lopt_levs['Energy'] == next_row_lev]                
                self.display_lopt_lev(selected_lev)    
            else:
                self.lopt_level_panel.Hide()
                self.lopt_line_panel.Hide()
                self.sizer_8.Layout()
                
        except IndexError:  # virtual line selected
            self.lopt_level_panel.Hide()
            self.lopt_line_panel.Hide()
            self.sizer_8.Layout()
                  
    def on_col_clicked(self, event):  
        #self.lopt_lev_ojlv.SortBy(1)
        # self.lopt_lev_ojlv.RebuildGroups()
        event.Skip()
        
    def on_strans_lev_search(self, event):
        """User has entered text to search the STRANS levs for."""
        self.search_listview(event, self.strans_lev_ojlv)
        
    def on_strans_line_search(self, event): 
        """User has entered text to search the STRANS lines for."""
        search_bar = event.GetEventObject()
        try:
            self.search_listview(event, self.strans_lines_ojlv)
            search_bar.SetBackgroundColour(wx.WHITE)
        except:
            search_bar.Clear()
            search_bar.SetBackgroundColour(wx.RED)
            
    def on_New(self, event):  
        """Opens new project dialog and sets all variables and creates all needed files for a new project."""
        self.new_proj = newProject(self)
        self.new_proj.ShowModal()  
        
        #try:
        self.template_config_file = 'config/template.ini'
        self.new_config = configparser.ConfigParser()
        self.new_config.read(self.template_config_file)
        
        self.new_config.set('files', 'strans_lev_file', self.new_proj.main_element_lev_file)
        self.new_config.set('files', 'strans_lin_file', self.new_proj.linelist_file)
        self.new_config.set('files', 'df_file', self.new_proj.df_file)
        self.plot_df_file = self.new_proj.project_file_name + '_plot.pkl'
        self.new_config.set('files', 'plot_file', self.plot_df_file)              
        self.new_config.set('files', 'lopt_lev_comments_file', self.new_proj.project_file_name + '_lev_comments.pkl')
       
        other_lev_files = ''
        for file in self.new_proj.other_element_lev_files:
            other_lev_files += f'{file["shortname"]},{file["filename"]}\n'
        self.new_config.set('files', 'other_lev_files', other_lev_files)            
        self.new_config.set('tame', 'project_title', self.new_proj.project_name)
        self.new_config.set('tame', 'main_element_name', self.new_proj.main_element_name)            
        self.new_config.set('lopt', 'fixed_levels', '')  # no fixed levels for new project            
        self.new_config_file = self.new_proj.project_file_name + '.ini'
        self.frame_statusbar.SetStatusText('Creating new project')
        
        with open(self.new_config_file, 'w') as configfile:
            self.new_config.write(configfile)

        self.frame_statusbar.SetStatusText('Writing project configuration')  
        
        if self.new_proj.plot_files:  # if plot files have been selected by the user
            self.plot_df = pd.concat([pd.read_csv(f, skiprows=4, delim_whitespace=True, names=['wavenumber', f'{f.split("/")[-1].split(".")[0]}']) for f in self.new_proj.plot_files], ignore_index=True)
            self.frame_statusbar.SetStatusText('Creating line plot database')

        else:
            self.plot_df = pd.DataFrame(columns=['wavenumber', 'file'])  # blank dataframe
            
        self.plot_df.to_pickle(self.plot_df_file) 
        self.frame_statusbar.SetStatusText('Saving line plot database to file')
        
        self.project_config_file = self.new_config_file
        
        self.main_config.set('project', 'project_config', self.project_config_file)
        self.save_main_config()
        self.frame_statusbar.SetStatusText('Saving project configuration')
        
        self.load_project()          
        self.frame_statusbar.SetStatusText('New project created successfully')
        self.strans_lines_ojlv.DeleteAllItems()
        self.lopt_lev_ojlv.DeleteAllItems()
        self.main_panel.ChangeSelection(0)
            
        # except AttributeError:  # if the user cancelled the new project process.
        #     print('uh oh')
        #     event.Skip()
  
    def on_about(self, event):  
        """Displays the about dialog."""
        description = """
        Term Analysis Made Easy (TAME) is designed to make the process of term analysis more user-friendly and promote the standardised storage and transfer of project files.
        
        TAME uses a new implementation of the STRANS code (and has retained the name for historical reasons) to find matching lines in a linelist from a list of user-supplied energy levels. These matched lines are passed to the least-squares optimisation program (LOPT) to calculate optimised energy levels.
        
        To download the latest TAME version or submit bug reports and improvement suggestions, please visit the project homepage. 
        """     
        aboutInfo = wx.adv.AboutDialogInfo()
        aboutInfo.SetName("Term Analysis Made Easy (TAME)")
        aboutInfo.SetVersion(TAME_VERSION_STRING)
        aboutInfo.SetDescription(description)
        aboutInfo.SetCopyright("(C) 2021")
        aboutInfo.SetWebSite("https://github.com/Christian-Clear/TASS")
        aboutInfo.AddDeveloper("TAME: Christian Clear \nLOPT: Alexander Kramida")
        wx.adv.AboutBox(aboutInfo)
 
    def on_lost_lines(self, event):
        """Shows the lostLinesDialog and removes the user designation from any selected lines."""
        lost_lines_dialog = lostLines(self)
        
        if lost_lines_dialog.ShowModal() == wx.ID_OK:
            checked_lines = lost_lines_dialog.checked_lines
            
            for line_wn in checked_lines:
                self.update_df_cell(float(line_wn), 'user_desig', '')
    
  
class fixedLevels(fixedLevelsDialog):
    """Dialog class for the fixed level selection."""
    def __init__(self, *args, **kwds):
        """Initilise and set the checkboxes of the current fixed levels."""
        fixedLevelsDialog.__init__(self, *args, **kwds)  
        strans_levs = self.GetParent().strans_levs
        self.fixed_levels = [x for x in self.GetParent().lopt_fixed_levels if x != '']
        self.fixed_levels = [x['label'] for x in strans_levs if x['label'] in self.fixed_levels]  # needed in case a user changes the label of a fixed level.
        self.fixed_level_lc.EnableCheckBoxes(True)
        self.fixed_level_lc.DeleteAllItems()
                
        for level in strans_levs:
            self.fixed_level_lc.Append([level['label'], f"{level['energy']:.4f}"])

        ticked_lev_idx = [strans_levs.index(x) for x in strans_levs if x['label'] in self.fixed_levels]
        with wx.EventBlocker(self):  # stops item checked event from firing
            for idx in ticked_lev_idx:
                self.fixed_level_lc.CheckItem(idx)    
        
        
            
    def on_fixed_lev_checked(self, event):
        """User has checked the checkbox of a level. Update fixed_levels"""
        level = self.fixed_level_lc.GetItem(event.GetIndex(), 0).GetText()
        self.fixed_levels.append(level)
               
    def on_fixed_lev_unchecked(self, event):
        """User has unhecked the checkbox of a level. Update fixed_levels"""
        level = self.fixed_level_lc.GetItem(event.GetIndex(), 0).GetText()
        self.fixed_levels = [x for x in self.fixed_levels if x != level]
    
    def on_fixed_lev_ok(self, event):
        """User has clicked OK button. Check that at least one level has been selected to be fixed."""
                
        if self.fixed_levels == []:  # no fixed level(s) selected
            wx.MessageBox('At least one level must be fixed.', 'Fixed Level Error', 
                          wx.OK | wx.ICON_EXCLAMATION)
        else:
            self.EndModal(wx.ID_OK)
            
            
class lostLines(lostLinesDialog):
    """Dialog class for the lost lines selection. Restores lines that have been removed from the LOPT fit by the user, 
    that would otherwise not be visible anywhere in TAME."""
    
    def __init__(self, *args, **kwds):
        """Populate the list ctrl with any lines that have a user specified designation that is not in the current
        strans input levels file."""
        lostLinesDialog.__init__(self, *args, **kwds)  
        
        self.lost_lines_lc.EnableCheckBoxes(True)
        self.lost_lines_lc.DeleteAllItems()
        self.checked_lines = []
        
        strans_levs = [x['label'] for x in self.GetParent().strans_levs]
        lines = self.GetParent().df.loc[self.GetParent().df.user_desig.str.len() > 0 ].values.tolist()
        
        
        for line in lines:
            wn = line[0]
            user_desig = dict(line[8])
            
            if user_desig['element_name'] != self.GetParent().main_element_name:  
                if user_desig['upper_level'] not in strans_levs or user_desig['lower_level'] not in strans_levs:
                    user_desig_str = f"{user_desig['element_name']}: {user_desig['upper_level']} - {user_desig['lower_level']}" 
                    self.lost_lines_lc.Append([wn, user_desig_str])
        
    def on_lost_lines_checked(self, event):
        """User has checked the checkbox of a line. Update checked_lines"""
        line = self.lost_lines_lc.GetItem(event.GetIndex(), 0).GetText()
        self.checked_lines.append(line)
               
    def on_lost_lines_unchecked(self, event):
        """User has unhecked the checkbox of a line. Update checked_lines"""
        line = self.lost_lines_lc.GetItem(event.GetIndex(), 0).GetText()
        self.checked_lines = [x for x in self.checked_lines if x != line]      
    
    def on_lost_lines_ok(self, event):           
        self.EndModal(wx.ID_OK)
        
        
             
class newProject(newProjectDialog):
    def __init__(self, *args, **kwds):
        newProjectDialog.__init__(self, *args, **kwds)
        self.button_9.Disable()  # greys out back button on first window
        self.other_element_lev_files = []

        self.other_element_levs_ojlv.SetColumns([
            ColumnDefn("Filename", "left", 500, 'filename'),          
            ColumnDefn("Element Short Name", "left", 50, 'shortname', isSpaceFilling=True, isEditable=True)])
        
        self.other_element_levs_ojlv.SetEmptyListMsg('')
        self.other_element_levs_ojlv.cellEditMode = self.other_element_levs_ojlv.CELLEDIT_SINGLECLICK
        
        self.plot_files = None  # in case the user does not specify plot files
        
    def on_project_file_btn(self, event):
        
        with wx.DirDialog (None, "Select Project Folder", "", wx.DD_DEFAULT_STYLE) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind

                
            self.project_path = fileDialog.GetPath() + '/'
            self.project_path_tc.SetValue(str(self.project_path))   
        
    def on_main_lev_file_btn(self, event):
        with wx.FileDialog(self, "Select Main Element Level File", wildcard="level files (*.lev)|*.lev",
                       defaultDir=self.project_path,
                       style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            self.main_element_lev_file = fileDialog.GetPath()
            self.main_lev_file_tc.SetValue(str(self.main_element_lev_file))

    def on_linelist_file_btn(self, event):  # wxGlade: NewProjectFrame.<event_handler>
        with wx.FileDialog(self, "Select Linelist File", wildcard="linelist files (*.lin)|*.lin",
                       defaultDir=self.project_path,
                       style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            self.linelist_file = fileDialog.GetPath()
            self.linelist_file_tc.SetValue(str(self.linelist_file))

    def on_plot_files_btn(self, event):  # wxGlade: NewProjectFrame.<event_handler>
        with wx.FileDialog(self, "Select linelist file", wildcard="plot files (*.asc)|*.asc",
                       defaultDir=self.project_path,
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
        self.main_element_name = self.element_name_tc.GetValue()
        self.project_file_name = self.project_path + self.project_name.replace(' ', '_')
        self.df_file = self.project_file_name + '.pkl'
        
        if not self.project_name or not self.main_element_name or not self.linelist_file:# or not self.plot_files:
            wx.MessageBox('Please fill in all fields before continuing', 'Missing Project Parameters', 
                      wx.OK | wx.ICON_EXCLAMATION)
        else:
            self.new_project_nb.ChangeSelection(1)

    def on_other_lev_files_btn(self, event):
        with wx.FileDialog(self, "Select linelist file", wildcard="level files (*.lev)|*.lev",
                       defaultDir=self.project_path,
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
                

class LevhamsPopupMenu(wx.Menu):

    def __init__(self, parent):
        super(LevhamsPopupMenu, self).__init__()
        self.parent = parent

        ami = wx.MenuItem(self, wx.NewId(), 'Add level to line(s) matching')
        self.Append(ami)
        self.Bind(wx.EVT_MENU, self.OnAdd, ami)

    def OnAdd(self, event):
        for level in self.parent.levhams_output_ojlv.GetSelectedObjects():
            if 'avg_energy' in level.keys():
                self.parent.strans_levs.insert(0, {'label': '', 'j':0.0 , 'energy':level['avg_energy'] , 'parity':0})  # inserts blank line at head of the table
                self.parent.display_strans_levs()
                self.parent.main_panel.ChangeSelection(0)
            

        
                
class preferenceDialog(propertiesDialog):
    def __init__(self, *args, **kwds):   
        propertiesDialog.__init__(self, *args, **kwds)        
        self.parent = self.GetParent()
        
        self.strans_wn_discrim.SetValue(f'{self.parent.strans_wn_discrim}')
        self.lopt_delwn_discrim.SetValue(f'{self.parent.star_discrim}')
        self.lopt_plot_width.SetValue(f'{self.parent.lopt_plot_width}')
        
    def on_strans_wn_discrim(self, event):  
        if not self.parent.is_float(self.strans_wn_discrim.GetValue()):
            wx.MessageBox('This field must be a number', 'Incorrect Parameter', 
                          wx.OK | wx.ICON_EXCLAMATION)
            self.strans_wn_discrim.SetValue(f'{self.parent.strans_wn_discrim}')

    def on_lopt_delwn_discrim(self, event):  
        if not self.parent.is_float(self.lopt_delwn_discrim.GetValue()):
            wx.MessageBox('This field must be a number', 'Incorrect Parameter', 
                          wx.OK | wx.ICON_EXCLAMATION)
            self.lopt_delwn_discrim.SetValue(f'{self.parent.strans_wn_discrim}')

    def on_lopt_plot_width(self, event):  
        if not self.parent.is_float(self.lopt_plot_width.GetValue()):
            wx.MessageBox('This field must be a number', 'Incorrect Parameter', 
                          wx.OK | wx.ICON_EXCLAMATION)
            self.lopt_plot_width.SetValue(f'{self.parent.strans_wn_discrim}')
         
        
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
