#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import wx
import numpy as np
import pandas as pd
import os
from TAME_GUI import mainWindow
import os.path
import bisect
import configparser


class MyFrame(mainWindow):
    def __init__(self, *args, **kwds):
        mainWindow.__init__(self, *args, **kwds)
             
        self.load_main_config()
        self.load_project()  
        #print(self.df.head())
        

        
         
    def display_strans_levs(self):
        """Writes values from a list to the strans_level_list listctrl"""
        self.strans_level_list.DeleteAllItems()
        for lev in self.strans_levs:
            self.strans_level_list.Append(lev)

            
    def display_strans_lines(self):
        """Writes lines with designations from self.df to the strans_lines_list listctrl"""   
        self.strans_lines_list.DeleteAllItems()
        lines = self.df.loc[self.df.main_desig.str.len() > 0 ].values.tolist()

        for line in lines:
            main_level_string = ''
            other_level_string = ''
            for lev_pair in line[6]:
                if not main_level_string:
                    sep = ''
                else:
                    sep = ',     '
                main_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['even_level'] + ' - ' + lev_pair['odd_level']
                
            for lev_pair in line[7]:
                if not other_level_string:
                    sep = ''
                else:
                    sep = ',     '
                other_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['even_level'] + ' - ' + lev_pair['odd_level']
                
            
            line[0] = f'{line[0]:.4f}'
            line[1] = f'{line[1]:.0f}' 
            line[2] = f'{line[2]:.0f}' 
            line[3] = f'{line[3]:.0f}' 
            line[5] = f'{line[5]:.4f}' 
            line[6] = main_level_string
            line[7] = other_level_string
            self.strans_lines_list.Append(line)
                    
    
    def create_df(self, lines_file):
        """Creates a new pandas DataFrame from a list of lines in 'lines_file' and saves to a pickle file."""
        self.df = pd.read_csv(lines_file)  # create new dataframe from the input lines file 
        self.df['main_desig'] = np.empty((len(self.df), 0)).tolist()  # append column of empty lists.
        self.df['other_desig'] = np.empty((len(self.df), 0)).tolist()  # append column of empty lists.
        self.save_df()
        
    
    def save_df(self):
        """Saves the main pandas DataFrame self.df to the pickle file."""
        self.df.to_pickle(self.df_file)    
    
    
    def main_strans(self, strans_levs):
        """Runs strans for the main element under study"""
                     
        self.df['main_desig'] = np.empty((len(self.df), 0)).tolist()  # replaces any values in main_desig column with empty lists
        desig_list = self.df[['wavenumber', 'main_desig']].values.tolist()
          
        matched_lines = self.strans(strans_levs, desig_list, self.main_element_name)
            
        self.df.update(matched_lines)  
        self.display_strans_lines()
        
        
    def other_strans(self, other_lev_list):
        """Runs strans for all other elements that could be present in the linelist"""
                     
        self.df['other_desig'] = np.empty((len(self.df), 0)).tolist()  # replaces any values in other_desig column with empty lists
        desig_list = self.df[['wavenumber', 'other_desig']].values.tolist()
        
        for other_lev in other_lev_list:
            element_name, level_file = other_lev.split(',')
            strans_levs = np.loadtxt(level_file, dtype=str, delimiter=',', skiprows=1) 
            matched_lines = self.strans(strans_levs, desig_list, element_name)
            
        self.df.update(matched_lines)  
        self.display_strans_lines()
        
        
    def strans(self, strans_levs, desig_list, element_name):
        """Creates list of all possible transitions between levels of opposite parity that obey
        the J selection rule. The list is then compared to all lines in the self.db database and lines with 
        wavenumbers that within self.strans_wn_discrim are assigned the labels of the even and odd level."""
        
        strans_levs_even = [x for x in strans_levs if x[-1]=='1']  # NOTE: may want to move this to __init__
        strans_levs_odd = [x for x in strans_levs if x[-1]=='0']
        
        strans_levs_even = sorted(strans_levs_even, key=lambda x: x[1])
        strans_levs_odd = sorted(strans_levs_odd, key=lambda x: x[1])
        
        
        for even_lev in strans_levs_even:    
            j_even = float(even_lev[1])
            
            if j_even == 0.0:  # J selection rule J != 0 to 0
                left_j = bisect.bisect_left(KeyList(strans_levs_odd, key=lambda x: float(x[1])), j_even + 1)
            else:  # the other J selection 
                left_j = bisect.bisect_left(KeyList(strans_levs_odd, key=lambda x: float(x[1])), j_even - 1)
                
            right_j = bisect.bisect_right(KeyList(strans_levs_odd, key=lambda x: float(x[1])), j_even + 1)
                            
            for odd_lev in strans_levs_odd[left_j:right_j]:                
                label_even = even_lev[0]
                label_odd = odd_lev[0]
                energy_even = float(even_lev[2])
                energy_odd = float(odd_lev[2])                    
                match_wn = abs(energy_even - energy_odd)
                    
                left = bisect.bisect_left(KeyList(desig_list, key=lambda x: x[0]), match_wn - self.strans_wn_discrim)
                right = bisect.bisect_left(KeyList(desig_list, key=lambda x: x[0]), match_wn + self.strans_wn_discrim)
                    
                for matched_line in desig_list[left:right]:
                    matched_line[1].append({'even_level': label_even, 'odd_level':label_odd, 'element_name': element_name})
                    
        return desig_list
        
        
        
    def load_main_config(self):
        """Reads the main TAME config file and sets variables accordingly"""
        self.main_config_file = 'tame.ini'
        self.main_config = configparser.ConfigParser()
        self.main_config.read(self.main_config_file)
        
        self.strans_wn_discrim = self.main_config.getfloat('strans', 'wn_discrim')
        self.project_config_file = self.main_config.get('project', 'project_config')
        
        
    def load_project_config(self):
        """Reads the project config file and sets variables accordingly"""
        self.project_config = configparser.ConfigParser()
        self.project_config.read(self.project_config_file)
        
        self.strans_lev_file = self.project_config.get('files', 'strans_lev_file')
        self.strans_lin_file = self.project_config.get('files', 'strans_lin_file')
        self.df_file = self.project_config.get('files', 'df_file')
        self.other_lev_list = self.project_config.get('files', 'other_lev_files').split('\n')
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
        self.load_df()               
        self.strans_levs = np.loadtxt(self.strans_lev_file, dtype=str, delimiter=',', skiprows=1)        
        self.display_strans_levs() 
        self.SetTitle(f"Term Analysis Made Easy (TAME) - {self.project_title}")
        
        
    def create_project_config(self):
        """Creates config file for a new project"""
        self.project_config = configparser.ConfigParser()
        
    
    def save_project(self):
        """Save changes to project .ini file, .lev file, df file"""
        
        self.save_df()
        
        with open(self.main_config_file, 'w') as configfile:
            self.main_config.write(configfile)


### Event-driven functions ###            
    
    def on_partial_strans(self, event):  # Run >> Strans(partial)
        self.frame_statusbar.SetStatusText('Running Strans ...')
        self.main_strans(self.strans_levs)
        self.frame_statusbar.SetStatusText('Strans Complete')
        
    def on_full_strans(self, event):  
        self.frame_statusbar.SetStatusText('Running Full Strans ...')
        self.main_strans(self.strans_levs)
        self.other_strans(self.other_lev_list)
        self.frame_statusbar.SetStatusText('Strans Complete')

            
    def on_strans_del(self, event):  
        """Delete levels from self.strans_lev_file. This will be a permanent change."""
        # item = -1
        # for i in range(self.strans_level_list.GetSelectedItemCount()):
        #     item = self.strans_level_list.GetNextSelected(item)
        #     print(item)
        if wx.MessageBox(f'Are you sure you want to delete these {self.strans_level_list.GetSelectedItemCount()} levels?', 'Delete Levels?', 
                      wx.YES_NO | wx.NO_DEFAULT | wx.ICON_EXCLAMATION) == wx.YES:
            print('Deleted!')
        else:
            return
        
    def on_Export_Linelist(self, event):  
        with wx.FileDialog(self, "Export Matched Linelist", wildcard="Linelist Files (*.lin)|*.lin",
                       style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            filename = fileDialog.GetPath()
            try:
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
                            main_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['even_level'] + ' - ' + lev_pair['odd_level']
                            
                        for lev_pair in line[7]:
                            if not other_level_string:
                                sep = ''
                            else:
                                sep = ',     '
                            other_level_string += sep + lev_pair['element_name'] + ': ' + lev_pair['even_level'] + ' - ' + lev_pair['odd_level']
                            
                        
                        line[0] = f'{line[0]:.4f}'
                        line[1] = f'{line[1]:.0f}' 
                        line[2] = f'{line[2]:.0f}' 
                        line[3] = f'{line[3]:.0f}' 
                        line[5] = f'{line[5]:.4f}' 
                        line[6] = main_level_string
                        line[7] = other_level_string
                        
                        file.writelines(','.join(line) + '\n')                    
                
            except:
                wx.MessageBox(f'Unsupported .ini file. Please select a TAME project file.', 'Unsupported File', 
                      wx.OK | wx.ICON_EXCLAMATION)
            
        
    def on_Open(self, event):  # File >> Open
        save_dlg = wx.MessageBox(f'Do you want to save changes to {self.project_title}', 'Save Changes?', 
                                 wx.YES_NO | wx.CANCEL | wx.CANCEL_DEFAULT | wx.ICON_INFORMATION)
        if save_dlg == wx.YES:
            self.save_project()
        elif save_dlg == wx.CANCEL:
            return

        
        with wx.FileDialog(self, "Open TAME project file", wildcard="project files (*.ini)|*.ini",
                       style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:

            if fileDialog.ShowModal() == wx.ID_CANCEL:
                return     # the user changed their mind
    
            self.project_config_file = fileDialog.GetPath()
            try:
                self.load_project()
                self.strans_lines_list.DeleteAllItems()
                self.main_config.set('project', 'project_config', self.project_config_file)
                self.frame_statusbar.SetStatusText('Project loaded successfully')
                
            except:
                wx.MessageBox(f'Unsupported .ini file. Please select a TAME project file.', 'Unsupported File', 
                      wx.OK | wx.ICON_EXCLAMATION)
                
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
        self.SetTopWindow(self.frame)
        self.frame.Show()
        return True


if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()