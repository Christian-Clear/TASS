#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
#
# generated by wxGlade 1.0.1 on Wed Mar  3 14:42:30 2021
#

import wx

# begin wxGlade: dependencies
# end wxGlade

# begin wxGlade: extracode
from ObjectListView import ObjectListView, ColumnDefn, GroupListView
# end wxGlade


class mainWindow(wx.Frame):
    def __init__(self, *args, **kwds):
        # begin wxGlade: mainWindow.__init__
        kwds["style"] = kwds.get("style", 0) | wx.DEFAULT_FRAME_STYLE | wx.STAY_ON_TOP
        wx.Frame.__init__(self, *args, **kwds)
        self.SetSize((821, 801))
        self.SetTitle("Term Analysis Made Easy (TAME)")

        # Menu Bar
        self.frame_menubar = wx.MenuBar()
        wxglade_tmp_menu = wx.Menu()
        wxglade_tmp_menu.Append(wx.ID_ANY, "New Project", "")
        wxglade_tmp_menu.AppendSeparator()
        item = wxglade_tmp_menu.Append(wx.ID_ANY, "Open ...", "")
        self.Bind(wx.EVT_MENU, self.on_Open, item)
        wxglade_tmp_menu.AppendSeparator()
        item = wxglade_tmp_menu.Append(wx.ID_ANY, "Save", "")
        self.Bind(wx.EVT_MENU, self.on_Save, item)
        item = wxglade_tmp_menu.Append(wx.ID_ANY, "Save As ...", "")
        self.Bind(wx.EVT_MENU, self.on_Save_As, item)
        wxglade_tmp_menu.AppendSeparator()
        wxglade_tmp_menu_sub = wx.Menu()
        item = wxglade_tmp_menu_sub.Append(wx.ID_ANY, "Matched Linelist (.lin)", "")
        self.Bind(wx.EVT_MENU, self.on_Export_Linelist, item)
        wxglade_tmp_menu.Append(wx.ID_ANY, "Export ...", wxglade_tmp_menu_sub, "")
        wxglade_tmp_menu.AppendSeparator()
        item = wxglade_tmp_menu.Append(wx.ID_ANY, "Exit", "")
        self.Bind(wx.EVT_MENU, self.on_Exit, item)
        self.frame_menubar.Append(wxglade_tmp_menu, "File")
        wxglade_tmp_menu = wx.Menu()
        item = wxglade_tmp_menu.Append(wx.ID_ANY, "Strans (partial)", "")
        self.Bind(wx.EVT_MENU, self.on_partial_strans, item)
        item = wxglade_tmp_menu.Append(wx.ID_ANY, "Strans (full)", "")
        self.Bind(wx.EVT_MENU, self.on_full_strans, item)
        wxglade_tmp_menu.AppendSeparator()
        item = wxglade_tmp_menu.Append(wx.ID_ANY, "LOPT", "")
        self.Bind(wx.EVT_MENU, self.on_lopt, item)
        self.frame_menubar.Append(wxglade_tmp_menu, "Run")
        self.SetMenuBar(self.frame_menubar)
        # Menu Bar end

        self.frame_statusbar = self.CreateStatusBar(1)
        self.frame_statusbar.SetStatusWidths([-1])

        self.panel_1 = wx.Panel(self, wx.ID_ANY)

        sizer_1 = wx.BoxSizer(wx.VERTICAL)

        self.main_panel = wx.Notebook(self.panel_1, wx.ID_ANY)
        sizer_1.Add(self.main_panel, 1, wx.EXPAND, 0)

        self.STRANS = wx.Panel(self.main_panel, wx.ID_ANY)
        self.main_panel.AddPage(self.STRANS, "STRANS")

        strans_main_sizer = wx.BoxSizer(wx.VERTICAL)

        self.window_1 = wx.SplitterWindow(self.STRANS, wx.ID_ANY)
        self.window_1.SetMinimumPaneSize(20)
        strans_main_sizer.Add(self.window_1, 1, wx.EXPAND, 0)

        self.window_1_pane_1 = wx.Panel(self.window_1, wx.ID_ANY, style=wx.BORDER_SUNKEN | wx.TAB_TRAVERSAL)

        sizer_2 = wx.BoxSizer(wx.VERTICAL)

        label_1 = wx.StaticText(self.window_1_pane_1, wx.ID_ANY, "Input Levels")
        sizer_2.Add(label_1, 0, wx.ALL, 2)

        self.strans_lev_ojlv = ObjectListView(self.window_1_pane_1, wx.ID_ANY, style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        sizer_2.Add(self.strans_lev_ojlv, 1, wx.EXPAND, 0)

        sizer_6 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_2.Add(sizer_6, 0, wx.EXPAND, 0)

        self.strans_add_levels_btn = wx.Button(self.window_1_pane_1, wx.ID_ANY, "Add Levels")
        sizer_6.Add(self.strans_add_levels_btn, 1, wx.EXPAND, 0)

        self.strans_del_levels_btn = wx.Button(self.window_1_pane_1, wx.ID_ANY, "Delete Levels")
        sizer_6.Add(self.strans_del_levels_btn, 1, 0, 0)

        self.strans_save_levels_btn = wx.Button(self.window_1_pane_1, wx.ID_ANY, "   Save Level Changes   ")
        sizer_6.Add(self.strans_save_levels_btn, 1, 0, 0)

        self.window_1_pane_2 = wx.Panel(self.window_1, wx.ID_ANY, style=wx.BORDER_SUNKEN | wx.TAB_TRAVERSAL)

        sizer_5 = wx.BoxSizer(wx.VERTICAL)

        label_2 = wx.StaticText(self.window_1_pane_2, wx.ID_ANY, "Identified Lines")
        sizer_5.Add(label_2, 0, wx.ALL, 2)

        self.strans_lines_ojlv = ObjectListView(self.window_1_pane_2, wx.ID_ANY, style=wx.LC_REPORT|wx.SUNKEN_BORDER)
        sizer_5.Add(self.strans_lines_ojlv, 1, wx.EXPAND, 0)

        self.LOPT = wx.Panel(self.main_panel, wx.ID_ANY)
        self.main_panel.AddPage(self.LOPT, "LOPT")

        sizer_3 = wx.BoxSizer(wx.VERTICAL)

        self.window_2 = wx.SplitterWindow(self.LOPT, wx.ID_ANY)
        self.window_2.SetMinimumPaneSize(20)
        sizer_3.Add(self.window_2, 1, wx.EXPAND, 0)

        self.window_2_pane_1 = wx.Panel(self.window_2, wx.ID_ANY)

        sizer_7 = wx.BoxSizer(wx.HORIZONTAL)

        self.lopt_lev_ojlv = GroupListView(self.window_2_pane_1, wx.ID_ANY, style=wx.LC_REPORT|wx.SUNKEN_BORDER, showItemCounts = False)
        sizer_7.Add(self.lopt_lev_ojlv, 1, wx.EXPAND, 0)

        self.window_2_pane_2 = wx.Panel(self.window_2, wx.ID_ANY)

        sizer_8 = wx.BoxSizer(wx.HORIZONTAL)

        sizer_8.Add((0, 0), 0, 0, 0)

        self.LEVHAMS = wx.Panel(self.main_panel, wx.ID_ANY)
        self.main_panel.AddPage(self.LEVHAMS, "LEVHAMS")

        sizer_4 = wx.BoxSizer(wx.VERTICAL)

        sizer_4.Add((0, 0), 0, 0, 0)

        self.LEVHAMS.SetSizer(sizer_4)

        self.window_2_pane_2.SetSizer(sizer_8)

        self.window_2_pane_1.SetSizer(sizer_7)

        self.window_2.SplitVertically(self.window_2_pane_1, self.window_2_pane_2)

        self.LOPT.SetSizer(sizer_3)

        self.window_1_pane_2.SetSizer(sizer_5)

        self.window_1_pane_1.SetSizer(sizer_2)

        self.window_1.SplitVertically(self.window_1_pane_1, self.window_1_pane_2)

        self.STRANS.SetSizer(strans_main_sizer)

        self.panel_1.SetSizer(sizer_1)

        self.Layout()

        self.Bind(wx.EVT_BUTTON, self.on_strans_add, self.strans_add_levels_btn)
        self.Bind(wx.EVT_BUTTON, self.on_strans_del, self.strans_del_levels_btn)
        self.Bind(wx.EVT_BUTTON, self.on_strans_save, self.strans_save_levels_btn)
        self.Bind(wx.EVT_CLOSE, self.on_Exit, self)
        # end wxGlade

    def on_Open(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_Open' not implemented!")
        event.Skip()

    def on_Save(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_Save' not implemented!")
        event.Skip()

    def on_Save_As(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_Save_As' not implemented!")
        event.Skip()

    def on_Export_Linelist(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_Export_Linelist' not implemented!")
        event.Skip()

    def on_Exit(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_Exit' not implemented!")
        event.Skip()

    def on_partial_strans(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_partial_strans' not implemented!")
        event.Skip()

    def on_full_strans(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_full_strans' not implemented!")
        event.Skip()

    def on_lopt(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_lopt' not implemented!")
        event.Skip()

    def on_strans_add(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_strans_add' not implemented!")
        event.Skip()

    def on_strans_del(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_strans_del' not implemented!")
        event.Skip()

    def on_strans_save(self, event):  # wxGlade: mainWindow.<event_handler>
        print("Event handler 'on_strans_save' not implemented!")
        event.Skip()

# end of class mainWindow

class MyApp(wx.App):
    def OnInit(self):
        self.frame = mainWindow(None, wx.ID_ANY, "")
        self.SetTopWindow(self.frame)
        self.frame.Show()
        return True

# end of class MyApp

if __name__ == "__main__":
    app = MyApp(0)
    app.MainLoop()
