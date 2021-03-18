#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import six
from ObjectListView import GroupListView, FastObjectListView



class GroupListViewTame(GroupListView):
    """Custom class for TASS application. Removes bugs with group sorting order when the user has selected an
    ALwaysSortByColumn. Chnages to the original GroupListView are highlighted with # CC."""

    
    def SortGroups(self, groups=None, ascending=None):
        """
        Sort the given collection of groups in the given direction (defaults to ascending).

        The model objects within each group will be sorted as well
        """
        if groups is None:
            groups = self.groups
        if ascending is None:
            ascending = self.sortAscending

        # If the groups are locked, we sort by the sort column, otherwise by the grouping column.
        # The primary column is always used as a secondary sort key.
        if self.GetAlwaysGroupByColumn():
            sortCol = self.GetSortColumn()
        else:
            sortCol = self.GetGroupByColumn()

        # Sorting event wasn't handled, so we do the default sorting
        def _getLowerCaseKey(group):
            try:
                return group.key.lower()
            except:
                return group.key
        
        if self.sortColumnIndex == self.alwaysGroupByColumnIndex:  # CC - added to stop sorting of group column when other column headers pressed
            
            if six.PY2:
                groups.sort(key=_getLowerCaseKey, reverse=(not ascending))
            else:
                groups = sorted(
                    groups,
                    key=_getLowerCaseKey,
                    reverse=(
                        not ascending))

        # Sort the model objects within each group.
        for x in groups:
            self._SortObjects(x.modelObjects, sortCol, self.GetPrimaryColumn())
            
        self.groups = groups  # CC - otherwise the listctrl isn't updated
            
        
    def _HandleColumnClick(self, evt):
        """
        The user has clicked on a column title
        """
        if self.GetAlwaysGroupByColumn == None:  # CC - therefore no GroupBy column has been set
            # If they click on a new column, we have to rebuild our groups
            if evt.GetColumn() != self.sortColumnIndex:
                self.groups = None

        FastObjectListView._HandleColumnClick(self, evt)
        
        
    def SetObjects(self, modelObjects, preserveSelection=False):
        """
        Set the list of modelObjects to be displayed by the control.
        """
        # Force our groups to be rebuilt, if we are supposd to be showing them
        if self.showGroups:
            self.groups = None
        else:
            self.groups = list()
            
        
        FastObjectListView.SetObjects(self, modelObjects, preserveSelection)        
        
        
        # CC - added to ensure that groups are sorted straight after being added to the GLV.
        self.SortBy(self.alwaysGroupByColumnIndex, True)
        self._FormatAllRows()
        
    def SearchColumn(self, searchColumn, searchString):
        """New function to search column for search string"""
        self._FindByTyping(searchColumn, searchString)
            
            