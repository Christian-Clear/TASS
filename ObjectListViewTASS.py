#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ObjectListView import ObjectListView
import itertools
import wx


class ObjectListViewTame(ObjectListView):

    def _FindByTyping(self, searchColumn, prefix):
        """
        Select the first row passed the currently focused row that has a string representation
        that begins with 'prefix' in the given column
        """
        start = max(self.GetFocusedRow(), 0)

        # If the user is starting a new search, we don't want to consider the
        # current row
        if len(prefix) == 1:
            start = (start + 1) % self.GetItemCount()

        # A binary search on a sorted column can handle any number of rows. A linear
        # search cannot. So we impose an arbitrary limit on the number of rows to
        # consider. Above that, we don't even try
        if self.GetItemCount() > self.MAX_ROWS_FOR_UNSORTED_SEARCH:
            self._SelectAndFocus(0)
            return

        # Do a linear, wrapping search to find the next match. To wrap, we consider
        # the rows in two partitions: start to the end of the collection, and then
        # from the beginning to the start position. Expressing this in other languages
        # is a pain, but it's elegant in Python. I just love Python :)
            
        for i in itertools.chain(range(start, self.GetItemCount()), range(0, start)):
            #self.__rows += 1
            model = self.GetObjectAt(i)
            if model is not None:
                strValue = searchColumn.GetStringValue(model)
                if strValue.lower().startswith(prefix.lower()):
                    self._SelectAndFocus(i)
                    return True
        return False

