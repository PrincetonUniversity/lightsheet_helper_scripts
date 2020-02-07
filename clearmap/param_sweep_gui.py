#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 16:47:13 2020

@author: wanglab
"""
	
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter.messagebox import showerror

class MyFrame(Frame):
    def __init__(self):
        Frame.__init__(self)
        self.master.title("Example")
        self.master.rowconfigure(5, weight=1)
        self.master.columnconfigure(5, weight=1)
        self.grid(sticky=W+E+N+S)

        self.button = Button(self, text="Browse", command=self.load_file, width=10)
        self.button.grid(row=1, column=0, sticky=W)

    def load_file(self):
        fname = askopenfilename(filetypes=(("All files", "*.*") ))
        if fname:
            try:
                print("here it comes: %s" % self.settings["template"].set(fname))
            except:                     # <- naked except is a bad idea
                showerror("Open Source File", "Failed to read file\n'%s'" % fname)
            return


if __name__ == "__main__":
    MyFrame().mainloop()

