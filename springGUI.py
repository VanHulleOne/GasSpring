# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 09:26:40 2019

springGUI.py

@author: Myself
"""

from tkinter import Tk, Label, Button

class SpringSelector:
    def __init__(self, master):
        self.master = master
        master.title('Spring Selector')
        
        self.label = Label(master, text='Helps select a spring')
        self.label.pack()
        
        self.greet_button = Button(master, text='Greet', command=self.greet)
        self.greet_button.pack()
        
        self.close_button = Button(master, text='Close', command=master.quit)
        self.close_button.pack()
        
    def greet(self):
        print('Greetings')
        
root = Tk()
selector = SpringSelector(root)
root.mainloop()