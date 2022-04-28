# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 09:26:40 2019

springGUI.py

@author: Myself
"""

from abc import abstractmethod
import tkinter as tk
from collections import namedtuple

import springs as sp

validation = None

DONT_MOVE = False
clickedSpring = None

OutputNT = namedtuple('OutputNT', 'heading attr')

outputCols = [
              OutputNT('Name', 'name'),
              OutputNT('OD', 'OD'),
              OutputNT('ID', 'ID'),
              OutputNT('FL', 'freeLength'),
              OutputNT('Rate', 'rate'),
              OutputNT('F@L1', lambda s: s.getForce(variableWidgets['Length 1'].value)),
              OutputNT('F@L2', lambda s: s.getForce(variableWidgets['Length 2'].value)),
              OutputNT('SSL', 'safeSolidLength'),
              OutputNT('Max Def', 'maxDeflection'),
              OutputNT('log10(FOS)', lambda s: s.getLifeFOS(variableWidgets['Min Life'].value)),
              ]


def entered(event):
    global DONT_MOVE
    DONT_MOVE = True
    
def left(event):
    global DONT_MOVE
    DONT_MOVE = False

def attrFilter(springs, widget):
    return sp.rangeAttrFilter(springs, widget.attr, widget.low, widget.high)

def attrMinMax(springs, attr):
    if springs:
        return min(s.__getattribute__(attr) for s in springs), max(s.__getattribute__(attr) for s in springs)
    return 0,0

def lengthMinMax(springs, _):
    if springs:
        return min(s.minLength for s in springs), max(s.freeLength for s in springs)
    else:
        return 0,0

def lengthFilter(springs, widget):
    if widget.value:
        return sp.lengthFilter(springs, widget.value)
    return springs

def forceMinMax(springs, _):
    if springs:
        return 0, max(s.maxDeflection * s.rate for s in springs)
    else:
        return 0,0

def forceFilter(springs, widget):
    lengthWidget = variableWidgets[widget.dependant]
    lengthWidget.update(springs)
    length = lengthWidget.value
    if length:
        return sp.forceFilter(springs, length, widget.low, widget.high)
    return springs

def strokeFilter(springs, widget):
    if widget.value:
        minF1 = variableWidgets['Force 1'].low
        maxF1 = variableWidgets['Force 1'].high
        minF2 = variableWidgets['Force 2'].low
        maxF2 = variableWidgets['Force 2'].high
        
        return sp.strokeFilter(springs, widget.value, minF1, maxF1, minF2, maxF2)
    return springs
    
def checkFilter(springs, attr, widget):    
    if all(v=='0' for v in widget.values):
        return springs
    filteredSprings = []
    for value in widget.values:
        if value != '0':
            filteredSprings.extend(sp.selectionFilter(springs, attr, value))
    return filteredSprings

def safesolidLengthFilter(springs, attr, widget):
    filteredSprings = []
    for value in widget.values:
        if value != '0':
            filteredSprings.extend(sp.selectionFilter(springs, attr, value=='True')) # Note the ==, not assignment
    return filteredSprings if filteredSprings else springs

def lifeFilter(springs, widget):
    l1 = variableWidgets['Length 1'].value
    l2 = variableWidgets['Length 2'].value
    stroke = variableWidgets['Stroke'].value
    if stroke:
        minF1 = variableWidgets['Force 1'].low
        maxF1 = variableWidgets['Force 1'].high
        minF2 = variableWidgets['Force 2'].low
        maxF2 = variableWidgets['Force 2'].high
        return sp.lifeFilter(springs, widget.value, stroke=stroke, forces = (minF1, maxF1, minF2, maxF2))
    return sp.lifeFilter(springs, widget.value, length1=l1, length2=l2)

def label_format(text, value):
    return text + ' {:0.3}'.format(float(value))

def updateEvent(event):
    selector.updateVarFrames()

def updateNoEvent():
    selector.updateVarFrames()
    
def makeEntry(parent):
    entry = tk.Entry(parent, width=6, validate='key', validatecommand=(validation, '%P'))
    entry.bind('<FocusOut>', updateEvent)
    entry.bind('<KeyRelease-Return>', updateEvent)
    entry.bind('<Enter>', entered)
    entry.bind('<Leave>', left)
    return entry
    
class Widget:
    def __init__(self, label, filterFunc,*, minMaxFunc=attrMinMax, attr=None):
        self.label = label
        self.filterFunc = filterFunc
        self.minMaxFunc = minMaxFunc
        self.attr=attr
        self.frame = None
        
    def buildFrame(self, parent):
        self.frame = tk.Frame(parent, highlightthickness=1, highlightbackground='black')
        self.lab_title = tk.Label(self.frame, text=self.label)
        self.lab_remaining = tk.Label(self.frame, text='')
        self.frame._maker = self
        
        self.lab_title.grid(row=0, column=2)
        
        self.buildVars()
        
        return self.frame
    
    @abstractmethod
    def buildVars(self):
        pass
    
    @abstractmethod
    def update(self, springs):
        pass
    
    def __repr__(self):
        return self.__class__.__name__ + '({})'.format(self.label)

class RadioWidget(Widget):
    def __init__(self, label, filterFunc, options=None):
        super().__init__(label, filterFunc)
        self.options = options
        self.value = None
        
    def buildVars(self):
        self.valueVar = tk.StringVar(value=self.options[0])
        
        self.radios = [tk.Radiobutton(self.frame, text=option,
                                      variable=self.valueVar, value=option, command=updateNoEvent)
                       for option in self.options]
        
        for col, button in enumerate(self.radios):
            button.grid(row=2, column=col)
            button.bind('<Enter>', entered)
            button.bind('<Leave>', left)
        
        self.lab_remaining.grid(row=2, column=len(self.radios), padx=10)

    
    def update(self, springs):
        self.springs = springs
        self.value = self.valueVar.get()
        
        self.filteredSprings = self.filterFunc(springs, self)
        
        self.lab_remaining['text'] = len(self.filteredSprings)
        return self.filteredSprings


class CheckWidget(Widget):
    def __init__(self, label, filterFunc, attr=None, options=None):
        super().__init__(label, filterFunc)
        self.attr = attr
        self.options = options
        self.values = []
        
    def buildVars(self):        
        self.valueVars = [tk.StringVar(value=option) for option in self.options]
        
        self.checkButts = [tk.Checkbutton(self.frame, text=option,
                                      variable=var, command=updateNoEvent,
                                      onvalue=option, offvalue=None)
                       for option, var in zip(self.options, self.valueVars)]
        
        for col, button in enumerate(self.checkButts):
            button.grid(row=2, column=col)            
            button.bind('<Enter>', entered)
            button.bind('<Leave>', left)
            button.deselect()
        
        self.lab_remaining.grid(row=2, column=len(self.checkButts), padx=10)

    
    def update(self, springs):
        self.springs = springs
        self.values = [var.get() for var in self.valueVars]
        
        self.filteredSprings = self.filterFunc(springs, self.attr, self)
        
        self.lab_remaining['text'] = len(self.filteredSprings)
        return self.filteredSprings

class SingleWidget(Widget):
    def __init__(self, label, filterFunc,*, dummy=None, **kwargs):
        super().__init__(label, filterFunc, **kwargs) 
        self.minN = 0
        self.maxN = 0
        self.value = None
        
    def buildVars(self):        
        self.lab_min = tk.Label(self.frame, text=label_format('Min:', self.minN))
        self.lab_max = tk.Label(self.frame, text=label_format('Max:', self.maxN))
        self.ent_value = makeEntry(self.frame)
        
        self.lab_min.grid(row=1, column=0)
        self.ent_value.grid(row=1, column=1)
        self.lab_max.grid(row=1, column=2, padx=10)
        self.lab_remaining.grid(row=1, column=3, padx=10)

    
    def update(self, springs):
        self.springs = springs
        self.minN, self.maxN = self.minMaxFunc(springs, self.attr)
        self.lab_min['text'] = label_format('Min:', self.minN)
        self.lab_max['text'] = label_format('Max:', self.maxN)
       
        try:
            self.value = float(self.ent_value.get())
        except ValueError:
            self.value = None

        filteredSprings = self.filterFunc(springs, self)
            
        self.lab_remaining['text'] = len(filteredSprings)
        
        return filteredSprings


class RangeWidget(Widget):
    def __init__(self, label, filterFunc,*, dependant=None, **kwargs):
        super().__init__(label, filterFunc, **kwargs)
        self.dependant = dependant        
        self.minN = 0
        self.maxN = 0
    
    def buildVars(self):
        self.lab_min = tk.Label(self.frame, text=label_format('Min:', self.minN))
        self.lab_max = tk.Label(self.frame, text=label_format('Max:', self.maxN))
        self.ent_low = makeEntry(self.frame)
        self.ent_high = makeEntry(self.frame)
        
        self.lab_min.grid(row=1, column=0)
        self.ent_low.grid(row=1, column=1)
        self.lab_max.grid(row=1, column=2, padx=10)
        self.ent_high.grid(row=1, column=3)
        self.lab_remaining.grid(row=1, column=4, padx=10)
        
        self.frame._maker = self

        return self.frame
    
    def update(self, springs):
        self.springs = springs
        self.minN, self.maxN = self.minMaxFunc(springs, self.attr)
        self.lab_min['text'] = label_format('Min:', self.minN)
        self.lab_max['text'] = label_format('Max:', self.maxN)
       
        try:
            self.low = float(self.ent_low.get())
        except ValueError:
            self.low = self.minN
        try:
            self.high = float(self.ent_high.get())
        except ValueError:
            self.high = self.maxN
        
        filteredSprings = self.filterFunc(springs, self)
            
        self.lab_remaining['text'] = len(filteredSprings)
        
        return filteredSprings

variableWidgets = {'OD': RangeWidget('OD (in)', attrFilter, attr='OD'),
              'ID': RangeWidget('ID (in)', attrFilter, attr='ID'),
              'Free Length': RangeWidget('Free Length (in)', attrFilter, attr='freeLength'),
              'Force 1': RangeWidget('Force 1 (lbs)', forceFilter, minMaxFunc=forceMinMax, dependant='Length 1'),
              'Force 2': RangeWidget('Force 2 (lbs)', forceFilter, minMaxFunc=forceMinMax, dependant='Length 2'),
              'Length 1':  SingleWidget('Length 1 (in)', lengthFilter, minMaxFunc=lengthMinMax),
              'Length 2':  SingleWidget('Length 2 (in)', lengthFilter, minMaxFunc=lengthMinMax),
              'Material': CheckWidget('Material', checkFilter, attr='material', options = ['SST', 'MW', 'HD']),
              'Ends': CheckWidget('Ends', checkFilter, attr='ends', options = ['C', 'CG']),
              'Stroke': SingleWidget('Stroke', strokeFilter, attr='maxDeflection'),
              'SafeSolidLength': CheckWidget('Safe Solid Length', safesolidLengthFilter, attr='safeSolidLength', options=['True', 'False']),
              'Min Life': RadioWidget('Min Life', lifeFilter, options = ['10^5', '10^6', 'Inf']),
              }

def only_floats(num):
    if num == '.' or num == '':
        return True
    try:
        float(num)
        return True
    except Exception:
        pass
    return False

# sampleSpring = sp.Spring(True)

class SpringSelector:
    def __init__(self, master):
        global dummy_frame, validation
        self.master = master
        self.master.title('Spring Selector')
        self.moving_widget = None
        
        self.res_window = tk.Text(self.master)
        scrollbar = tk.Scrollbar(master)
        scrollbar.grid(row=1, column=2, sticky='ns')
        self.res_window.config(yscrollcommand=scrollbar.set)
        scrollbar.configure(command=self.res_window.yview)
        
        self.var_frame = tk.Frame(self.master, highlightthickness=4, highlightbackground='black')
        self.dummy_frame = tk.Frame(self.var_frame, width=50, height=44, highlightthickness=4)#, highlightbackground='red')
        validation = self.var_frame.register(only_floats)
        
        self.frames = [rw.buildFrame(self.var_frame) for rw in variableWidgets.values()]
        
        self.updateVarFrames()
        for row, item in enumerate(self.frames):            
            item.grid(row=row, column=0, sticky='ew')
        
        self.headingText = tk.Text(self.master, height=1, width=len(outputCols)*8+8)
        self.headingText.insert(tk.END, '\t'.join(col.heading for col in outputCols))
        self.headingText.configure(state='disabled')
        
        self.offset = None
            
        self.var_frame.grid(row=0, column=0, sticky='nw', rowspan=2)
        self.headingText.grid(row=0, column=1, sticky='nwe')
        self.res_window.grid(row=1, column=1, sticky='nswe')
        
        self.master.rowconfigure(0, weight=0)
        self.master.rowconfigure(1, weight=1)
        
        self.master.columnconfigure(0, weight=0)
        self.master.columnconfigure(1, weight=1)
        
        self.isOver = False
        
        self.var_frame.bind('<Enter>', self.over)
        self.var_frame.bind('<Leave>', self.not_over)
        
        self.master.bind('<Button-1>', self.on_drag_start)
        self.master.bind('<B1-Motion>', self.on_drag_motion)
        self.master.bind('<ButtonRelease-1>', self.on_drag_stop)
        
        self.res_window.bind('<ButtonRelease-1>', self.highlightCurrentLine)
        self.res_window.tag_configure('curr_line', background='#e9e9e9')
        
    def highlightCurrentLine(self, event):
        global clickedSpring
        line = self.res_window.get('current linestart', 'current lineend')
        if line:
            springName = line.split()[0]
            clickedSpring = sp.springDict[springName]
            print(clickedSpring)
        
        
        self.res_window.tag_remove('curr_line', 1.0, 'end')
        try:
            _ = self.res_window.get(tk.SEL_FIRST, tk.SEL_LAST)
        except tk.TclError:
            self.res_window.tag_add('curr_line', 'insert linestart', 'insert lineend+1c')
        
        
    def showResults(self, springs):
        self.res_window.configure(state='normal')
        self.res_window.delete('1.0', 'end')
        for s in springs:
            strings = []
            for heading, attr in outputCols:
                if isinstance(attr, str):
                    strings.append(s.__getattribute__(attr))
                else:
                    try:
                        strings.append('{:0.3}'.format(attr(s)))
                    except TypeError:
                        strings.append('')
                    except ValueError:
                        strings.append(attr(s))
            self.res_window.insert(tk.END, '\t'.join(str(data) for data in strings) + '\n')
        self.res_window.configure(state='disabled')
        
    def updateVarFrames(self):
        springs = sp.springs
        for frame in self.frames:
            try:
                springs = frame._maker.update(springs)
            except Exception as e:
                print('Error:')
                print(e)
                raise e
        self.showResults(springs)

    
    def over(self, event):
        self.isOver = True
    
    def not_over(self, event):
        self.isOver = False
        
    def on_drag_stop(self, event):
        
        if DONT_MOVE or not self.moving_widget:
            return
        self.frames[self.offset] = self.moving_widget
        self.moving_widget.configure(highlightbackground='black')
        self.moving_widget = None
        self.updateVarFrames()
        for row, frame in enumerate(self.frames):
            frame.grid(row=row, column=0)
        self.dummy_frame.grid_remove()
        

    def on_drag_motion(self, event):
        if DONT_MOVE or not self.moving_widget:
            return
        if self.isOver:
            mouse_y = self.var_frame.winfo_pointery()
            subframe_height = self.frames[0].winfo_height()
            parent_y = self.var_frame.winfo_rooty()
            local_y = mouse_y - parent_y
            curr_offset = local_y // subframe_height
            curr_offset = clamp(curr_offset, 0, len(self.frames)-1)
            
            if curr_offset != self.offset or self.dummy_frame not in self.frames:
                self.frames[self.offset], self.frames[curr_offset] = self.frames[curr_offset], self.dummy_frame
                self.offset = curr_offset
                for row, frame in enumerate(self.frames):
                    frame.grid(row=row, column=0)

            x = self.moving_widget.winfo_x() - self.moving_widget._drag_start_x + event.x
            y = self.moving_widget.winfo_y() - self.moving_widget._drag_start_y + event.y
            self.moving_widget.place(x=x, y=y)
      
    def on_drag_start(self, event):
        
        if DONT_MOVE:
            return
        if self.isOver:
            mouse_y = self.var_frame.winfo_pointery()
            subframe_height = self.frames[0].winfo_height()
            parent_y = self.var_frame.winfo_rooty()
            local_y = mouse_y - parent_y
            self.offset = clamp(local_y // subframe_height, 0, len(self.frames)-1)
                        
            self.moving_widget = self.frames[self.offset]
            self.moving_widget.configure(highlightbackground='red')
            self.moving_widget._drag_start_x = event.x
            self.moving_widget._drag_start_y = event.y
            self.moving_widget.lift()
            
def clamp(n, minn, maxn):
    if n < minn:
        return minn
    if n > maxn:
        return maxn
    return n

        
root = tk.Tk()
selector = SpringSelector(root)
root.mainloop()
