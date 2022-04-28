# -*- coding: utf-8 -*-
"""
Created on Mon Jun 21 14:17:01 2021

testing.py

@author: lvanhulle
"""

# import tkinter as tk

# class MyApp(tk.Tk):
#     def __init__(self, *args, **kwargs):
#         tk.Tk.__init__(self, *args, **kwargs)
#         self.text = tk.Text(self)
#         self.text.pack(side="top", fill="both", expand=True)
#         self.text.tag_configure("current_line", background="#e9e9e9")
#         self._highlight_current_line()

#     def _highlight_current_line(self, interval=100):
#         '''Updates the 'current line' highlighting every "interval" milliseconds'''
#         self.text.tag_remove("current_line", 1.0, "end")
#         self.text.tag_add("current_line", "insert linestart", "insert lineend+1c")
#         self.after(interval, self._highlight_current_line)

# app = MyApp()
# app.mainloop()

import tkinter as tk

class CurrentHighlightedLineText(tk.Text):

    """Text widget with current line highlighted"""

    def __init__(self, root, *args, **kwargs):
        self.text = tk.Text.__init__(self, root, *args, **kwargs)

        self.tag_configure('currentLine', background='#e9e9e9')
        self.bind('<Key>', lambda _: self.highlightCurrentLine())
        self.bind('<Button-1>', lambda _: self.highlightCurrentLine())
        self.highlightCurrentLine(delay=0)

    def highlightCurrentLine(self, delay=10):

        def delayedHighlightCurrentLine():
            self.tag_remove('currentLine', 1.0, "end")
            self.tag_add('currentLine', 'insert linestart', 'insert lineend+1c')
            print(text.index("insert linestart", "insert lineend"))
        # This bound function is called before the cursor actually moves.
        # So delay checking the cursor position and moving the highlight 10 ms.

        self.after(delay, delayedHighlightCurrentLine)


if __name__ == "__main__":
    root = tk.Tk()

    text = CurrentHighlightedLineText(root)
    text.grid(row=0, column=0, sticky='nesw')

    root.grid_rowconfigure(0, weight=1)
    root.grid_columnconfigure(0, weight=1)

    root.mainloop()