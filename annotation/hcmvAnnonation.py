try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True


import sys


def vp_start_gui():
    '''Starting point when module is the main routine.'''
    global val, w, root

    root = tk.Tk()

   

    top = Toplevel1 (root)
    root.mainloop()



def create_Toplevel1(root, *args, **kwargs):
    
    '''Starting point when module is imported by another program.'''
    global w, w_win, rt
    rt = root
    w = tk.Toplevel (root)
    VirHosFilt_support.set_Tk_var()
    top = Toplevel1 (w)
    VirHosFilt_support.init(w, top, *args, **kwargs)
    return (w, top)

def destroy_Toplevel1():
    global w
    w.destroy()
    w = None

class Toplevel1:
    def __init__(self, top=None):
        
        

        


        '''This class configures and populates the toplevel window.
           top is the toplevel containing window.'''
        _bgcolor = '#d9d9d9'  # X11 color: 'gray85'
        _fgcolor = '#000000'  # X11 color: 'black'
        _compcolor = '#d9d9d9' # X11 color: 'gray85'
        _ana1color = '#d9d9d9' # X11 color: 'gray85' 
        _ana2color = '#ececec' # Closest X11 color: 'gray92' 
        self.style = ttk.Style()
        if sys.platform == "win32":
            self.style.theme_use('winnative')
        self.style.configure('.',background=_bgcolor)
        self.style.configure('.',foreground=_fgcolor)
        self.style.configure('.',font="TkDefaultFont")
        self.style.map('.',background=
            [('selected', _compcolor), ('active',_ana2color)])

        w=567
        h=320
        ws = root.winfo_screenwidth() # width of the screen
        hs = root.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/2) - (w/2)
        y = (hs/2) - (h/2)

        top.geometry('%dx%d+%d+%d' % (w, h, x, y)) 
        top.title("Annotation tool")
        top.configure(highlightcolor="black")

        inputFileLabel = tk.Label(top)
        inputFileLabel.configure(text="Input file")
        inputFileLabel.place(x=20,y=20,width=70,height=20)

        inputFileEntry = tk.Entry(top)
        inputFileEntry.place(x=20,y=40,width=300,height=30)
        inputFileEntry.insert(0,"Please select a genome....")

        inputFileButton = tk.Button(top)
        inputFileButton.place(x=330,y=40,width=100,height=30)
        inputFileButton.configure(text="Open file")

        runButton = tk.Button(top)
        runButton.place(x=440,y=40,width=100,height=30)
        runButton.configure(text="Run")

        logFileLabel = tk.Label(top)
        logFileLabel.place(x=20,y=100,height=20,width=100)
        logFileLabel.configure(text="Log window")

        logFrame = tk.Frame(top)
        logFrame.place(x=20, y=120, height=180, width=520)
        logFrame.configure(relief='groove')
        logFrame.configure(borderwidth="2")
        logFrame.configure(relief='groove')
        logFrame.configure(width=125)
        logArea = tk.Text(top,state='disabled')
        logArea.place(x=25,y=125,height=170, width=510)
        logArea.configure(background="white",borderwidth=5)
        logArea.configure(selectbackground="#c4c4c4")




if __name__ == '__main__':
    vp_start_gui()
