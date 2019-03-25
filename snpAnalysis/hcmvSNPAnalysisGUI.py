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

        w=900
        h=430
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
        inputFileEntry.insert(0,"Please select a file....")

        inputFileButton = tk.Button(top)
        inputFileButton.place(x=330,y=40,width=100,height=30)
        inputFileButton.configure(text="Open file")

        posToPlotLabel = tk.Label(top)
        posToPlotLabel.configure(text="Positions to plot")
        posToPlotLabel.place(x=470,y=20,width=110,height=20)

        posToPlotEntry = tk.Entry(top)
        posToPlotEntry.place(x=470,y=40,width=300,height=30)
        posToPlotEntry.insert(0,"Please select a file....")

        posToPlotButton = tk.Button(top)
        posToPlotButton.place(x=780,y=40,width=100,height=30)
        posToPlotButton.configure(text="Open file")


        cdsFileLabel = tk.Label(top)
        cdsFileLabel.configure(text="CDS file")
        cdsFileLabel.place(x=20,y=80,width=60,height=20)

        cdsFileEntry = tk.Entry(top)
        cdsFileEntry.place(x=20,y=100,width=300,height=30)
        cdsFileEntry.insert(0,"Please select a file....")

        cdsFileButton = tk.Button(top)
        cdsFileButton.place(x=330,y=100,width=100,height=30)
        cdsFileButton.configure(text="Open file")


        gffFileLabel = tk.Label(top)
        gffFileLabel.configure(text="Annotation file")
        gffFileLabel.place(x=470,y=80,width=110,height=20)

        gffFileEntry = tk.Entry(top)
        gffFileEntry.place(x=470,y=100,width=300,height=30)
        gffFileEntry.insert(0,"Please select a file....")

        gffFileButton = tk.Button(top)
        gffFileButton.place(x=780,y=100,width=100,height=30)
        gffFileButton.configure(text="Open file")


        freqCutoffLabel = tk.Label(top)
        freqCutoffLabel.configure(text="Freq cutoff")
        freqCutoffLabel.place(x=20,y=140,width=80,height=20)

        freqCutoffEntry = tk.Entry(top,justify='right')
        freqCutoffEntry.place(x=20,y=160,width=100,height=30)
        freqCutoffEntry.insert(0,"0.1")

        numThreadsLabel = tk.Label(top)
        numThreadsLabel.configure(text="Num. threads")
        numThreadsLabel.place(x=150,y=140,width=90,height=20)

        numThreadsEntry = tk.Entry(top,justify='right')
        numThreadsEntry.place(x=150,y=160,width=100,height=30)
        numThreadsEntry.insert(0,"8")
        
        runButton = tk.Button(top)
        runButton.place(x=330,y=160,width=100,height=30)
        runButton.configure(text="Run")
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        #runButton = tk.Button(top)
        #runButton.place(x=440,y=40,width=100,height=30)
        #runButton.configure(text="Run")


        logFileLabel = tk.Label(top)
        logFileLabel.place(x=20,y=200,height=20,width=100)
        logFileLabel.configure(text="Log window")

        logFrame = tk.Frame(top)
        logFrame.place(x=20, y=220, height=180, width=860)
        logFrame.configure(relief='groove')
        logFrame.configure(borderwidth="2")
        logFrame.configure(relief='groove')
        logFrame.configure(width=125)
        logArea = tk.Text(top,state='disabled')
        logArea.place(x=25,y=225,height=170, width=850)
        logArea.configure(background="white",borderwidth=5)
        logArea.configure(selectbackground="#c4c4c4")




if __name__ == '__main__':
    vp_start_gui()
