#!/usr/bin/python
installationDirectory = "/home3/scc20x/Software/mySoftware/GRACy/"


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

from PIL import ImageTk, Image
import sys
import tkFont
import os



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

        def launchQC():
            os.system("python "+installationDirectory+"qualityCheck/readsQualityCheck.py "+installationDirectory+" &")

        def denovoAssembly():
            os.system("python "+installationDirectory+"assembly/hcmvAssembly.py "+installationDirectory+" &")

        def genotyping():
            os.system("python "+installationDirectory+"genotyping/hcmvGenotypingGUI.py "+installationDirectory+" &")

        def annotation():
            os.system("python "+installationDirectory+"annotation/hcmvAnnonation.py "+installationDirectory+" &")

        def snpAnalysis():
            os.system("python "+installationDirectory+"snpAnalysis/hcmvSNPAnalysisGUI.py "+installationDirectory+" &")

        def dbSubmission():
            os.system("python "+installationDirectory+"databaseSubmission/hcmvDataSubmission.py "+installationDirectory+" &")




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

        w=351
        h=624
        ws = root.winfo_screenwidth() # width of the screen
        hs = root.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/2) - (w/2)
        y = (hs/2) - (h/2)

        top.geometry('%dx%d+%d+%d' % (w, h, x, y))
        top.title("GRACy")
        top.configure(highlightcolor="black")

        #image = ImageTk.PhotoImage(file=installationDirectory+"resources/Medicon_virus.jpg")

        gmail=ImageTk.PhotoImage(file=installationDirectory+'resources/GRACyMainPage_Small.jpg')
        self.lab=tk.Label(image=gmail)
        self.lab.photo=gmail
        self.lab.pack()

        readsQCButton = tk.Button(top, command=launchQC)
        readsQCButton.place(x=80,y=37,height=30,width=120)
        readsQCButton.configure(text= "Filtering",bg="#0069B3",fg="white",font=("Maiandra GD",14,'bold'),bd=0,highlightbackground="#0069B3",highlightcolor="#0069B3",activebackground='#0069B3')

        assemblyButton = tk.Button(top,command=denovoAssembly)
        assemblyButton.place(x=150,y=98,height=30,width=120)
        assemblyButton.configure(text= "Assembly",bg="#0069B3",fg="white",font=("Maiandra GD",14,'bold'),bd=0,highlightbackground="#0069B3",highlightcolor="#0069B3",activebackground='#0069B3')

        gentypingButton = tk.Button(top,command=genotyping)
        gentypingButton.place(x=80,y=158,height=30,width=150)
        gentypingButton.configure(text= "Genotyping",bg="#0069B3",fg="white",font=("Maiandra GD",14,'bold'),bd=0,highlightbackground="#0069B3",highlightcolor="#0069B3",activebackground='#0069B3')


        snpAnalysisButton = tk.Button(top,command=snpAnalysis)
        snpAnalysisButton.place(x=120,y=218,height=30,width=150)
        snpAnalysisButton.configure(text= "SNP analysis",bg="#0069B3",fg="white",font=("Maiandra GD",14,'bold'),bd=0,highlightbackground="#0069B3",highlightcolor="#0069B3",activebackground='#0069B3')


        annotationButton = tk.Button(top,command = annotation)
        annotationButton.place(x=70,y=278,height=30,width=150)
        annotationButton.configure(text= "Annotation",bg="#0069B3",fg="white",font=("Maiandra GD",14,'bold'),bd=0,highlightbackground="#0069B3",highlightcolor="#0069B3",activebackground='#0069B3')

        dbSubmissionButton = tk.Button(top,command = dbSubmission)
        dbSubmissionButton.place(x=120,y=341,height=30,width=150)
        dbSubmissionButton.configure(text= "Submission",bg="#0069B3",fg="white",font=("Maiandra GD",14,'bold'),bd=0,highlightbackground="#0069B3",highlightcolor="#0069B3",activebackground='#0069B3')





if __name__ == '__main__':
    vp_start_gui()
