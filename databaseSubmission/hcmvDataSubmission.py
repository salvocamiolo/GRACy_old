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
        h=500
        ws = root.winfo_screenwidth() # width of the screen
        hs = root.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/2) - (w/2)
        y = (hs/2) - (h/2)

        top.geometry('%dx%d+%d+%d' % (w, h, x, y)) 
        top.title("Data submission tool")
        top.configure(highlightcolor="black")

        runToSubmitFileLabel = tk.Label(top)
        runToSubmitFileLabel.configure(text="Run to submit")
        runToSubmitFileLabel.place(x=20,y=20,width=90,height=20)

        runToSubmitFileEntry = tk.Entry(top)
        runToSubmitFileEntry.place(x=20,y=40,width=300,height=30)
        runToSubmitFileEntry.insert(0,"Please select a file....")

        runToSubmitButton = tk.Button(top)
        runToSubmitButton.place(x=330,y=40,width=100,height=30)
        runToSubmitButton.configure(text="Open file")

        fastqFolderLabel = tk.Label(top)
        fastqFolderLabel.configure(text="Fastq folder")
        fastqFolderLabel.place(x=20,y=80,width=80,height=20)

        fastqFolderEntry = tk.Entry(top)
        fastqFolderEntry.place(x=20,y=100,width=300,height=30)
        fastqFolderEntry.insert(0,"Please select a folder....")

        fastqFolderButton = tk.Button(top)
        fastqFolderButton.place(x=330,y=100,width=100,height=30)
        fastqFolderButton.configure(text="Open folder")

        fastqInfoLabel = tk.Label(top)
        fastqInfoLabel.configure(text="Fastq info table")
        fastqInfoLabel.place(x=20,y=140,width=110,height=20)

        fastqInfoEntry = tk.Entry(top)
        fastqInfoEntry.place(x=20,y=160,width=300,height=30)
        fastqInfoEntry.insert(0,"Please select a file....")

        fastqInfoFileButton = tk.Button(top)
        fastqInfoFileButton.place(x=330,y=160,width=100,height=30)
        fastqInfoFileButton.configure(text="Open file")


        projectInfoLabel = tk.Label(top)
        projectInfoLabel.configure(text="Project info table")
        projectInfoLabel.place(x=470,y=20,width=110,height=20)

        projectInfoEntry = tk.Entry(top)
        projectInfoEntry.place(x=470,y=40,width=300,height=30)
        projectInfoEntry.insert(0,"Please select a file....")

        projectInfoButton = tk.Button(top)
        projectInfoButton.place(x=780,y=40,width=100,height=30)
        projectInfoButton.configure(text="Open file")


        sampleInfoLabel = tk.Label(top)
        sampleInfoLabel.configure(text="Samples info Table")
        sampleInfoLabel.place(x=470,y=80,width=130,height=20)

        sampleInfoEntry = tk.Entry(top)
        sampleInfoEntry.place(x=470,y=100,width=300,height=30)
        sampleInfoEntry.insert(0,"Please select a file....")

        sampleInfoButton = tk.Button(top)
        sampleInfoButton.place(x=780,y=100,width=100,height=30)
        sampleInfoButton.configure(text="Open file")


        usernameLabel = tk.Label(top)
        usernameLabel.configure(text="Username")
        usernameLabel.place(x=470,y=140,width=75,height=20)

        usernameEntry = tk.Entry(top,justify='right')
        usernameEntry.place(x=470,y=160,width=200,height=30)
        usernameEntry.insert(0,"myENA_Username")

        passwordLabel = tk.Label(top)
        passwordLabel.configure(text="Password")
        passwordLabel.place(x=680,y=140,width=75,height=20)

        passwordEntry = tk.Entry(top,justify='right')
        passwordEntry.place(x=680,y=160,width=200,height=30)
        passwordEntry.insert(0,"*************")

        createProjectchkValue = tk.BooleanVar() 
        createProjectchkValue.set(True)
        createProjectCheckButton = tk.Checkbutton(top,variable=createProjectchkValue)
        createProjectCheckButton.place(x=20,y=230,height=20,width=150)
        createProjectCheckButton.configure(text="Create new Project")
        #createProjectCheckButton.bind('<Button-1>',selectKraken)

        createSampleschkValue = tk.BooleanVar() 
        createProjectchkValue.set(True)
        createSamplesCheckButton = tk.Checkbutton(top,variable=createProjectchkValue)
        createSamplesCheckButton.place(x=200,y=230,height=20,width=150)
        createSamplesCheckButton.configure(text="Create new Project")

        testSubmissionchkValue = tk.BooleanVar() 
        testSubmissionchkValue.set(True)
        testSubmissionCheckButton = tk.Checkbutton(top,variable=testSubmissionchkValue)
        testSubmissionCheckButton.place(x=380,y=230,height=20,width=150)
        testSubmissionCheckButton.configure(text="Test submission")


        runButton = tk.Button(top)
        runButton.place(x=780,y=220,width=100,height=30)
        runButton.configure(text="Submit")
       

        logFileLabel = tk.Label(top)
        logFileLabel.place(x=20,y=270,height=20,width=100)
        logFileLabel.configure(text="Log window")

        logFrame = tk.Frame(top)
        logFrame.place(x=20, y=290, height=180, width=860)
        logFrame.configure(relief='groove')
        logFrame.configure(borderwidth="2")
        logFrame.configure(relief='groove')
        logFrame.configure(width=125)
        logArea = tk.Text(top,state='disabled')
        logArea.place(x=25,y=295,height=170, width=850)
        logArea.configure(background="white",borderwidth=5)
        logArea.configure(selectbackground="#c4c4c4")


        

if __name__ == '__main__':
    vp_start_gui()
