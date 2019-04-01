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

        w=880
        h=550
        ws = root.winfo_screenwidth() # width of the screen
        hs = root.winfo_screenheight() # height of the screen
        # calculate x and y coordinates for the Tk root window
        x = (ws/2) - (w/2)
        y = (hs/2) - (h/2)

        top.geometry('%dx%d+%d+%d' % (w, h, x, y)) 
        top.title("GRACy")
        top.configure(highlightcolor="black")

        image = ImageTk.PhotoImage(file="./resources/Medicon_virus.png")

        gmail=ImageTk.PhotoImage(file='./resources/Medicon_virus_times.jpg')
        self.lab=tk.Label(image=gmail)
        self.lab.photo=gmail
        self.lab.pack()

        readsQCButton = tk.Button(top)
        readsQCButton.place(x=535,y=30,height=30,width=180)
        readsQCButton.configure(text= "Reads QC",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98")

        assemblyButton = tk.Button(top)
        assemblyButton.place(x=640,y=265,height=30,width=180)
        assemblyButton.configure(text= "De novo assembly",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98")

        gentypingButton = tk.Button(top)
        gentypingButton.place(x=600,y=412,height=30,width=150)
        gentypingButton.configure(text= "Genotyping",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98")

        annotationButton = tk.Button(top)
        annotationButton.place(x=200,y=495,height=30,width=150)
        annotationButton.configure(text= "Annotation",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98")

        dbSubmissionButton = tk.Button(top)
        dbSubmissionButton.place(x=65,y=167,height=30,width=150)
        dbSubmissionButton.configure(text= "DB submission",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98")

        snpAnalysisButton = tk.Button(top)
        snpAnalysisButton.place(x=65,y=358,height=30,width=150)
        snpAnalysisButton.configure(text= "SNP analysis",bg="#568F98",fg="white",font=("Times",14,'bold'),bd=0,highlightbackground="#568F98",highlightcolor="#568F98")


        #def make_button(imageFile,xpos,ypos,h,w):
        #    image = ImageTk.PhotoImage(file=imageFile)
        ##    b = tk.Button(top)
        #    b.config(image=image)
        #    b.image = image
        #    b.place(x=xpos,y=ypos,height=h,width=w)
            
        
        #make_button("portion.png",20,20,80,80)
        #make_button("portion.png",150,150,160,160)



if __name__ == '__main__':
    vp_start_gui()
