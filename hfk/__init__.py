from hfk._hfk import HFKhat
from sys import stdout

from Tkinter import Tk, Frame, Canvas, Button, IntVar, Text
from Tkinter import END, ACTIVE, DISABLED

__all__ = ['HFK', 'TkHFK']

class HFK:

    """
    Container for the mod 2 HFK^ invariant of a link.
    
    The ranks are saved in the attribute dictionary hfk_ranks.  For an
    n-component link the keys are tuples of multigradings (A, M1, ...,  Mn)
    and the value is the rank of the corresponding HFK hat group.
    (Currently, this only works for knots.)
    """
    
    def __init__(self, Xlist, Olist, name=None):
        """
        Initialize an HFK object.
        
        The Xlist and Olist arguments are sequences describing a grid
        presentation of a link.  The HFK ranks are not computed during
        initialization, since the computation can take a long time,
        and one might need to prepare for it.
        """
        self.Xlist = Xlist
        self.Olist = Olist
        self.name = name
        self.hfk_dict = None
        self.aborted = 0
                           
    def __repr__(self):
        """
        Return the diagram of Xs and Os, as a string.
        """
        size = len(self.Xlist)
        result = ""
        for i in range(size):
            for j in range(size):
                if self.Xlist[i] == j:
                    result += " X "
                elif self.Olist[i] == j:
                    result += " O "
                else:
                    result += " . "
            result += "\n"
        return result

    def HFK_ranks(self):
        """
        Print the matrix of HFK^ ranks.
         
        Calls compute_HFK_ranks if necessary.
        """
        if self.hfk_dict == None:
            self.compute_HFK_ranks()
            if self.aborted:
                print "Computation aborted!"
                self.aborted = 0
                self.hfk_dict = None
        if self.hfk_dict:
            print self.HFK_ranks_str()

    def compute_HFK_ranks(self):
        """
        Fill in the hfk_dict dictionary.
        """
        self.hfk_dict = HFKhat(self.Xlist, self.Olist, self.progress)
        if len(self.hfk_dict) == 0:
            self.Mmax = self.Mmin = self.Amax = self.Amin = 0
        else:
            gradings = self.hfk_dict.keys()
            self.Mmax = max([k[0] for k in gradings])
            self.Mmin = min([k[0] for k in gradings])
            self.Amax = max([k[1] for k in gradings])
            self.Amin = min([k[1] for k in gradings])
        self.finished()

    def finished(self):
        """
        Overide to set a semaphore when the computation is run in
        a separate thread.
        """
        return
    
    def HFK_ranks_str(self):
        """
        Return string displaying the array of ranks of mod 2 HFK^ groups.
        """
        matrix = "\n"
        if self.Amax - self.Amin + self.Mmax - self.Mmin == 0:
            matrix += "0"
        for a in range(self.Amax, self.Amin - 1, -1):
            for m in range(self.Mmin, self.Mmax + 1):
                try:
                    matrix += "%3d"%self.hfk_dict[(m,a)]
                except:
                    matrix += "  0"
            matrix += "\n"
        return matrix


    def progress(self, message, percent):
        """
        Callback function to be passed to hfk._hfk.HFKhat.

        Override this method to show progress your own way.
        """
        if percent == 0:
            print message + "     ",
            stdout.flush()
        elif percent == -1:
            print message
        else:
            print "\b\b\b\b\b%3d%%"%percent,
            stdout.flush()
        return self.aborted

class TkHFK(HFK):

    """
    Subclass of HFK that uses a TkInter progress dialog.
    """

    def compute_HFK_ranks(self):
        """
        Fill in the hfk_dict dictionary.
        """
        self.window = Tk()
        if self.name:
            self.window.title(self.name)
        else:
            self.window.title('HFK')
        self.bar = Progressbar(self.window)
        self.bar.pack(padx=10, pady=10)
        self.text = Text(self.window, font='Courier 16',
                         width=40, height=18, padx=10)
        self.text.pack(fill='x', padx=10, pady=10)
        button = Button(self.window, text="Stop", default=ACTIVE,
                        command=self.abort) 
        button.pack()
        self.window.after(10, lambda: HFK.compute_HFK_ranks(self))
        self.done = IntVar(self.window)
        self.done.set(0)
        self.window.waitvar(self.done)
        button.config(state=DISABLED)

    def HFK_ranks(self):
        if self.hfk_dict == None:
            self.compute_HFK_ranks()
            if self.aborted:
                self.text.insert(END, 'Computation aborted!\n')
                self.aborted = 0
                self.hfk_dict=None
        if self.hfk_dict:
            self.text.insert(END, 'X=%s\n'%str(self.Xlist).replace(', ',','))
            self.text.insert(END, 'O=%s\n'%str(self.Olist).replace(', ',','))
            self.text.insert(END, '\nMatrix of HFK^ ranks:\n')
            self.text.insert(END, self.HFK_ranks_str()+'\n')

    def finished(self):
        self.done.set(1)

    def abort(self):
        self.aborted=1
        
    def progress(self, message, percent):
        if percent < 0:
            self.bar.set("Done.")
            self.text.insert(END, message+'\n')                
        else:
            self.bar.set(message, percent)
        self.window.update()
        return self.aborted


class Progressbar(Frame):
    def __init__(self, parent, width=300, height=18):
        self.parent = parent
        Frame.__init__(self, parent, relief='ridge', bd=3)
        self.canvas = Canvas(self, width=width, height=height, highlightthickness=0)
        self.canvas.pack(fill='both', expand=1)
        self.bar = self.canvas.create_rectangle(0, 0, 0, 0, fill='orchid1', width=0)
        self.text = self.canvas.create_text(width/2, 2+height/2, text='')
        self.set(0,'')

    def set(self, message=None, percent=0):
        percent = max(percent, 0)
        self.percent = min(percent, 100)
        if message == None:
            message = str(self.percent) + '%'
        self.canvas.coords(self.bar, 0, 0,
                           self.canvas.winfo_width()*self.percent/100,
                           self.canvas.winfo_height())
        self.canvas.itemconfigure(self.text, text=message)
