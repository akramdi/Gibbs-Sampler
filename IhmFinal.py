# -*- coding: UTF-8 -*-

##import filedialog
from tkinter import *
from GibbsFinalVersion import GibbsSampler,readFasta
class IhmGibbs:
    
    #éléments graphiques
    root=[] #c'est la fenetre
    mainFrame=[] #c'est le cadran qui contient le champs text qui affiche les séquences
    optionFrame=[] #c'est le cadran qui affiche les inputs/parametres de départ
    resultText=[] #c'est la zone de texte qui va afficher les sequences
    
    fileToLoadEntry=[] #le champ qui contient le fichier fasta à charger
    nbOfRunsEntry=[] #le champ qui contient le nombre de runs
    iterationsEntry=[] #le nombre d'itération par runs
    motifWidthEntry=[] #champ qui contient la taille du motif
    residuesChoice=[] #choix des résidus acides aminés ou nucléotides
    #fin des éléments
    
    sequences=[] #les sequences édités
    
    residuesChoicesDic = {
        1: "ARNDCQEGHILKMFPSTWYV",
        2: "ACTG",
        }
    
    
    def __init__(self):
        self.root=Tk()
        self.root.geometry('1200x500')
        
        title='GIBBS SAMPLING'
        self.root.title(title)
        
        self.optionFrame=Frame(self.root)
        self.optionFrame.pack()
        
        self.mainFrame=Frame(self.root)
        self.mainFrame.pack()

        
        #faire les scrollbars
        xscrollbar = Scrollbar(self.mainFrame, orient=HORIZONTAL)
        xscrollbar.grid(row=1, column=0, sticky=N+S+E+W)
        
        yscrollbar = Scrollbar(self.mainFrame)
        yscrollbar.grid(row=0, column=1, sticky=N+S+E+W)
        
        #faire la textbox des résultats
        self.resultText = Text(self.mainFrame, wrap=NONE,
                xscrollcommand=xscrollbar.set,
                yscrollcommand=yscrollbar.set,width=100)
        self.resultText.grid(row=0, column=0)
        
        xscrollbar.config(command=self.resultText.xview)
        yscrollbar.config(command=self.resultText.yview)
        
        #text d'introduction
        self.resultText.insert(END, 'Hello this is a GIBBS sampling app')
        
        #faire les box entry pour nom de fichiers
        label = Label(self.optionFrame, text="Fasta file :")
        label.pack(side=LEFT)
        self.fileToLoadEntry = Entry( self.optionFrame)
        self.fileToLoadEntry.insert(0, "hth.fst")
        self.fileToLoadEntry.pack(side=LEFT)
        #bouton pour prévisualiser le fichier fasta
        b = Button(self.optionFrame, text ="Preview", command = self.readFastaFileAndShow )
        b.pack(side=LEFT)
         #bouton parcourir
        button_brow= Button(self.optionFrame, text ="browse", command = self.browser)
        button_brow.pack(side=LEFT)
         
         #entry pour le nombre de runs
        label = Label(self.optionFrame, text="Nb Runs :")
        label.pack(side=LEFT)
        self.nbOfRunsEntry = Entry( self.optionFrame)
        self.nbOfRunsEntry.insert(0, "1")
        self.nbOfRunsEntry.pack(side=LEFT)
         
        label = Label(self.optionFrame, text="Iterations :")
        label.pack(side=LEFT)
        self.iterationsEntry = Entry( self.optionFrame)
        self.iterationsEntry.insert(0, "200")
        self.iterationsEntry.pack(side=LEFT)
        
        label = Label(self.optionFrame, text="Width(W) :")
        label.pack(side=LEFT)
        self.motifWidthEntry = Entry( self.optionFrame)
        self.motifWidthEntry.insert(0, "18")
        self.motifWidthEntry.pack(side=LEFT)   

        label = Label(self.optionFrame, text="Residues :")
        label.pack(side=LEFT)
        
        #radio button pour selectionner les residues
        self.residuesChoice=IntVar()
        residueFrame=Frame(self.optionFrame)
        residueFrame.pack(side=LEFT)
        R1 = Radiobutton(residueFrame, text="Acides amines", variable=self.residuesChoice, value=1)
        R1.pack()
        R1.select()
        R2 = Radiobutton(residueFrame, text="Nucleotides", variable=self.residuesChoice, value=2)
        R2.pack()
        
        
        b = Button(self.optionFrame, text ="RunGibbs", command = self.runGibbs)
        b.pack(side=LEFT)
        
    def run(self):
        self.root.mainloop()
    
    def runGibbs(self):
        self.readFastaFileAndShow();
        g=GibbsSampler(
                self.sequences, 
                int(self.motifWidthEntry.get()), 
                self.residuesChoicesDic[self.residuesChoice.get()]
                )
        
        g.optimizeAlignment(int(self.iterationsEntry.get()), int(self.nbOfRunsEntry.get()))
        self.resultText.delete(1.0, END)
        self.resultText.insert(END,self.prettyPrintSequences(g.sequences,int(self.motifWidthEntry.get())))
    
    def readFastaFileAndShow(self):
        self.sequences=readFasta(self.fileToLoadEntry.get())
        self.resultText.delete(1.0, END)
        self.resultText.insert(END,self.prettyPrintSequences(self.sequences,int(self.motifWidthEntry.get())))
    
    def prettyPrintSequences(self,sequences,w=18):
        i=0
        maxStart=max([k['start'] for k in sequences])
        res=""
        for seq in sequences:
            startPos=seq['start']
            endpos=seq['start']+w
            res+= str(i).rjust(2,'0')+" "+seq['nom'][3:20] +"   "+ seq['seq'][:startPos].rjust(maxStart,'-')+"["+seq['seq'][startPos:endpos]+"]"+seq['seq'][endpos:]+"\n"
            i+=1
        return res
    def browser (self):
        link= filedialog.askopenfilename(parent=self.optionFrame,title='Please select ')
        self.fileToLoadEntry.delete(0,END)
        self.fileToLoadEntry.insert(0,link)

#main
ihm=IhmGibbs()
ihm.run()
