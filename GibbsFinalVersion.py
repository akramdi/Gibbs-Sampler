# -*- coding: UTF-8 -*-

import random
import math


class GibbsSampler:
    """ Classe définissant une liste de séquences. Chaque séquence est un dictionnaire:
        Sequencelist=[seq1, seq2,...,seqi,..., seqn].
        seqi ={'seq':string aa ou nuc, 'nom': le nom de la seq, 'start': int représentant position du debut du motif}

        toutes les séquences sont caracterisées par :
            - MotifWidth : la longueur du motif à trouver
            - residues : les résidues des sequences (aa ou nuc)
            - pseudocount total: somme des pseudocounts de chaque résidu """


    sequences = []
    motifWidth = 0
    pseudoCounts = {}
    pseudoCountsTotal=0
    residues =""
    
    def __init__(self, sequences, motifWidth,residues):
        self.sequences = sequences
        self.motifWidth = motifWidth
        self.residues=residues
        self.calculatePseudoCounts()
        self.prettyPrintSequences()
        
    
    def prettyPrintSequences(self):
        """ Méthode permettant d'afficher joliement l'aligment local en mettant en evidence les motifs avec des '[]' """
        
        i=0
        temp=[]
        #récuperer la position du début du motif la plus elevée 
        for k in self.sequences:
            temp.append(k['start'])
        maxStart=max(temp)
        #pour chaque séquence..
        for seq in self.sequences:
            #..récuperer la position du début et fin du motif
            startPos=seq['start']
            endpos=seq['start']+self.motifWidth

            #L'idée est de faire glisser tous les segements [début séquence:début motif] vers la position max (Maxstart) en rajoutant '-' à sa doite..
            #.. +'[' pour marquer le début du motif + le segement du motif +']' pour marquer la fin du motif+ la fin de ma séquence.    
            print (str(i).rjust(3,'0')," ", seq['seq'][:startPos].rjust(maxStart,'-')+"["+seq['seq'][startPos:endpos]+"]"+seq['seq'][endpos:])
            i+=1
        
    def calculatePseudoCounts(self):
        """ Méthode permettant de calculer les pseudocount(bi) de chaque résidu selon la formule: pseudocount residu i= beta*fréquence i
            Renvoie un dictionnaire des pseudocounts de la forme : pseudoCounts={'residu i': valeur pseudocount i}"""


        #calcul nbr tot de résidus dans toutes les séquences:
        temp=[]
        for s in self.sequences:
            temp.append(len(s['seq']))

        total=float(sum(temp))
        #calcul beta : sqrt(N-1)
        beta = math.sqrt((len(self.sequences)-1))
        print ("total ", total)
        
        #calcul pseudocount
        ##compter le nbr tot de chaque résidu dans toutes les séquences:
        n=0
        for res in self.residues:
            temp2=[]
            for s in self.sequences:
                temp2.append(s['seq'].count(res))
            n=sum(temp2)    
            print (res," ", n)
            #calcul pseudocount du residu i = fréquence residu i * beta 
            self.pseudoCounts[res] = (n / total)*beta
        #calcul du pseudocount tot    
        self.pseudoCountsTotal=sum(self.pseudoCounts.values())

        print ("pseudo counts")
        print (self.pseudoCounts)
        
    def optimizeAlignment(self, maxIterations=50,runs=3):
        """ methode permettant de trouver le meilleur alignement local en mettant a jour les cle 'start' de chaque sequence
            avec la position du debut du motif la plus probable et affiche l'alignement local'"""
        
        print ("start finding motif")
        #initlialiser la valeur de F et la liste des positions 
        bestF=0
        bestStartingPos = []

        
        #boucle mère
        for m in range(runs):
            #################################################### INITIALISATION #######################################################
            print ("Randomizing starting position for the ", m, "th time...", "best score : ",bestF)
            self.randomizeStartingPositions()
            
            i = 0
            ##################################################### ITERATIONS ##########################################################
            while i<maxIterations :
                if i%10==0: #afficher ce qui suit tous les 10 tours, pour observer la convergence
                    print (str(i).rjust(4,'0'), "th Iteration ;", "bestScore : ",bestF)
                    #print "best starting positions ", bestStartingPos

                ####################PREDICTION STEP :
                    #choisir une séquence z (le choix se fait selon l'ordre de la liste)
                for indexZ in range(len(self.sequences)):
                    #définir la séquence z 
                    sequenceZ=self.sequences[indexZ]
                    #créer liste contenant toutes les seq SAUF la séquence z, en supprimant la séquence z après avoir copié toutes les séquences.
                    sequencesWithoutSequenceZ=self.sequences[:]#copie
                    del(sequencesWithoutSequenceZ[indexZ])#suppression 
                    
                    #calculer les matrices de fréquences et d'occurrences
                    occurrenceMatrix,freqMatrix=self.calculateFreqMatrix(sequencesWithoutSequenceZ)
                    backgroundFreq=self.calculateBackgroundFreq(indexZ)

                #################### SAMPLING STEP
                    #print str(séquenceZ).rjust(3,'0'), "th teration "
                    #Echantillonner une position du début du motif dans la séquence Z
                    self.sampleNewPosition(freqMatrix,backgroundFreq, sequenceZ)
                
                ###Calculer F
                #calculer les matrices d'occurences et fréquences pour l'aligement trouvé
                fullOccMatrix,fullFreqMatrix=self.calculateFreqMatrix(self.sequences)
                fullBackgoundFreq=self.calculateBackgroundFreq(-1) #je donne un argument mock, je n'ai pas de sequence z à tester.
                #calculer F
                F = self.calculateF(fullFreqMatrix,fullBackgoundFreq,fullOccMatrix)
                
                #Tester si F est maximisé en comparant avec celui calculé juste avant
                if F > bestF:
                    # si oui, mettre à jour bestF
                    bestF = F
                    bestStartingPos = []
                    # récuperer les positions échantillonnées qui ont permis le calcul de ce F max
                    for seq in self.sequences:
                        bestStartingPos.append(seq['start'])
                     
            
                i+=1
            
        
        
        #une fois la convergence atteinte, mettre à jour les key 'start' des séquences avec les positions des motifs
        for k in range(len(self.sequences)):
            self.sequences[k]['start'] = bestStartingPos[k]

        #Afficher le F max et l'alignement local
        print ("the best score is : " ,bestF )
        self.prettyPrintSequences()
                
    
    def randomizeStartingPositions(self):
        """ Méthode permettant de générer aléatoirement, pour chaque séquence, une position du début du motif.
		La position i est comprise entre 0: (longueure de la séquence i - longueure motif)
		La méthode initialise la key 'start' de chaque séquence avec la position générée. """
        
        for seq in self.sequences:
            seqSize=len(seq['seq'])
            startpos=random.randint(0,seqSize - self.motifWidth)
            seq['start']=startpos
        return
    
    def calculateBackgroundFreq(self,indexZ):
        """Méthode permettant de calculer les fréquences de chaque résidu en dehors du motif(background) pour toutes les séquences (Y COMPRIS la sequence z).
                Prends en argument l'indice de la séquences z exclue
		Renvoie un dictionnaire :backgroundFreq = {"residu j" : fréquence du residu j dans le background} """
		
        
        #initialisation de la matrice des occurrence au background avec les pseudocounts de chaque résidu 
        backgroundOcc={}
        backgroundFreq={}
        for res in self.residues:
            backgroundOcc[res]=self.pseudoCounts[res]
        #calcul nbr total des résidus dans toutes les séquences
        temp=[]
        for seq in self.sequences:
            temp.append(len(seq['seq']))
        
        NbTotalRes=float(sum(temp))
        # calcul des occurrences de chaque residu en dehors du motif. Pour la séquences exclue z, le calcul se fait sur toute sa longueur.
        i=0
        for seq in self.sequences:
            for x in range(len(seq['seq'])): #parcours de chaque séquence
                #appliquer la suite seulement si je me trouve sur la séquence z, ou (avant le motif ou apres le motif pour les séquences autres que z)
                if (indexZ == i or (x < seq['start'] or x >= ( seq['start']+ self.motifWidth))):
                    backgroundOcc[seq['seq'][x]] += 1 # calculer cji 
            i+=1
        # calcul des fréquences de chaque résidu à partir des occurrences: q0j =(cij+bi)/(nbtotresidus+somme pseudocount)
        for res in self.residues:
            backgroundFreq[res]=backgroundOcc[res]/(NbTotalRes+self.pseudoCountsTotal)
            
        return backgroundFreq




    def calculateFreqMatrix(self, sequences):
        """ Méthode permettant de calculer les fréquences de chaque résidu j à chaque position i du motif pour toutes les séquences.
	    Prend en argument la liste de sequences SANS la sequence exclue z.
	    Renvoie :
	     - un dictionnaire :occMatrix (cij)= {"residu j": [liste des occurrences du residu j a chaque postion i des motifs de toutes les séquences]}
	     - un dictionnaire :FreqMatrix (qij)= {"residu j": [liste des frequences du residu j a chaque postion des motifs de toutes les séquences]}
                """

        freqMatrix={} #q ij
        occurrenceMatrix={} #c ij
        
        #initialiser chaque key du dictionnaire occurrenceMatrix par une liste de 0 de longueur self.motifWidth
        for res in self.residues:
            occurrenceMatrix[res]=[0]* self.motifWidth
            freqMatrix[res]=[0]* self.motifWidth
            
        #compter les residues à chaque position du motif et pour toutes les séquences prises en argument:	
        for res in self.residues:
            for j in range(self.motifWidth): 
                for s in sequences:
                    position = s['start'] + j
                    if s['seq'][position] == res:
                        occurrenceMatrix[res][j] += 1
                    #calcluer la frequence de res juste après :qij= (cij+ pseudocount i)/(N-1 + somme des pseudocounts)
                    freqMatrix[res][j] = (occurrenceMatrix[res][j] + self.pseudoCounts[res])/ (float(len(sequences)-1) + self.pseudoCountsTotal)
        
        return occurrenceMatrix, freqMatrix

    
            
    def sampleNewPosition(self,freqMatrix,backgroundFreq,sequenceZ):
        """ Méthode permettant d'échantillonner la position du début du motif dans la sequence z
		Prends en argument : freqMatrix, backgroundFreq, la sequenceZ. 
		Calcule les poids Ax= Qx/Px et les garde dans la liste 'probas',


		===>>>
		Selectionne un Ax aléatoirement (apres normalisation de la distribution ), donc une position du début du motif dans la sequence z.
		 
		
		Qx : produit des frequences des residus du segement x
		Px: produit de ces residus au background """
        
        probas = []
        #parcours toutes les positions possibles du debut du motif dans la séquence z : les segements x allant de 0: L-w+1

        #pour chaque début de segement x..
        for r in range(len(sequenceZ['seq']) - self.motifWidth + 1):
            
            Qx = 1 
            Px = 1
            #...parcourir ce segement de son debut jusqu'a w...
            for x in range(self.motifWidth):
                #...Et récuperer les fréquence du residu à cette position du motif et multiplier par la freq de l'itération précédente 
                Qx*= freqMatrix[ sequenceZ['seq'][r+x] ][x]
                #Récuperer freq de ce residu au background
                Px*= backgroundFreq[sequenceZ['seq'][r+x]]
            
            probas.append(Qx /Px)
        
        sequenceZ["start"]=self.generateRandomPosition(probas)
        
        return 
        
    def generateRandomPosition(self, distribution):
        """ Méthode permettant de choisir une position selon lq distribution des poids des segements Ax """
        
        total = sum(distribution)
        #normaliser la distribution
        normalizedDistrib = [float(i)/total for i in distribution]
        #générer un reél entre 0:1
        selectedPosition = random.random()

        current = 0
        for p in range(len(normalizedDistrib)):
            #tester si le poid >= au réel génére selon une loi unifrome
            if normalizedDistrib[p] + current >= selectedPosition:
                return p #si condition vérifiée, retourner position du poid à choisir
            current += normalizedDistrib[p]
    
    def calculateF(self,freqMatrix,backgroundFreq,counts):
        """ Méthode permettant de calculer la fonction à maximiser : F 
	    F = sum (c ij * log (q ij/pj)), la somme se fait sur les i=0:motifWidth et les rédidues"""
        
        F = 0
        for residue in self.residues:
            for i in range(self.motifWidth):
                F += counts[residue][i] * math.log(
                    freqMatrix[residue][i] / backgroundFreq[residue], 2)
        return F
        
        
def readFasta(fiFasta):
    """ Fonction permettant de lire un fichier fasta.
            Prends en arguement un fichier fasta et renvoie une liste de dictionnaire , chaque dictionnaire est une séquence et contient 3 keys :
            ex : seq={'seq': sequence, 'start': a ou position du debut du motif, 'nom': le nom de la sequence}"""
    
    sequences =[]
    start =  True
    seq =''
    nom=''
    with open(fiFasta,'r') as f:
        for line in f:
            #si je suis sur la ligne du nom
                if line[0]=='>':
                    #si ce n'est pas la toute premiere ligne de nom du fichier
                    if not start:
                        #stoker la sequence lu juste avant
                        sequences.append({'seq': seq, 'start': 0, 'nom':nom})
                    else:
                        #si je suis sur la toute premiere ligne du fichier:
                        #marquer le passage
                        start = False
                    #dans les deux cas, le lis le nom et je cree une seq vide pour la concatenation 
                    seq =''
                    nom = line[1:-1]
            #si je suis sur les lignes des sequence:
                else:
                    seq +=line[:-1]
                    
    #apres le for pour la derniere sequence
        seq = seq[:-1]
        sequences.append({'seq': seq, 'start': 0, 'nom':nom})
    return sequences  

def test():
    seqFromFasta=readFasta("hth.fst")
    g=GibbsSampler(seqFromFasta, 18, "ARNDCQEGHILKMFPSTWYV")
    g.optimizeAlignment(3000, 3)

