from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from math import *
from scipy import linalg

################# IMPORT TP1 #################"""""

"Fonction de lecture du fichier .msh"
def lire_fichier_msh(nomfichier):
	fichier = open(nomfichier, "r")
	line = fichier.readline()
	first_line = line.split()
	
	#lecture first line
	nbn = int(first_line[0])
	nbe = int(first_line[1])
	nba = int(first_line[2])
	
	#Lecture des différents sommets (tableau associant un sommet à D coordonnées x, y, et z)
	coord = []
	refn =  []
	for i in range(nbn) :
		line = fichier.readline()
		x, y, ref = line.split()
		coord.append([])
		coord[i].append(float(x))
		coord[i].append(float(y))
		refn.append(ref)

	#Lecture des différents triangles (tableau associant un triangle à 3 sommets s1, s2, S3)
	tri = []
	refe = []
	for i in range(nbe) :
		line = fichier.readline()
		s1, s2, s3, ref = line.split()
		tri.append([int(s1)-1,int(s2)-1,int(s3)-1])
		refe.append(ref)
		
	#Lecture des différentes aretes (tableau associant une arete à 2 sommets s1 et s2)	
	aretes = []
	refa = []
	for i in range(nba) :
		line = fichier.readline()
		s1, s2, ref = line.split()
		aretes.append([int(s1)-1, int(s2)-1])
		refa.append(ref)
	
	return nba, nbe, nbn, coord, tri, aretes, refe, refn, refa
	
"Calcul du pas pour l'ensemble du maillage (correspondant au diamètre le plus grand de tous les triangles"
def calcul_pas(nbn, nbe, tri, coord):
	dmax = 0
	for i in range(nbe):
		s1 = tri[i][0]
		s2 = tri[i][1]
		s3 = tri[i][2]
		
		#Un peu de pythagore pour les distances (il faudrait aussi considérer les vecteurs unitaires)
		l1 = sqrt( (coord[s2][0]-coord[s1][0])**2 +  (coord[s2][1]-coord[s1][1])**2  ) #s1s2
		l2 = sqrt( (coord[s3][0]-coord[s2][0])**2 +  (coord[s3][1]-coord[s2][1])**2  ) #s2s3
		l3 = sqrt( (coord[s1][0]-coord[s3][0])**2 +  (coord[s1][1]-coord[s3][1])**2  ) #s3s1
		
		current_d = max(l1, l2, l3)
		if(current_d > dmax):
			dmax = current_d
		
	return dmax

"Calcul de la qualité du maillage à partir des pas et rayons des cercles inscrits de chaque triangle"
def calcul_qualite(nbn, nbe, tri, coord):
	dmax = 0
	for i in range(nbe):
	
		s1 = tri[i][0]
		s2 = tri[i][1]
		s3 = tri[i][2]
		
		#calcul de h_t, le pas local
		l1 = sqrt( (coord[s2][0]-coord[s1][0])**2 +  (coord[s2][1]-coord[s1][1])**2  ) #s1s2
		l2 = sqrt( (coord[s3][0]-coord[s2][0])**2 +  (coord[s3][1]-coord[s2][1])**2  ) #s2s3
		l3 = sqrt( (coord[s1][0]-coord[s3][0])**2 +  (coord[s1][1]-coord[s3][1])**2  ) #s3s1
		
		h_t = max(l1, l2, l3)

		#calcul de r_t le rayon du cercle inscrit
		d_t = 0.5*(l1 + l2 + l3)
		mes = abs(0.5 * ((coord[s2][0]-coord[s1][0])*(coord[s3][1]-coord[s1][1]) - (coord[s3][0]-coord[s1][0]) *  (coord[s2][1]-coord[s1][1])  ))
		
		r_t = mes/d_t
		
		#calcul de la qualite du triangle local
		if(r_t > 0):
			current_d = (sqrt(3)/6) * (h_t/r_t)
		else :
			print("ALERTE val 0" + str(mes))
			current_d = 0
		if(current_d > dmax):
			dmax = current_d 

		
	return dmax
		
########################## Projet ###################

#On n'a pas la fonction exacte u

def fct_uE(x,y,v=500):
	"La température extérieure de uE(x,y), au bord de Fourier-Robin"
	return 0.001*(v**2);
	
def fct_f(x,y):
	"La fonction source de chaleur f(x,y) = pi**2/4 * sin(pi/2*x) + (pi**2/4 * x**3 - pi**2 * x - 2) * cos (pi/2 * y)"
	return (pi**2)/4.0 * sin(pi/2.0*x) + ((pi**2)/4.0 * (x**2) - (pi**2) * x - 2) * cos (pi/2.0 * y)

def fct_kappa(x=0,y=0,cond_mat=0.188):
	"La fonction de conductivité k(x,y). Dépend du matériau, ici l'amiante par défaut."
	return 0.188;

def fct_e(x,y):
	"L'épaisseur du bouclier en mètres"
	return 1.5

def fct_alpha(x=0,y=0):
	"Le facteur de transfert a(x,y) = 1e+8 ici (au bord de Fourier Robin)"
	return fct_kappa(x,y)/fct_e(x,y)

def coeffelem_P1_rigid(coord, tri, indice_triangle):
	"La matrice kl = (k1 ij), i et j < 3, pour le triangle T"
	
	#On note les 3 sommets de Tl
	M1 = tri[indice_triangle][0]
	M2 = tri[indice_triangle][1]
	M3 = tri[indice_triangle][2]
	
	#pour simplifier un peu
	x1 = coord[M1][0]
	x2 = coord[M2][0]
	x3 = coord[M3][0]
	
	y1 = coord[M1][1]
	y2 = coord[M2][1]
	y3 = coord[M3][1]
	
	k = np.zeros((3,3)) #matrice de 3x3
	
	#fct_kappa reprend le calcul de kappa
	kappa = fct_kappa()
	
	#mes(T) (on reprend depuis le TP1)
	mes_t = 0.5*abs((x2-x1)*(y3-y2)-(x3-x2)*(y2-y1))
	
	#définir les k
	k[0,0] = (fct_kappa()/(4*mes_t)) * ((x2-x3)**2 + (y2-y3)**2)
	
	k[1,1] = (fct_kappa()/(4*mes_t)) * ((x3-x1)**2 + (y3-y1)**2)
	
	k[2,2] = (fct_kappa()/(4*mes_t)) * ((x1-x2)**2 + (y1-y2)**2)

	k[0,1] = (fct_kappa()/(4*mes_t)) * ((-(x1-x3)*(x2-x3)) - ((y1-y3)*(y2-y3)))
	k[1,0] = k[0,1]
	
	k[0,2] = (fct_kappa()/(4*mes_t)) * ((-(x3-x2)*(x1-x2)) - ((y3-y2)*(y1-y2)))
	k[2,0] = k[0,2]
	
	k[1,2] = (fct_kappa()/(4*mes_t)) * ((-(x2-x1)*(x3-x1)) - ((y2-y1)*(y3-y1)))
	k[2,1] = k[1,2]
	
	return k
	
def coeffelem_P1_source(coord, tri, indice_triangle):
	"Le vecteur fl = (fl i), i<3 pour le triangle T"
	
	#On note les 3 sommets de Tl
	M1 = tri[indice_triangle][0]
	M2 = tri[indice_triangle][1]
	M3 = tri[indice_triangle][2]

	#pour simplifier un peu
	x1 = coord[M1][0]
	x2 = coord[M2][0]
	x3 = coord[M3][0]
	
	y1 = coord[M1][1]
	y2 = coord[M2][1]
	y3 = coord[M3][1]

	#mes(T) (on reprend depuis le TP1)
	mes_t = 0.5*abs((x2-x1)*(y3-y2)-(x3-x2)*(y2-y1))
	
	#On prend ici le milieu du triangle, soit la moyenne des coordonnées des 3 sommets (là où positionnait l'information de l'indice au TP1)
	x = (x1 + x2 + x3)/3.0
	y = (y1 + y2+ y3)/3.0
	f = fct_f(x,y)
	
	#et le vecteur
	v = [1,1,1]
	retour = []
	
	for i in range(0,3):
		retour.append(mes_t/3.0*f*v[i])
	#on multiplie le tout
	return retour

def coeffelem_P1_transf(coord, M1, M2):
	"Le vecteur ea = (ea i), i < 2, pour l'arete A"
	#La mesure en 1d d'une arete est sa longueur
	#pour simplifier un peu
	x1 = coord[M1][0]
	x2 = coord[M2][0]
	
	y1 = coord[M1][1]
	y2 = coord[M2][1]

	mes_a = sqrt( abs((x2-x1)**2 +  (y2-y1)**2 ))
		
	partie = mes_a/2.0

	#point milieu
	x = (x1 + x1)/2.0
	y = (y1 + y2)/2.0

	#alpha est constant ici
	alpha = fct_alpha(x,y)
	
	#et le vecteur
	v = [1,1]
	e = []
	for i in range(0,2):
		e.append(partie * alpha * fct_uE(x,y) * v[i])
	return e
	
def coeffelem_P1_poids(coord, M1, M2):
	"La matrice pa == (pa ij), ij < 2, pour l'arete A"
	
	#La mesure en 1d d'une arete est sa longueur
	#pour simplifier un peu
	x1 = coord[M1][0]
	x2 = coord[M2][0]
	
	y1 = coord[M1][1]
	y2 = coord[M2][1]

	mes_a = sqrt( abs((x2-x1)**2 +  (y2-y1)**2 ))
		
	partie = mes_a/6.0

	#point milieu
	x = (x1 + x1)/2.0
	y = (y1 + y2)/2.0
	
	#alpha est constant ici
	alpha = fct_alpha(x,y)
	
	#et le vecteur/matrice
	v = [[2,1],[1,2]]
	
	p = [[0,0],[0,0]]
	for i in range(0,2):
        	for j in range(0,2):
        	 p[i][j] = partie * alpha * v[i][j]
	return p

def assemblage_EF_P1(coord, tri, aretes, nbn, nbe, nba,refa):
	"Affectation de la matrice EF - P1 A et du second membre F"
	
	#Etape 1 : Mise à zéros
	A = np.zeros((nbn,nbn)) #matrice de n*n
	F = np.zeros(nbn)
	
	#Etape 2 : Addition des termes volumiques
	for i in range(0,nbe):
	
		k = coeffelem_P1_rigid(coord, tri, i)
		f = coeffelem_P1_source(coord, tri, i)
		I_un = tri[i][0]
		I_deux = tri[i][1]
		I_trois = tri[i][2]

		F[I_un] += f[0]
		F[I_deux] += f[1]
		F[I_trois] += f[2]
		
		#ligne 1
		A[I_un][I_un] += k[0][0]
		A[I_un][I_deux] += k[0][1]
		A[I_un][I_trois] += k[0][2]
		
		#ligne 2
		A[I_deux][I_un] += k[1][0]
		A[I_deux][I_deux] += k[1][1]
		A[I_deux][I_trois] += k[1][2]	
		
		#ligne 2
		A[I_trois][I_un] += k[2][0]
		A[I_trois][I_deux] += k[2][1]
		A[I_trois][I_trois] += k[2][2]
	
	K = np.copy(A) #on sauvegarde la matrice à cette étape
	
	#Etape 3 : Addition des termes de bord de Fourier/Robin
	#ici, nba == gammaF si on ne prend que les bords de Fourier/Robin, soit ceux marqués à 1
	for i in range(nba):
		if(int(refa[i]) == 1):
			I_un = aretes[i][0]
			I_deux = aretes[i][1]
			e = coeffelem_P1_transf(coord, I_un, I_deux)
			p = coeffelem_P1_poids(coord, I_un, I_deux)
				
			F[I_un] += e[0]
			F[I_deux] += e[1]
			
			#ligne 1
			A[I_un][I_un] += p[0][0]
			A[I_un][I_deux] += p[0][1]
			
			#ligne 2
			A[I_deux][I_un] += p[1][0]
			A[I_deux][I_deux] += p[1][1]
		
	return A, F, K
	
def validation_pen(coord, tri, aretes, nbn, nbe, nba, refa):
	A, F, K = assemblage_EF_P1(coord, tri, aretes, nbn, nbe, nba, refa)
	print("=========== validation_pas_a_pas ===========")
	print("nbn =")
	print(nbn)
	print("nbe =")
	print(nbe)
	print("nba =")
	print(nba)
	print("A =")
	print(A)
	print("F =")
	print(F)
	
	Uh = np.linalg.solve(A,F)
	print("Uh =")
	print(Uh)

	print("___---***   RESULTATS  ***---___")
	print("--------------------------------")
	print("{    min(Uh) : %.2f" % min(Uh))
	print("{    max(Uh) : %.2f" % max(Uh))
	print("{    max(Uh) : %.2f" % mean(Uh))
	print("{       h    : %.3f" % calcul_pas(nbn, nbe, tri, coord))
	print("{       Q    : %.3f" % calcul_qualite(nbn, nbe, tri, coord))
	
	print("--------------------------------")
	print("=========== validation_pas_a_pas ===========\n")

nomfichier = input("Nom du fichier : ")
nba, nbe, nbn, coord, tri, aretes, refe, refn, refa = lire_fichier_msh(nomfichier)
validation_pen(coord, tri, aretes, nbn, nbe, nba, refa)
