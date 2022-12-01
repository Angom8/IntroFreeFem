from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from math import *

"Fonction de lecture du fichier .msh"
def lire_fichier_msh(nomfichier):
	fichier = open(nomfichier, "r")
	line = fichier.readline()
	first_line = line.split()
	
	#lecture first line
	mdn = int(first_line[0])
	mde = int(first_line[1])
	mda = int(first_line[2])
	
	#Lecture des différents sommets (tableau associant un sommet à D coordonnées x, y, et z)
	coord = []
	refn =  []
	for i in range(mdn) :
		line = fichier.readline()
		x, y, ref = line.split()
		coord.append([])
		coord[i].append(float(x))
		coord[i].append(float(y))
		refn.append(ref)

	#Lecture des différents triangles (tableau associant un triangle à 3 sommets s1, s2, S3)
	tri = []
	refe = []
	for i in range(mde) :
		line = fichier.readline()
		s1, s2, s3, ref = line.split()
		tri.append([int(s1)-1,int(s2)-1,int(s3)-1])
		refe.append(ref)
		
	#Lecture des différentes aretes (tableau associant une arete à 2 sommets s1 et s2)	
	aretes = []
	refa = []
	for i in range(mda) :
		line = fichier.readline()
		s1, s2, ref = line.split()
		aretes.append([int(s1)-1, int(s2)-1])
		refa.append(ref)
	
	return mda, mde, mdn, coord, tri, aretes, refe, refn, refa

"Traçage sur matplotlib du maillage brut à partir des coordonnées en utilisant triplot"
def trace_maillage(mdn, mde, mda, coord, tri, aretes):
	xs = []
	ys = []
	for i in range(mdn):
		xs.append(coord[i][0])
		ys.append(coord[i][1])
	plt.figure(2)
	plt.triplot(xs, ys,triangles=tri,scalex=4,scaley=4)

"Affichage des indicdes pour chaque sommet, arrete et triangle"
def trace_maillage_ind(mdn, mda, mde, coord, tri, aretes):

	#Pour les triangles et aretes, on fait la moyenne des coordonnées pour bien placer au milieu
	#Aretes
	for i in range(mda):
		s1 =  aretes[i][0]
		s2 = aretes[i][1]
		x = (coord[s1][0] + coord[s2][0])/2
		y = (coord[s1][1] + coord[s2][1])/2
		
		plt.text(x, y, str(i), color='red') #=> afficher du texte dans le plot

	#Sommets
	for i in range(mdn):
		plt.text(coord[i][0], coord[i][1], str(i), color='green') #=> afficher du texte dans le plot

	#Triangles
	for i in range(mde):
		s1 =  tri[i][0]
		s2 = tri[i][1]
		s3 = tri[i][2]
		x = (coord[s1][0] + coord[s2][0] + coord[s3][0])/3
		y = (coord[s1][1] + coord[s2][1]+ coord[s3][1])/3
		plt.text(x, y, str(i), color='blue')
		
"Affichage des références via les tableaux liés pour chaque sommet, arete et triangle"
def trace_maillage_ref(mdn, mda, mde, coord, tri, aretes, refe, refn, refa):

	#Pour les triangles et aretes, on fait la moyenne des coordonnées pour bien placer au milieu
	#Aretes
	for i in range(mda):
		s1 =  aretes[i][0]
		s2 = aretes[i][1]
		x = (coord[s1][0] + coord[s2][0])/2
		y = (coord[s1][1] + coord[s2][1])/2
		
		plt.text(x, y, refa[i], color='red') #=> afficher du texte dans le plot
		
	#Sommet
	for i in range(mdn):
		plt.text(coord[i][0], coord[i][1], refn[i], color='green')
		
	#Triangle
	for i in range(mde):
		s1 =  tri[i][0]
		s2 = tri[i][1]
		s3 = tri[i][2]
		x = (coord[s1][0] + coord[s2][0] + coord[s3][0])/3
		y = (coord[s1][1] + coord[s2][1]+ coord[s3][1])/3
		plt.text(x, y, refe[i], color='blue')
		
"Calcul du pas pour l'ensemble du maillage (correspondant au diamètre le plus grand de tous les triangles"
def calcul_pas(mdn, mde, tri, coord):
	dmax = 0
	for i in range(mde):
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
def calcul_qualite(mdn, mde, tri, coord):
	dmax = 0
	for i in range(mde):
	
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

nomfichier = input("Nom du fichier : ")
mda, mde, mdn, coord, tri, aretes, refe, refn, refa = lire_fichier_msh(nomfichier)
affichage = input("Affichage des indices ou des références ? (ind/ref/else)\n")
trace_maillage(mdn, mde, mda, coord, tri, aretes)

if(affichage == "ind"):
	trace_maillage_ind(mdn, mda, mde, coord, tri, aretes)
elif(affichage == "ref"):
	trace_maillage_ref(mdn, mda,mde, coord, tri, aretes, refe, refn, refa)
	
pas =  calcul_pas(mdn, mde, tri, coord)
qualite = calcul_qualite(mdn, mde, tri, coord)

plt.title("Maillage avec " + str(mde) + " triangles et " + str(mdn) + " sommets.\nPas : "+ str(pas) + ", Qualité : "+ str(qualite))
plt.show()
	

