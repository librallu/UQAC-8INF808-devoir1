#ifndef __ENTETE_H_
#define __ENTETE_H_

#include <cstdio>
#include <cstdlib> 
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <ctime>  
#include <cmath>
#include <vector>
#include <limits>
#include <cfloat>
using namespace std;

enum eProb	{ALPINE, BANANE, EGGHOLDER};

/* POUR TRAITER MAXCUT: ENLEVER LES COMMENTAIRES ET RETIRER L'AUTRE ENREGISTREMENT DE tProblem
struct tArc
{
	int Ni;										//Noeud d�part de l'arc NB: 0 � NbNoeud-1 
	int Nj;										//Noeud arriv�e de l'arc NB: 0 � NbNoeud-1
	int Poids;									//Poids de l'arc (+1, 0, -1)
};

struct tProblem									//**D�finition du probl�me de MAXCUT:
{
	eProb	Fonction;							//**Nom de la fonction ou du probl�me � traiter
	int		D;									//**Nbre de Variables (Nombre de noeuds pour le MAX CUT)
	double	Xmin;								//**Domaine des variables: valeur minimale 
	double	Xmax;								//**Domaine des variables: valeur maximale
	std::string Nom;							//**Nom du fichier de donn�es
	int NbNoeud;								//**Nbre de noeuds dans l'instance NB: num�rot�s de 0 � NbNoeud-1
	int NbArc;									//**Nbre d'arcs dans l'instance  NB: num�rot�s de 0 � NbArc-1
	std::vector <tArc> Arc;						//**D�finition de chaque arc. NB: Tableaux de 0 � NbArc-1. 
};*/

struct tProblem									//**D�finition pour fonction continue:
{
	eProb	Fonction;							//**Nom de la fonction ou du probl�me � traiter
	int		D;									//**Nbre de Variables (dimensions)
	double	Xmin;								//**Domaine des variables: valeur minimale 
	double	Xmax;								//**Domaine des variables: valeur maximale
};

struct tPosition								//**D�finition de la position d'une particule: 
{
	std::vector<double>		X;					//**Position actuelle de la particule pour chacune des dimensions
	double					FctObj;				//**Valeur de la fonction objectif
};

struct tParticule								//**D�finition d'une solution: 
{
	tPosition				Pos;				//**Position actuelle de la particule
	std::vector<double>		V;					//**Vitesse actuelle de la particule pour chacune des dimensions
	std::vector<tParticule*> Info;				//**Liste des informatrices de la particule
	tPosition				BestPos;			//**Meilleure position de la particule (exp�rience)
};

struct tPSO
{
	int		Taille;						//**Taille de l'essaim (nombre de particules)
	double	C1;							//**Coefficient de confiance en soi
	double	C2;							//**Coefficient de confiance en son exp�rience
	double	C3;							//**Coefficient de confiance en sa meilleure informatrice
	int		NbInfo;						//**Nombre d'informatrices pour chaque particule
	double	Vmax;						//**Vitesse maximale								//### Pour MAXSAT
	int		Iter;						//**Compteur du nombre de g�n�rations				
	double	Precision;					//**Niveau de pr�cision souhait� pour les probl�mes � variables continues	//### Pour fonctions continues
	int		CptEval;					//**COMPTEUR DU NOMBRE DE SOLUTIONS EVALUEES. SERT POUR CRITERE D'ARRET.
	int		NB_EVAL_MAX;				//**CRITERE D'ARRET: MAXIMUM "NB_EVAL_MAX" EVALUATIONS.
};

#endif
