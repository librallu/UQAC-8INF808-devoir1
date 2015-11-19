#include "Entete.h"

//%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT: %%%%%%%%%%%%%%%%%%%%%%%%% 
//Le fichier de probleme (.txt) doit se trouver dans le r�pertoire courant du projet 
//pour une ex�cution � l'aide du compilateur. Indiquer les arguments du programme dans les 
//propri�t�s du projet - d�bogage - arguements.
//Sinon, utiliser le r�pertoire execution.

//*****************************************************************************************
// Prototype des fonctions
//*****************************************************************************************
void InitialisationIntervalleVariable(tProblem &unProb);
tPosition InitialisationEssaim(std::vector<tParticule> &unEssaim, tProblem unProb, tPSO &unPSO);
void InitialisationInformatrices(std::vector<tParticule> &unEssaim, tPSO &unPSO);
void InitialisationPositionEtVitesseAleatoire(tParticule &Particule, tProblem unProb);
void EvaluationPosition(tPosition &Pos, tProblem unProb, tPSO &unPSO);
void DeplacerEssaim(std::vector<tParticule> &unEssaim, tProblem unProb, tPSO &unPSO, tPosition &Meilleure);
void DeplacerUneParticule(tParticule &Particule, tProblem unProb, tPSO &unPSO);
tParticule* TrouverMeilleureInformatrice(tParticule &Particule, tPSO &unPSO);
double AleaDouble(double a, double b);
void AfficherSolutions(std::vector<tParticule> unEssaim, int Debut, int Fin, int Iter, tProblem unProb);
void AfficherUneSolution(tPosition P, int Iter, tProblem unProb);
void AfficherResultats (tPosition uneBest, tProblem unProb, tPSO unPSO);
void LibererMemoireFinPgm (std::vector<tParticule> &unEssaim, tProblem unProb, tPSO unPSO);
/* ### POUR MAXCUT
void LectureProblemeMAXCUT(std::string FileName, tProblem & unProb, tPSO &unPSO);
void AfficherProblemeMAXCUT (tProblem unProb);
*/
//******************************************************************************************
// Fonction main
//*****************************************************************************************
int main(int NbParam, char *Param[])
{
	tProblem LeProb;					//**D�finition de l'instance de probl�me
	tPSO LePSO;							//**D�finition des param�tres du PSO
	std::vector<tParticule> Essaim;		//**Ensemble de solutions 
	tPosition Best;						//**Meilleure solution depuis le d�but de l'algorithme
	
	//string NomFichier;					//### Pour MAXCUT
	
	//**Lecture des param�tres
	LePSO.Taille		= atoi(Param[1]);
	LePSO.C1			= atof(Param[2]);
	LePSO.C2			= atof(Param[3]);
	LePSO.C3			= atof(Param[4]);
	LePSO.NB_EVAL_MAX	= atoi(Param[5]);
	//NomFichier.assign(Param[6]);			//### Pour MAXCUT
	LePSO.Iter			= 0;
	LePSO.CptEval		= 0;
		
	srand((unsigned) time(NULL));			//**Precise un germe pour le generateur aleatoire
	cout.setf(ios::fixed|ios::showpoint);
	
	//**Sp�cifications du probl�me � traiter
	LeProb.Fonction = ALPINE;				//**Sp�cifie le probl�me trait�
	InitialisationIntervalleVariable(LeProb);
	
	//**Lecture du fichier de MAXCUT
	//LectureProblemeMAXCUT(NomFichier, LeProb, LePSO);			//### Pour MAXCUT
	//AfficherProblemeMAXCUT(LeProb);							//### Pour MAXCUT
	
	//**Dimension du tableaux de l'essaim selon le nombre de particules
	Best = InitialisationEssaim(Essaim, LeProb, LePSO);
	InitialisationInformatrices(Essaim, LePSO);

	//AfficherSolutions(Essaim, 0, LePSO.Taille, LePSO.Iter, LeProb);
	AfficherUneSolution(Best, LePSO.Iter, LeProb);

	//**Boucle principale du PSO
	while (LePSO.CptEval < LePSO.NB_EVAL_MAX) 	//**NE PAS ENLEVER LA CONDITION SUR LE NOMBRE D'�VALUATION
	{																									
		LePSO.Iter++;
		DeplacerEssaim(Essaim, LeProb, LePSO, Best);  
		//AfficherSolutions(Essaim, 0, LePSO.Taille, LePSO.Iter, LeProb);
		AfficherUneSolution(Best, LePSO.Iter, LeProb);
	};

	AfficherResultats (Best, LeProb, LePSO);		//**NE PAS ENLEVER
	LibererMemoireFinPgm(Essaim, LeProb, LePSO);

	system("PAUSE");
	return 0;
}

//**-----------------------------------------------------------------------
//**D�termine l'intervalle de recherche selon la fonction choisie
void InitialisationIntervalleVariable(tProblem &unProb)
{
	switch(unProb.Fonction)
	{
		case ALPINE:	unProb.Xmin = -10.0;	unProb.Xmax = 10.0;	unProb.D = 2; break;
		case BANANE:	unProb.Xmin = -10.0;	unProb.Xmax = 10.0;	unProb.D = 2; break;
		default:		unProb.Xmin = 0.0;		unProb.Xmax = 0.0;	unProb.D = 0; break; 
	}
}

//**-----------------------------------------------------------------------
//**Dimension des vecteurs de l'essaim, initialisation des particules et retourne la meilleure solution rencontr�e
tPosition InitialisationEssaim(std::vector<tParticule> &unEssaim, tProblem unProb, tPSO &unPSO)
{
	int i;
	tPosition Meilleure;

	unEssaim.resize(unPSO.Taille);		//**Le tableau utilise les indices de 0 � Taille-1.
	
	for (i=0; i<unPSO.Taille; i++)
	{
		unEssaim[i].Pos.X.resize(unProb.D);
		unEssaim[i].V.resize(unProb.D);
		InitialisationPositionEtVitesseAleatoire(unEssaim[i], unProb);
		EvaluationPosition(unEssaim[i].Pos, unProb, unPSO);
		//Initialisation de la meilleure position de la particule � sa position initiale
		unEssaim[i].BestPos.X.resize(unProb.D);
		unEssaim[i].BestPos = unEssaim[i].Pos;
		//Conservation de la meilleure solution initiale
		if (i == 0)	
			Meilleure = unEssaim[i].Pos;
		else
			if (unEssaim[i].Pos.FctObj < Meilleure.FctObj)
				Meilleure = unEssaim[i].Pos;
	}
	//Retourne la meilleure solution renontr�e jusqu'� maintenant
	return(Meilleure);
}

//**-----------------------------------------------------------------------
//**D�termine le groupement d'informatrices des particules
void InitialisationInformatrices(std::vector<tParticule> &unEssaim, tPSO &unPSO)
{
	int i;
	
	//***** � D�FINIR PAR L'�TUDIANT *****
	//M�thode BIDON: une particule a comme informatrice seulement elle-m�me
	unPSO.NbInfo = 1;		//� D�TERMINER: le nombre d'informatrices par particules
	for(i=0; i<unPSO.Taille; i++)
	{
		//Dimension du vecteur d'informatrices
		unEssaim[i].Info.resize(unPSO.NbInfo);

		//S'ajoute elle-m�me
		unEssaim[i].Info[0] = & unEssaim[i];
	}
}

//**-----------------------------------------------------------------------
//**Initialise al�atoirement le vecteur de positions et de vitesses
void InitialisationPositionEtVitesseAleatoire(tParticule &Particule, tProblem unProb)
{
	int d;

	for(d=0; d<unProb.D; d++)
	{
		Particule.Pos.X[d] = AleaDouble(unProb.Xmin, unProb.Xmax);
		Particule.V[d] = AleaDouble(unProb.Xmin, unProb.Xmax);
	}
}

//**-----------------------------------------------------------------------
//**�value la position de la particule avec la fonction choisie
void EvaluationPosition(tPosition &Pos, tProblem unProb, tPSO &unPSO)
{
	double xd, som1=0.0, som2=0.0, valeur=0.0, p=1;
	int d;

	switch(unProb.Fonction)
	{
		case ALPINE: //Alpine: Min 0 en (0,0 ... 0)
				for(d=0; d<unProb.D; d++) valeur = valeur + fabs(Pos.X[d]*(sin(Pos.X[d])+0.1));
				break;
		case BANANE: //Rosenbrock ou fonction Banane. Min 0 en (1, ... 1)
				for(d=0; d<unProb.D-1; d++)
				{
					xd = 1-Pos.X[d];
					valeur += xd*xd;
					xd = Pos.X[d]*Pos.X[d]-Pos.X[d+1];
					valeur += 100*xd*xd;
				}
				break;
		default: valeur = FLT_MAX;
	}
	Pos.FctObj = valeur;
	unPSO.CptEval ++;
}

//-----------------------------------------------------------------------
void DeplacerEssaim(std::vector<tParticule> &unEssaim, tProblem unProb, tPSO &unPSO, tPosition &Meilleure)
{
	int i;
	//D�placement de l'essaim pour chaque particule----------------------
	for(i=0; i<unPSO.Taille; i++)
	{
		DeplacerUneParticule(unEssaim[i], unProb, unPSO);
		
		//�valuation de la nouvelle position-----------------------------
		EvaluationPosition(unEssaim[i].Pos, unProb, unPSO);
		
		//Mise � jour de la meilleure position de la particule-----------
		if(unEssaim[i].Pos.FctObj <= unEssaim[i].BestPos.FctObj)
		{
			unEssaim[i].BestPos = unEssaim[i].Pos;
			//M�morisation du meilleur r�sultat atteint jusqu'ici------------
			if(unEssaim[i].BestPos.FctObj < Meilleure.FctObj)
				Meilleure = unEssaim[i].BestPos;
		}
	}
}

//-----------------------------------------------------------------------
//D�placement d'une seule particule
void DeplacerUneParticule(tParticule &Particule, tProblem unProb, tPSO &unPSO)
{
	tParticule* MeilleureInfo;
	int d;

	//Meilleure informatrice de la particule-----------------------------
	MeilleureInfo = TrouverMeilleureInformatrice(Particule, unPSO);

	//Calcul de la nouvelle vitesse--------------------------------------
	for(d=0; d<unProb.D; d++)
		Particule.V[d] =	unPSO.C1 * Particule.V[d] + 
							AleaDouble(0,unPSO.C2) * (Particule.BestPos.X[d] - Particule.Pos.X[d]) + 
							AleaDouble(0,unPSO.C3) * (MeilleureInfo->BestPos.X[d] - Particule.Pos.X[d]);

	//Mise � jour de la nouvelle position--------------------------------
	for(d=0; d<unProb.D; d++)
		Particule.Pos.X[d] += Particule.V[d];

	//Confinement d'intervalle pour la valeur des positions--------------
	for(int d=0; d<unProb.D; d++)
	{
		if(Particule.Pos.X[d] < unProb.Xmin)
		{
			Particule.Pos.X[d] = unProb.Xmin;
			Particule.V[d] = 0; //Remise � z�ro de la vitesse
		}
		if(Particule.Pos.X[d] > unProb.Xmax)
		{
			Particule.Pos.X[d] = unProb.Xmax;
			Particule.V[d] = 0; //Remise � z�ro de la vitesse
		}
	}
}

//-----------------------------------------------------------------------
//Trouve la meilleure informatrice d'une particule
tParticule* TrouverMeilleureInformatrice(tParticule &Particule, tPSO &unPSO)
{
	int Rang;
	int k;
	double Valeur;

	Rang = 0;
	Valeur = Particule.Info[0]->BestPos.FctObj;
	for(k=1; k<unPSO.NbInfo; k++)
	{
		if(Particule.Info[k]->BestPos.FctObj < Valeur)
		{
			Valeur = Particule.Info[k]->BestPos.FctObj;
			Rang = k;
		}
	}
	return Particule.Info[Rang];
}

//**-----------------------------------------------------------------------
//**Renvoie un double al�atoire entre a et b
double AleaDouble(double a, double b)
{
	double resultat = (double)rand()/RAND_MAX;
	return (resultat*(b-a) + a);
}

//**-----------------------------------------------------------------------
//DESCRIPTION: Fonction qui affiche le d�tail des solutions (de Debut jusqu'a Fin-1) de l'essaim
void AfficherSolutions(std::vector<tParticule> unEssaim, int Debut, int Fin, int Iter, tProblem unProb)
{
	int i,j;
	cout << "AFFICHAGE ESSAIM (Iteration #" << Iter << ")" << endl;
	for (i=Debut; i<Fin; i++)
	{
		cout << "SOL: " << i << "\t";
		for(j=0; j<unProb.D; j++)
		{
			cout << setprecision(6) << setw(10) << unEssaim[i].Pos.X[j];
		}
		cout <<  "\tFctObj: " << setw(10) << unEssaim[i].Pos.FctObj;
		cout << endl;
		cout << "BEST POS: " << i << "\t";
		for(j=0; j<unProb.D; j++)
		{
			cout << setprecision(6) << setw(10) << unEssaim[i].BestPos.X[j];
		}
		cout <<  "\tFctObj: " << setw(10) << unEssaim[i].BestPos.FctObj;
		cout << endl;
	}
}

//**-----------------------------------------------------------------------
//Fonction qui affiche le d�tail d'une solution
void AfficherUneSolution(tPosition Pos, int Iter, tProblem unProb)
{
	int j;
	cout << "ITER #" << Iter << " \t";

	cout << "SOL: ";
	for(j=0; j<unProb.D; j++)
	{
		cout << setprecision(6) << setw(10) << Pos.X[j];
	}
	cout <<  "\tFctObj: " << setw(10) << Pos.FctObj;
	cout << endl;
}

//**-----------------------------------------------------------------------
//DESCRIPTION: Fonction affichant les r�sultats du PSO
void AfficherResultats (tPosition uneBest, tProblem unProb, tPSO unPSO)
{
	cout.setf(ios::fixed|ios::showpoint);
	cout << endl << endl;
	cout << "===============================================================================" <<endl;
	cout << "                                    RESULTATS FINAUX" << endl;
	cout << "===============================================================================" <<endl;
	cout << "FONCTION: " ;
	switch (unProb.Fonction)
	{
		case ALPINE: cout << "ALPINE"; break;
		case BANANE: cout << "BANANE"; break;
		default: cout << " a definir...";
	}
	cout << endl; 
	cout << "PARAMETRES" << endl;
	cout << "Taille : " << unPSO.Taille << "\tC1: " << setprecision(6) << unPSO.C1 
		 << "\tC2: " << setprecision(6) << unPSO.C2 << "\tC3: " << setprecision(6) << unPSO.C3 
		 << "\tNbInfo: " << setprecision(2) << unPSO.NbInfo << endl; 
	cout << "===============================================================================" << endl; 
	cout << "Nombre d iterations effectuees : " << unPSO.Iter << endl;
	cout << "Nombre Total d'evaluations : " << unPSO.CptEval << "/" << unPSO.NB_EVAL_MAX << endl;
	cout << "Meilleure solution trouvee: "  << setprecision(6) << uneBest.FctObj << endl; 
	cout << endl << "===============================================================================" << endl;
}

//**-----------------------------------------------------------------------
//**Liberation de la m�moire allou�e dynamiquement
void LibererMemoireFinPgm (std::vector<tParticule> &unEssaim, tProblem unProb, tPSO unPSO)
{
	int i;

	for (i=0; i<unPSO.Taille; i++)
	{
		unEssaim[i].Pos.X.clear();
		unEssaim[i].BestPos.X.clear();
		unEssaim[i].V.clear();
		unEssaim[i].Info.clear();
	}
	unEssaim.clear();
}

//**-----------------------------------------------------------------------
//**Lecture du Fichier probleme MAXSAT et initialiation de la structure Problem
//**NB: Lors de la lecture, le num�ro des noeuds est transform� (0 � NbNoeud-1 AU LIEU de 1 � NbNoeud)
/* ### POUR MAXCUT
void LectureProblemeMAXCUT(std::string FileName, tProblem & unProb, tPSO &unPSO)
{
	ifstream	Fichier;
	int x, i, j;

	unProb.Nom = FileName;
	Fichier.open(unProb.Nom.c_str(), ios::in) ;
	if (Fichier.fail())
		cout << "Erreur � l'ouverture" << endl ;
	else
	{
		Fichier >> unProb.NbNoeud >> unProb.NbArc; // Lecture du nombre de noeuds (dimension) et du nombre d'arcs
		unProb.D = unProb.NbNoeud;

		//Allocation dynamique en fonction de la taille du probleme	
		unProb.Arc.resize(unProb.NbArc);	
		
		for (x=0; x<unProb.NbArc; x++)
		{
			Fichier >> i >> j >> unProb.Arc[x].Poids;
			unProb.Arc[x].Ni = i-1; //Transformation du num�ro des arcs : 0 � NbNoeud-1  (au lieu de 1 � NbNoeud)
			unProb.Arc[x].Nj = j-1;
		}
		
		Fichier.close() ;
		if (Fichier.fail())
			cout << "Erreur � la fermeture" << endl ;
	}
}
//**-----------------------------------------------------------------------
//Fonction d'affichage � l'�cran permettant de voir si les donn�es du fichier probl�me ont �t� lues correctement
void AfficherProblemeMAXCUT (tProblem unProb)
{
	int i;
	cout << "*******************************************************" << endl;
	cout << "NOM FICHIER:\t"		<<  unProb.Nom	<< endl;
	cout << "NOMBRE NOEUDS:\t"		<<  unProb.NbNoeud	<< endl;
	cout << "NOMBRE ARCS:\t"		<<  unProb.NbArc << endl;
	cout << endl << "ARC: ";
	for (i=0; i<unProb.NbArc; i++)
	{
		cout << endl << "#" << i << ": ";
		cout << "\t" << unProb.Arc[i].Ni << "\t" << unProb.Arc[i].Nj;
		cout << "\tPoids= " << setw(2) << unProb.Arc[i].Poids;
	}
	cout << endl << "*******************************************************";
	cout << endl << endl;
}
*/
