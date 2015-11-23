#include "Entete.h"

//%%%%%%%%%%%%%%%%%%%%%%%%% IMPORTANT: %%%%%%%%%%%%%%%%%%%%%%%%% 
//Le fichier de probleme (.txt) doit se trouver dans le répertoire courant du projet 
//pour une exécution à l'aide du compilateur. Indiquer les arguments du programme dans les 
//propriétés du projet - débogage - arguements.
//Sinon, utiliser le répertoire execution.

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
void LectureProblemeMAXCUT(std::string FileName, tProblem & unProb, tPSO &unPSO);
void AfficherProblemeMAXCUT (tProblem unProb);

// RAJOUT
void displayEssaim(std::vector<tParticule> &unEssaim) {
	std::cout << '[';
	for ( int i = 0 ; (unsigned)i < unEssaim.size(); i++ ) {
		std::cout << '(';
		std::cout << unEssaim[i].Pos.X[0];
		std::cout << ',';
		std::cout << unEssaim[i].Pos.X[1];
		std::cout << ',';
		std::cout << unEssaim[i].Pos.FctObj;
		std::cout << ')';
		if ( (unsigned)i < unEssaim.size()-1 )
			std::cout << ',';
	}
	std::cout << ']';
}

/*void LectureProblemeMAXCUT(std::string FileName, tProblem & unProb, tPSO &unPSO){
	std::fstream file;
	file.open(FileName.c_str(), std::fstream::in);
	int ni, nj, poids;
	tArc tmp_arc;
	if ( file.is_open() ){
		file >> unProb.NbNoeud >> unProb.NbArc;
		for ( int i = 0 ; i < unProb.NbArc ; i++ ){
			file >> ni >> nj >> poids;
			tmp_arc.Ni = ni;
			tmp_arc.Nj = nj;
			tmp_arc.Poids = poids;
			unProb.Arc.push_back(tmp_arc);
		}
		file.close();
	} else {
		std::cout << "ERREUR lecture fichier" << std::endl;
	}
}

void AfficherProblemeMAXCUT (tProblem unProb){
	std::cout << "=====" << std::endl;
	for ( int i = 0 ; i < unProb.NbArc ; i++ ){
		std::cout << unProb.Arc[i].Ni << " " << unProb.Arc[i].Nj << " " << unProb.Arc[i].Poids << std::endl;
	}
	std::cout << "=====" << std::endl;
}*/

/* ### POUR MAXCUT

*/
//******************************************************************************************
// Fonction main
//*****************************************************************************************
int main(int NbParam, char *Param[])
{
	tProblem LeProb;					//**Définition de l'instance de problème
	tPSO LePSO;							//**Définition des paramètres du PSO
	std::vector<tParticule> Essaim;		//**Ensemble de solutions 
	tPosition Best;						//**Meilleure solution depuis le début de l'algorithme
	
	string NomFichier;					//### Pour MAXCUT
	
	//**Lecture des paramètres
	LePSO.Taille		= atoi(Param[1]);
	LePSO.C1			= atof(Param[2]);
	LePSO.C2			= atof(Param[3]);
	LePSO.C3			= atof(Param[4]);
	LePSO.NB_EVAL_MAX	= atoi(Param[5]);
	NomFichier.assign(Param[6]);			//### Pour MAXCUT
	LePSO.Iter			= 0;
	LePSO.CptEval		= 0;
		
	srand((unsigned) time(NULL));			//**Precise un germe pour le generateur aleatoire
	cout.setf(ios::fixed|ios::showpoint);
	
	//**Spécifications du problème à traiter
	LeProb.Fonction = MAXCUT;				//**Spécifie le problème traité
	InitialisationIntervalleVariable(LeProb);
	
	//**Lecture du fichier de MAXCUT
	LectureProblemeMAXCUT(NomFichier, LeProb, LePSO);			//### Pour MAXCUT
	AfficherProblemeMAXCUT(LeProb);							//### Pour MAXCUT
	
	//**Dimension du tableaux de l'essaim selon le nombre de particules
	Best = InitialisationEssaim(Essaim, LeProb, LePSO);
	InitialisationInformatrices(Essaim, LePSO);

	//AfficherSolutions(Essaim, 0, LePSO.Taille, LePSO.Iter, LeProb);
	AfficherUneSolution(Best, LePSO.Iter, LeProb);

	//**Boucle principale du PSO
	while (LePSO.CptEval < LePSO.NB_EVAL_MAX) 	//**NE PAS ENLEVER LA CONDITION SUR LE NOMBRE D'ÉVALUATION
	{																									
		LePSO.Iter++;
		DeplacerEssaim(Essaim, LeProb, LePSO, Best);  
		//~ AfficherSolutions(Essaim, 0, LePSO.Taille, LePSO.Iter, LeProb);
		AfficherUneSolution(Best, LePSO.Iter, LeProb);
	};
	//~ displayEssaim(Essaim);
	AfficherResultats (Best, LeProb, LePSO);		//**NE PAS ENLEVER
	LibererMemoireFinPgm(Essaim, LeProb, LePSO);

	system("PAUSE");
	return 0;
}

//**-----------------------------------------------------------------------
//**Détermine l'intervalle de recherche selon la fonction choisie
void InitialisationIntervalleVariable(tProblem &unProb)
{
	switch(unProb.Fonction)
	{
		case ALPINE:	unProb.Xmin = -10.0;	unProb.Xmax = 10.0;		unProb.D = 2; break;
		case BANANE:	unProb.Xmin = -10.0;	unProb.Xmax = 10.0;		unProb.D = 2; break;
		case EGGHOLDER:	unProb.Xmin = -512.0;	unProb.Xmax = 512.0;	unProb.D = 2; break;
		case MAXCUT:	unProb.Xmin = 0;		unProb.Xmax = 1;		unProb.D = unProb.NbNoeud; break;
		default:		unProb.Xmin = 0.0;		unProb.Xmax = 0.0;		unProb.D = 0; break; 
	}
}

//**-----------------------------------------------------------------------
//**Dimension des vecteurs de l'essaim, initialisation des particules et retourne la meilleure solution rencontrée
tPosition InitialisationEssaim(std::vector<tParticule> &unEssaim, tProblem unProb, tPSO &unPSO)
{
	int i;
	tPosition Meilleure;

	unEssaim.resize(unPSO.Taille);		//**Le tableau utilise les indices de 0 à Taille-1.
	
	for (i=0; i<unPSO.Taille; i++)
	{
		unEssaim[i].Pos.X.resize(unProb.D);
		unEssaim[i].V.resize(unProb.D);
		InitialisationPositionEtVitesseAleatoire(unEssaim[i], unProb);
		EvaluationPosition(unEssaim[i].Pos, unProb, unPSO);
		//Initialisation de la meilleure position de la particule à sa position initiale
		unEssaim[i].BestPos.X.resize(unProb.D);
		unEssaim[i].BestPos = unEssaim[i].Pos;
		//Conservation de la meilleure solution initiale
		if (i == 0)	
			Meilleure = unEssaim[i].Pos;
		else
			if (unEssaim[i].Pos.FctObj > Meilleure.FctObj)
				Meilleure = unEssaim[i].Pos;
	}
	//Retourne la meilleure solution renontrée jusqu'à maintenant
	return(Meilleure);
}

//**-----------------------------------------------------------------------
//**Détermine le groupement d'informatrices des particules
void InitialisationInformatrices(std::vector<tParticule> &unEssaim, tPSO &unPSO)
{
	int i,j;
	
	unPSO.NbInfo = 10;		//À DÉTERMINER: le nombre d'informatrices par particules
	for(i=0; i<unPSO.Taille; i++)
	{
		//Dimension du vecteur d'informatrices
		unEssaim[i].Info.resize(unPSO.NbInfo);

		// RAJOUT
		//construction du réseau d'informatrices
		for (j=0; j<unPSO.NbInfo; j++){
			//~ unEssaim[i].Info[j] = &unEssaim[(i+j+1)%(unEssaim.size())];
			unEssaim[i].Info[j] = &unEssaim[(int)AleaDouble(0.,(double)unEssaim.size())];
		}
	}
}

//**-----------------------------------------------------------------------
//**Initialise aléatoirement le vecteur de positions et de vitesses
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
//**Évalue la position de la particule avec la fonction choisie
void EvaluationPosition(tPosition &Pos, tProblem unProb, tPSO &unPSO)
{
	double xd, som1=0.0, som2=0.0, valeur=0.0, p=1;
	
	// ADDING
	double x1, x2;
	
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
		case EGGHOLDER: // Calcul de la fonction eggholder
				x1 = Pos.X[0];
				x2 = Pos.X[1];
				valeur = -(x2+47.)*sin(sqrt( fabs(x2+x1/2.+47) ));
				valeur -= x1*sin( sqrt(fabs( x1 - (x2+47) )) );
			break;
		case MAXCUT: // calcul de la fonction maxcut
				valeur = 0;
				//~ std::cout << "===" << std::endl;
				for ( int u = 0 ; u < unProb.NbNoeud ; u++ ) {
					for ( int e = 0 ; e < unProb.NbArc ; e++ ) {
						if ( unProb.Arc[e].Ni == u ) {
							if ( Pos.X[unProb.Arc[e].Ni] == 1-Pos.X[unProb.Arc[e].Nj] ) {
								valeur += unProb.Arc[e].Poids;
								//~ std::cout << unProb.Arc[e].Ni << " " << unProb.Arc[e].Nj << std::endl;	
							}
						}
					}
				}
				//~ std::cout << "===" << std::endl;
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
	//Déplacement de l'essaim pour chaque particule----------------------
	for(i=0; i<unPSO.Taille; i++)
	{
		DeplacerUneParticule(unEssaim[i], unProb, unPSO);
		
		//Évaluation de la nouvelle position-----------------------------
		EvaluationPosition(unEssaim[i].Pos, unProb, unPSO);
		
		//Mise à jour de la meilleure position de la particule-----------
		if(unEssaim[i].Pos.FctObj >= unEssaim[i].BestPos.FctObj)
		{
			unEssaim[i].BestPos = unEssaim[i].Pos;
			//Mémorisation du meilleur résultat atteint jusqu'ici------------
			if(unEssaim[i].BestPos.FctObj > Meilleure.FctObj)
				Meilleure = unEssaim[i].BestPos;
		}
	}
}

// RAJOUT
double s(double v) {
	return 1./(1.+exp(-v));
}

//-----------------------------------------------------------------------
//Déplacement d'une seule particule
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


	//Mise à jour de la nouvelle position--------------------------------
	if ( unProb.Fonction == MAXCUT ) {
		
		// maj de la vitesse si trop élevée
		if ( Particule.V[d] < -4. ) {
			Particule.V[d] = -4.; 
		}
		if ( Particule.V[d] > 4. ) {
			Particule.V[d] = 4.;
		}
		
		for(d=0; d<unProb.D; d++){
			if ( AleaDouble(0., 1.) < s(Particule.V[d]) ) {
				Particule.Pos.X[d] = 1.;
			} else {
				Particule.Pos.X[d] = 0.;
			}
		}
			
	} else {
		for(d=0; d<unProb.D; d++)
			Particule.Pos.X[d] += Particule.V[d];
		//Confinement d'intervalle pour la valeur des positions--------------
		for(int d=0; d<unProb.D; d++)
		{
			if(Particule.Pos.X[d] < unProb.Xmin)
			{
				Particule.Pos.X[d] = unProb.Xmin;
				Particule.V[d] = 0; //Remise à zéro de la vitesse
			}
			if(Particule.Pos.X[d] > unProb.Xmax)
			{
				Particule.Pos.X[d] = unProb.Xmax;
				Particule.V[d] = 0; //Remise à zéro de la vitesse
			}
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
		if(Particule.Info[k]->BestPos.FctObj > Valeur)
		{
			Valeur = Particule.Info[k]->BestPos.FctObj;
			Rang = k;
		}
	}
	return Particule.Info[Rang];
}

//**-----------------------------------------------------------------------
//**Renvoie un double aléatoire entre a et b
double AleaDouble(double a, double b)
{
	double resultat = (double)rand()/RAND_MAX;
	return (resultat*(b-a) + a);
}

//**-----------------------------------------------------------------------
//DESCRIPTION: Fonction qui affiche le détail des solutions (de Debut jusqu'a Fin-1) de l'essaim
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
//Fonction qui affiche le détail d'une solution
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
//DESCRIPTION: Fonction affichant les résultats du PSO
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
		case MAXCUT: cout << "MAXCUT"; break;
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
//**Liberation de la mémoire allouée dynamiquement
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
//**NB: Lors de la lecture, le numéro des noeuds est transformé (0 à NbNoeud-1 AU LIEU de 1 à NbNoeud)
/* ### POUR MAXCUT */
void LectureProblemeMAXCUT(std::string FileName, tProblem & unProb, tPSO &unPSO)
{
	ifstream	Fichier;
	int x, i, j;

	unProb.Nom = FileName;
	Fichier.open(unProb.Nom.c_str(), ios::in) ;
	if (Fichier.fail())
		cout << "Erreur a l'ouverture" << endl ;
	else
	{
		Fichier >> unProb.NbNoeud >> unProb.NbArc; // Lecture du nombre de noeuds (dimension) et du nombre d'arcs
		unProb.D = unProb.NbNoeud;

		//Allocation dynamique en fonction de la taille du probleme	
		unProb.Arc.resize(unProb.NbArc);	
		
		for (x=0; x<unProb.NbArc; x++)
		{
			Fichier >> i >> j >> unProb.Arc[x].Poids;
			unProb.Arc[x].Ni = i-1; //Transformation du numéro des arcs : 0 à NbNoeud-1  (au lieu de 1 à NbNoeud)
			unProb.Arc[x].Nj = j-1;
		}
		
		Fichier.close() ;
		if (Fichier.fail())
			cout << "Erreur a la fermeture" << endl ;
	}
}

//**-----------------------------------------------------------------------
//Fonction d'affichage à l'écran permettant de voir si les données du fichier problème ont été lues correctement
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

