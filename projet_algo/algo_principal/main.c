#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <time.h>
#include "fonctions.h"

bool convergence;
//definition des variables principales
int longueur_masque; //longueur motif/masque IE. l
int d; // nombre maximal de substitutions autorisées
int nb_fenetre; //nombre de fenêtres dans le masque IE. k
int nb_masques;
// ------------------------------------------------------------------------------------------------
int main()
{
	//définition des variables utilitaires
	srand(time(NULL));
	int i,j;

	TInfo_ensemble_sequences info_seq; //instanciation de la structure contenants les informations des sequences
	TInfo_ensemble_sequences* ptr_info;
	TEnsemble_Sequences* ptr_ensemble;

	ptr_info = malloc(sizeof(TInfo_ensemble_sequences));
	ptr_ensemble = malloc(sizeof(TEnsemble_Sequences));
	ptr_ensemble = ptr_info->tete_ensemble_seq;

	//début de l'algorithme

	importer_parametres(&longueur_masque, &d, &nb_fenetre, &nb_masques);
	printf("longueur_masque= %d, ",longueur_masque);printf("d= %d, ",d);
	printf("nb_fenetre= %d, ",nb_fenetre);printf("nb_masques= %d\n",nb_masques);
	int masque[longueur_masque];
	memset(masque, 0, sizeof masque);

	importer_sequences_fasta(&ptr_info, &ptr_ensemble );
	// afficher_sequences(&ptr_info, &ptr_ensemble );

	for ( i=0; i<= nb_masques ; i++ )
	{
		generation_masque(longueur_masque, &masque, nb_fenetre);
		// affiche masques
		// printf("essais n°:%d\n", i);
		// for(j=0; j<longueur_masque; j++){printf("%d\n", masque[j]);}

		// creation_dictionnaire();
		//
		// for () //pour chaque k-mere suffisemment représenté
		// {
		// 	convergence = FALSE;
		// 	calculer_PSSM_ref();
		// 	do
		// 	{
		// 		ameliorer_PSSM();
		// 	}
		// 	while( convergence = FALSE )
		//
		// 	raffiner_PSSM();
		// }

	}
	return(0);
}
