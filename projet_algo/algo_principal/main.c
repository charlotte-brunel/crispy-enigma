#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "fonctions.h"

bool convergence;
//definition des variables principales
int l; //longueur motif
int d; // nombre maximal de substitutions autorisées
int k; //nombre de fenêtres dans le masque
int nb_masques;

int main()
{
	//définition des variables utilitaires
	int i = 0;
	int masque[l];

	TInfo_ensemble_sequences info_seq; //instanciation de la structure contenants les informations des sequences
	TInfo_ensemble_sequences* ptr_info;
	TEnsemble_Sequences* ptr_ensemble;

	ptr_info = malloc(sizeof(TInfo_ensemble_sequences));
	ptr_ensemble = malloc(sizeof(TEnsemble_Sequences));
	ptr_ensemble = ptr_info->tete_ensemble_seq;

	//début de l'algorithme

	importer_parametres(&l,&d,&k,&nb_masques);
	printf("l= %d",l);printf("d= %d",d);printf("k= %d",k);printf("nb_masques= %d \n ",nb_masques);

	importer_sequences_fasta(&ptr_info, &ptr_ensemble );
	afficher_sequences(&ptr_info, &ptr_ensemble );

	// for ( i; i <= nb_masques; i++)
	// {
	// 	generation_masque(l, &masque[l], k);
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
	// }
	//desallocation mémoire
	// free(ptr_info); free(ptr_ensemble);

	return(0);
}
