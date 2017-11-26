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
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
int main()
{
	//définition des variables utilitaires
	srand(time(NULL));
	int i,j;

	int masque[longueur_masque];
	memset(masque, 0, sizeof masque);

	TInfo_ensemble_sequences structure_info_seq; //instanciation de la structure contenants les informations des sequences
	TPtr_info_ensemble_sequences ptr_info = malloc(sizeof(TInfo_ensemble_sequences));
	TPtr_ensemble_sequences ptr_ensemble = malloc(sizeof(TEnsemble_Sequences));
	ptr_info->tete_ensemble_seq = ptr_ensemble; //tete dictionnaire des séquences

	TPtr_ensemble_sequences element_generation_sequence=malloc(sizeof(TEnsemble_Sequences)); // On cr�� le premier �l�ment de la structure
  TPtr_ensemble_sequences tete_liste_pour_parcours_masque= element_generation_sequence;
  TPtr_ensemble_sequences p_generation_seq= element_generation_sequence;
  TPtr_ensemble_sequences tete_liste_pour_recup_motif= element_generation_sequence;
  TPtr_ensemble_sequences tete_liste_pour_insertion_motif= element_generation_sequence;
  TPtr_ensemble_sequences element_generation_sequence_next= NULL;
  ptr_liste_motif element_motif= malloc(sizeof(liste_chaine_motif));
  ptr_liste_motif tete_liste_motif= element_motif;


	TPtr_Cellkmer element_kmer=malloc(sizeof(TCellkmer));
	TPtr_Cellkmer tete_liste_kmer= element_kmer;
	TPtr_Cellkmer tete_liste_kmer2= element_kmer;
	TPtr_Cellkmer tete_liste_kmer3= element_kmer;
	TPtr_Cellkmer tete_liste_kmer4= element_kmer;
	element_kmer->suiv_kmer=NULL;
	TPtr_CellSequence element_sequence= malloc(sizeof(TCellSequence));
	TPtr_CellSequence tete_liste_sequence= element_sequence;
	TPtr_CellSequence tete_liste_sequence2= element_sequence;
	TPtr_CellSequence tete_liste_sequence3= element_sequence;
	element_sequence->suiv_sequence=NULL;
	TPtr_CellPos element_pos= malloc(sizeof(TCellPos));
	TPtr_CellPos tete_liste_pos= element_pos;
	TPtr_CellPos tete_liste_pos2= element_pos;
	TPtr_CellPos tete_liste_pos3= element_pos;
	element_pos->suiv_pos=NULL;
	TPtr_Cellkmer_selectionne element_liste_kmer_selectionne= malloc(sizeof(TCellkmer_selectionne));
  TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne= element_liste_kmer_selectionne;
  TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne2= element_liste_kmer_selectionne;
  TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne_pour_calcul= element_liste_kmer_selectionne;
  TPtr_Cell_Motif_PSSM element_liste_motif_PSSM= malloc(sizeof(TCell_Motif_PSSM));
  TPtr_Cell_Motif_PSSM tete_liste_motif_PSSM= element_liste_motif_PSSM;
  TPtr_Cell_Motif_PSSM tete_liste_motif_PSSM2= element_liste_motif_PSSM;
  TPtr_Cell_Motif_PSSM tete_liste_motif_PSSM_pour_calcul= element_liste_motif_PSSM;
  element_liste_motif_PSSM->suiv_motif=NULL;

  TPtr_Mot_Ameliorer_PSSM element_mot= (TMot_Ameliorer_PSSM*)malloc(sizeof(TMot_Ameliorer_PSSM));
  TPtr_Mot_Ameliorer_PSSM p_mot= element_mot;

  //MATRICE PSSM ET MOTIF CONSENSUS:
  double** matrice_PSSM[4][6];

//------------------------------------------------------------------------------------------------------------
	//début de l'algorithme

	importer_parametres(&longueur_masque, &d, &nb_fenetre, &nb_masques);
	printf("longueur_masque= %d, ",longueur_masque);printf("d= %d, ",d);
	printf("nb_fenetre= %d, ",nb_fenetre);printf("nb_masques= %d\n",nb_masques);

	importer_sequences_fasta(&ptr_info, &ptr_ensemble );
  afficher_sequences(&ptr_info, &ptr_ensemble );
	printf("done\n");

	for ( i=0; i<= nb_masques ; i++ )
	{
		generation_masque(longueur_masque, &masque, nb_fenetre);
		printf("done generation_masque \n");
		// affiche masques
		// printf("essais n°:%d\n", i);
		// for(j=0; j<longueur_masque; j++){printf("%d\n", masque[j]);}

		// creation_dictionnaire();
		parcours_masque(longueur_masque, &masque, nb_fenetre, ptr_info->nb_sequences, &tete_liste_pour_parcours_masque, &tete_liste_kmer, &tete_liste_sequence, &tete_liste_pos);
		printf("done parcours_masque \n");
		kmer_present_dans_chaque_sequence(ptr_info->nb_sequences, &tete_liste_kmer2, &tete_liste_sequence2, &tete_liste_pos2, &tete_liste_pour_recup_motif, &tete_liste_kmer_selectionne, &tete_liste_motif_PSSM);
		printf("done kmer_present_dans_chaque_sequence\n");
		affichage_dictionnaire_kmer(&tete_liste_kmer3, &tete_liste_sequence3, &tete_liste_pos3);
		printf("done affichage_dictionnaire_kmer\n");
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
