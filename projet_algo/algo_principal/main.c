#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <time.h>
#include "fonctions.h"

bool convergence;
int taille_motif; //longueur motif/masque IE. l
int d; // nombre maximal de substitutions autorisées
int nb_fenetres; //nombre de fenêtres dans le masque IE. k
int nb_masques; //nombre d'essais réalisés
/*************************************************************************************************************
 * * *                                    MAIN                                                           * * *
 *************************************************************************************************************/

int main(int argc, char *argv[])
{
	srand(time(NULL));
	int i, j;
	/********************************************
	*    Définition des variables principales   *
	*********************************************/
	int* masque;
	// memset(masque, 0, sizeof masque); //initialisation du masque
	double** matrice_PSSM;
	double** matrice_PSSM_nouv;
	char* motif_consensus =  (char*) calloc(taille_motif+1, sizeof(char));
	TPtr_Cellkmer tableau_kmer_QS[];

	/******************************
	*   Dictionnaire Sequences    *
	*******************************/
	TPtr_info_dictionnaire_sequences tete_info_dict_seq = malloc(sizeof(TInfo_dictionnaire_sequences)); //??? quelque chose avant malloc
	TPtr_dictionnaire_sequences tete_dictionnaire_sequences = malloc(sizeof(TDictionnaire_Sequences));
	tete_info_dict_seq->tete_dictionnaire_seq = tete_dictionnaire_sequences; //tete dictionnaire des séquences
	TPtr_dictionnaire_sequences p_dictionnaire_sequence = tete_dict_seq;

	/*******************************
	*    Dictionnaire KMERs       *
	********************************/
	TPtr_info_dictionnaire_kmer tete_info_dict_kmer = malloc(sizeof(TInfo_dictionnaire_kmer));
	tete_info_dict_kmer->nb_kmer = 0;
	tete_info_dict_kmer->tete_liste_kmer = NULL;

	TPtr_Cellkmer p_kmer = tete_info_dict_kmer->tete_liste_kmer; //on créer des pointeurs temporaires pour parcourir les listes
	TPtr_Cellkmer p_prec_kmer = p_kmer;

	/*************************************************************************************************************
	 * * *                                    DEBUT DE L'ALGORITHME                                          * * *
	 *************************************************************************************************************/
	importer_parametres( argv[1], &taille_motif, &d, &nb_fenetres, &nb_masques);
	importer_sequences_fasta( argv[2], &tete_info_dict_seq, &tete_dictionnaire_sequences );
  // afficher_sequences(&tete_info_dict_seq, &tete_dictionnaire_sequences );

	for ( i=0; i <= nb_masques; i++ )
	{
		masque = generation_masque(masque);
		parcours_masque_sur_seq( masque, tete_dictionnaire_sequences, tete_info_dict_kmer);
		// affichage_dictionnaire_kmer(tete_info_dict_kmer);

		p_dictionnaire_sequence = tete_dict_seq;
		p_kmer = tete_info_dict_kmer->tete_liste_kmer;
		while( p_kmer != NULL ) //pour chaque kmer entrée du dictionnaire de kmer
		{
			if ( p_kmer->nb_sequence == p_dictionnaire_sequence->nb_sequences )   // chaque kmer doit être sufisemment représenté ie. au moins présent une fois dans chaque séquence
			{
				convergence = FALSE;
				creation_PSSM(&matrice_PSSM);
				calcul_PSSM(p_kmer, &matrice_PSSM);
				// afficher_PSSM(matrice_PSSM);

				do //repeter l'amélioration de la PSSM jusqu'à convergence
				{
					creation_PSSM(&nouv_matrice_PSSM);
					recherche_motifs_maximisant_scores(p_kmer, matrice_PSSM, tete_info_dict_kmer); // revient à définir l'ensemble T
					calcul_PSSM_amelioree(p_kmer, &nouv_matrice_PSSM);

					if ( calcul_distance_PSSMs(&matrice_PSSM, &matrice_PSSM_nouv) > 0.8 )
					{
						copie_PSSM(&matrice_PSSM, &matrice_PSSM_nouv);
						liberation_PSSM(&matrice_PSSM_nouv);
					}else{  convergence = TRUE; }
				}while( convergence == FALSE );

				p_kmer->PSSM_consensus = matrice_PSSM;
				p_kmer->motif_consensus = identification_motif_consensus(matrice_PSSM);
				raffiner_version1(p_kmer);
				raffiner_version2(p_kmer);
			}
			p_prec_kmer = p_kmer;
			p_kmer = p_kmer->suiv_kmer;
		}
		creation_vecteur_kmer_pour_QuickSort(&tableau_kmer_QS);
		quick_sort_ST(&tableau_kmer_QS);
		
		generation_fichier_resultats(tableau_kmer_QS, i);
	}
	liberation_dictionnaire_kmer(tete_info_dict_kmer);
	liberation_dictionnaire_sequence(&tete_info_dict_seq, &tete_dictionnaire_sequences);
	free(motif_consensus)
	free(tete_info_dict_kmer);
	return(0);
}
