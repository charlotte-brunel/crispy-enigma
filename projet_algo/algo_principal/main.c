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
	int position= 0;
	int cpt_mot;
  int n_sequence;
  char mot[longueur_masque+1];
	int masque[longueur_masque];
	memset(masque, 0, sizeof masque);

	TInfo_ensemble_sequences structure_info_seq; //instanciation de la structure contenants les informations des sequences
	TPtr_info_ensemble_sequences ptr_info = malloc(sizeof(TInfo_ensemble_sequences));
	TPtr_ensemble_sequences ptr_ensemble = malloc(sizeof(TEnsemble_Sequences));
	ptr_info->tete_ensemble_seq = ptr_ensemble; //tete dictionnaire des séquences

	TPtr_ensemble_sequences element_dictionnaire_sequence = malloc(sizeof(TEnsemble_Sequences)); // On cr�� le premier �l�ment de la structure
  TPtr_ensemble_sequences tete_liste_pour_parcours_masque = element_dictionnaire_sequence;
  TPtr_ensemble_sequences p_dictionnaire_seq = element_dictionnaire_sequence;
  TPtr_ensemble_sequences tete_liste_pour_recup_motif = element_dictionnaire_sequence;
  TPtr_ensemble_sequences tete_liste_pour_insertion_motif = element_dictionnaire_sequence;
  TPtr_ensemble_sequences element_dictionnaire_sequence_next = NULL;
  ptr_liste_motif element_motif = malloc(sizeof(liste_chaine_motif));
  ptr_liste_motif tete_liste_motif = element_motif;


	TPtr_Cellkmer element_kmer = malloc(sizeof(TCellkmer));
	TPtr_Cellkmer tete_liste_kmer = element_kmer;
	TPtr_Cellkmer tete_liste_kmer2 = element_kmer;
	TPtr_Cellkmer tete_liste_kmer3 = element_kmer;
	TPtr_Cellkmer tete_liste_kmer4 = element_kmer;
	element_kmer->suiv_kmer = NULL;
	TPtr_CellSequence element_sequence = malloc(sizeof(TCellSequence));
	TPtr_CellSequence tete_liste_sequence = element_sequence;
	TPtr_CellSequence tete_liste_sequence2 = element_sequence;
	TPtr_CellSequence tete_liste_sequence3 = element_sequence;
	element_sequence->suiv_sequence = NULL;
	TPtr_CellPos element_pos = malloc(sizeof(TCellPos));
	TPtr_CellPos tete_liste_pos = element_pos;
	TPtr_CellPos tete_liste_pos2 = element_pos;
	TPtr_CellPos tete_liste_pos3 = element_pos;
	element_pos->suiv_pos = NULL;
	TPtr_Cellkmer_selectionne element_liste_kmer_selectionne = malloc(sizeof(TCellkmer_selectionne));
  TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne = element_liste_kmer_selectionne;
  TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne2 = element_liste_kmer_selectionne;
  TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne_pour_calcul = element_liste_kmer_selectionne;
  TPtr_Cell_Motif_PSSM element_liste_motif_PSSM = malloc(sizeof(TCell_Motif_PSSM));
  TPtr_Cell_Motif_PSSM tete_liste_motif_PSSM = element_liste_motif_PSSM;
  TPtr_Cell_Motif_PSSM tete_liste_motif_PSSM2 = element_liste_motif_PSSM;
  TPtr_Cell_Motif_PSSM tete_liste_motif_PSSM_pour_calcul = element_liste_motif_PSSM;
  element_liste_motif_PSSM->suiv_motif = NULL;

  TPtr_Mot_Ameliorer_PSSM element_mot = (TMot_Ameliorer_PSSM*)malloc(sizeof(TMot_Ameliorer_PSSM));
  TPtr_Mot_Ameliorer_PSSM p_mot = element_mot;

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

	for ( i = 0; i <= nb_masques; i++ )
	{
		printf("\nessais %d\n", i);
		generation_masque(longueur_masque, &masque, nb_fenetre);
		printf("generation_masque \n");
		// affiche masques
		// printf("essais n°:%d\n", i);
		// for(j=0; j<longueur_masque; j++){printf("%d\n", masque[j]);}

		// creation_dictionnaire();
		parcours_masque(longueur_masque, &masque, nb_fenetre, ptr_info->nb_sequences, &tete_liste_pour_parcours_masque, &tete_liste_kmer, &tete_liste_sequence, &tete_liste_pos);
		printf("parcours_masque \n");
		kmer_present_dans_chaque_sequence(ptr_info->nb_sequences, &tete_liste_kmer2, &tete_liste_sequence2, &tete_liste_pos2, &tete_liste_pour_recup_motif, &tete_liste_kmer_selectionne, &tete_liste_motif_PSSM);
		printf("kmer_present_dans_chaque_sequence\n");
		affichage_dictionnaire_kmer(&tete_liste_kmer3, &tete_liste_sequence3, &tete_liste_pos3);
		printf("affichage_dictionnaire_kmer\n");

		affichage_motif_selectionne(&tete_liste_kmer_selectionne2, &tete_liste_motif_PSSM2);
		printf("affichage_motif_selectionné\n");

		calcul_PSSM(&tete_liste_kmer_selectionne_pour_calcul, &tete_liste_motif_PSSM_pour_calcul, &file_info, &matrice_PSSM);

		do //repeter l'amélioration de la PSSM jusqu'à convergence
		{
			while (p_generation_seq != NULL) //pour chaque sequence
			{
				pos_max= -1;
				score_max= -100;
				//pour chaque mot:
				while (position<=25)
				{
					n_sequence= p_generation_seq->numero_sequence;
					for (cpt_mot=0; cpt_mot<=longueur_masque; cpt_mot++)
					{
						p_mot->mot[cpt_mot]= p_generation_seq-> sequence[position];
						if (cpt_mot==longueur_masque)
						{
							p_mot->mot[longueur_masque]= '\0';
							calcul_score(&p_mot, &matrice_PSSM, n_sequence, &p_generation_seq, longueur_masque);
						}
						position++;
					}
					position= position -5;
					//si le score de ce mot est supérieur au score max
					if (p_mot->score_mot> score_max)
					{
						score_max= p_mot-> score_mot;
						pos_max= position;
						strcpy(p_mot_selected->motif, p_mot->mot);
					}
					//p_mot->seq=n_sequence;
					TPtr_Mot_Ameliorer_PSSM p_nouv_mot=malloc(sizeof(TMot_Ameliorer_PSSM));
					p_mot->next_mot= p_nouv_mot;
					p_mot=p_nouv_mot;
				}
				//Post condition: le mot de la sequence courante le plus proche de la PSSM est identifié !
				position=0;
				TPtr_Cell_Motif_PSSM p_mot_selected_suiv= malloc(sizeof(TCell_Motif_PSSM));
				p_mot_selected-> suiv_motif= p_mot_selected_suiv;
				p_mot_selected=p_mot_selected_suiv;
				score_max=-100;
				pos_max=0;
				p_generation_seq= p_generation_seq-> next_sequence;
			}
			calcul_nouvelle_PSSM(&tete_mot_selected_calcul_PSSM, &matrice_PSSM_nouv, nb_sequence, &Ct);
			printf("Ct: %s \n", Ct);
			distance_PSSM=dist_PSSM(&matrice_PSSM, &matrice_PSSM_nouv, &distance_PSSM);
			if (distance_PSSM>0.8)
			{
				for (cpt=0; cpt<4; cpt++) //ancienne_matrice= nouv_matrice.
				{
  				for(k=0; k<taille_motif; k++)
  				{
						(matrice_PSSM)[cpt][k]= (matrice_PSSM_nouv)[cpt][k];
  				}
				}
				for (cpt=0; cpt<4; cpt++) //reinitialisation de matrice_PSSM_nouv à 0.
				{
  				for(k=0; k<taille_motif; k++)
  				{
						(matrice_PSSM_nouv)[cpt][k]= 0;
  				}
				}
				p_generation_seq=tete_generation_sequence;
				p_mot_selected=tete_mot_selected;
				p_mot=tete_mot;
			}
		}while(distance_PSSM>0.8);
		//Post-condition: tous les mots de longueur l maximisant le score sont contenu dans la structure chainée p_mot_selected.
		//Le motif consensus de cette structure a été calculé
		//Il faut maintenant calculer la distance de Hamming entre le motif consensus et tous les mots de longueur l contenu dans la structure p_mot
		p_mot_selected= tete_mot_selected;
		//RAFFINER - Version 1:
		distanceHammingSt1(&Ct, &p_mot_selected, &p_st1);
		// RAFFINER - Version 2:
		p_mot_selected= tete_mot_selected;
		//do{
			distanceHammingSt2(&Ct, &p_mot_selected, &p_st2);
			//T' (p_mot_selected_prim) <- mot mi de longueur l de Si, mi minimisant Dh(mi, Ct)
			TPtr_Mot_Ameliorer_PSSM p_mot_st2_prim=tete_mot_pour_st2_prim;
			p_generation_seq=tete_generation_sequence;
			distanceHammingSt2_prim(&Ct, &p_generation_seq, &p_mot_st2_prim, &p_st2_prim);


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
