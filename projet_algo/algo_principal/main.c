#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <time.h>
#include "fonctions.h"

bool convergence;
//definition des variables principales
int taille_motif; //longueur motif/masque IE. l
int d; // nombre maximal de substitutions autorisées
int nb_fenetre; //nombre de fenêtres dans le masque IE. k
int nb_masques;
/*************************************************************************************************************
 * * *                                    MAIN                                                           * * *
 *************************************************************************************************************/
int main()
{
	srand(time(NULL));
	//définition des variables utilitaires
	int i, j;
	int position = 0; // position de l'occ dans la sequence
	int pos_max;
	int cpt_mot;  //?
  int n_sequence; // numéro de séquence (cf structure dico seq)
  char Ct[taille_motif];
	int masque[taille_motif];
	memset(masque, 0, sizeof masque);
	double score_max;
	double distance_PSSM; //distance entre deux matrice PSSM
  int convergence = 0; //convergence = 0 tant que T(p_mot_selected) est différent de T'(p_mot_selected_prim) et est égal à 1 quand la convergence est atteinte
	int st1, st2, st2_prim;

	/*******************************
	 *        dictionnaire         *
	 *******************************/
//initialisation des structures
	TInfo_ensemble_sequences structure_info_seq; //instanciation de la structure contenants les informations des sequences
	TPtr_info_ensemble_sequences ptr_info = malloc(sizeof(TInfo_ensemble_sequences));
	TPtr_ensemble_sequences ptr_ensemble = malloc(sizeof(TEnsemble_Sequences));
	ptr_info->tete_ensemble_seq = ptr_ensemble; //tete dictionnaire des séquences

	// TPtr_ensemble_sequences element_dictionnaire_sequence = malloc(sizeof(TEnsemble_Sequences)); // On cr�� le premier �l�ment de la structure
	TPtr_ensemble_sequences element_dictionnaire_sequence = ptr_ensemble;
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
  TPtr_Cell_Motif_PSSM tete_mot_selected = malloc(sizeof(TCell_Motif_PSSM)); // on créé une nouvelle tete pour recuperer les mots pour ameliorer la PSSM
  TPtr_Cell_Motif_PSSM p_mot_selected = tete_mot_selected;
  TPtr_Cell_Motif_PSSM tete_mot_selected_calcul_PSSM = tete_mot_selected;
  element_liste_motif_PSSM->suiv_motif = NULL;
  TPtr_Mot_Ameliorer_PSSM tete_mot = (TMot_Ameliorer_PSSM*)malloc(sizeof(TMot_Ameliorer_PSSM));
  TPtr_Mot_Ameliorer_PSSM p_mot = tete_mot;
  TPtr_Mot_Ameliorer_PSSM tete_mot_pour_st2_prim = (TMot_Ameliorer_PSSM*)malloc(sizeof(TMot_Ameliorer_PSSM));

  //VARIABLES SCORE ST1 ST2:
  Ptr_st tete_st1 = malloc(sizeof(st));
  Ptr_st p_st1 = tete_st1;

  Ptr_st tete_st2 = malloc(sizeof(st));
  Ptr_st p_st2 = tete_st2;

  Ptr_st tete_st2_prim = malloc(sizeof(st));
  Ptr_st p_st2_prim = tete_st2_prim;


	//MATRICE PSSM ET MOTIF CONSENSUS:
  double** matrice_PSSM;
  double** matrice_PSSM_nouv;
  int nb_ligne = 4;

  matrice_PSSM = malloc(nb_ligne * sizeof(*matrice_PSSM));
  for(i=0 ; i<nb_ligne ; i++)
  {
    matrice_PSSM[i] = malloc(taille_motif * sizeof(*(matrice_PSSM[i]))); //On alloue des tableaux de 'taille2' variables.
  }
  for (i=0; i<nb_ligne; i++) //initialisation de la matrice à 0:
  {
    for(j=0; j<taille_motif; j++)
    {
      (matrice_PSSM)[i][j] = 0;
    }
  }
  matrice_PSSM_nouv = malloc(4 * sizeof(*matrice_PSSM_nouv));
  for(i=0 ; i<nb_ligne ; i++)
  {
    matrice_PSSM_nouv[i] = malloc(taille_motif * sizeof(*(matrice_PSSM_nouv[i]))); //On alloue des tableaux de 'taille2' variables.
  }
  for (i=0; i<nb_ligne; i++) //initialisation de la matrice à 0:
  {
    for(j=0; j<taille_motif; j++)
    {
      (matrice_PSSM_nouv)[i][j]= 0;
    }
  }
	/*************************************************************************************************************
	 * * *                                    DEBUT DE L'ALGORITHME                                          * * *
	 *************************************************************************************************************/
	importer_parametres(&taille_motif, &d, &nb_fenetre, &nb_masques);
	printf("taille_motif= %d, ",taille_motif);printf("d= %d, ",d);
	printf("nb_fenetre= %d, ",nb_fenetre);printf("nb_masques= %d\n",nb_masques);

	importer_sequences_fasta(&ptr_info, &ptr_ensemble );
  //afficher_sequences(&ptr_info, &ptr_ensemble );

	for ( i = 0; i <= nb_masques; i++ )
	{
		generation_masque(taille_motif, &masque, nb_fenetre);
		// creation_dictionnaire();
		parcours_masque(taille_motif, &masque, nb_fenetre, ptr_info->nb_sequences, &tete_liste_pour_parcours_masque, &tete_liste_kmer, &tete_liste_sequence, &tete_liste_pos);

		kmer_present_dans_chaque_sequence(ptr_info->nb_sequences, &tete_liste_kmer2, &tete_liste_sequence2, &tete_liste_pos2, &tete_liste_pour_recup_motif, &tete_liste_kmer_selectionne, &tete_liste_motif_PSSM);

		affichage_dictionnaire_kmer(&tete_liste_kmer3, &tete_liste_sequence3, &tete_liste_pos3);

		affichage_motif_selectionne(&tete_liste_kmer_selectionne2, &tete_liste_motif_PSSM2);

	// // On calculera la PSSM seulement pour les kmers qui sont présent dans plus de 7 sequences:
	// 	while (tete_liste_kmer4 != NULL)
	// 	{
	// 		if (tete_liste_kmer4->nb_sequence >= 7)
	// 		{
	// 			calcul_PSSM(&tete_liste_kmer_selectionne_pour_calcul, &tete_liste_motif_PSSM_pour_calcul, &matrice_PSSM, taille_motif);
	// 			do //repeter l'amélioration de la PSSM jusqu'à convergence
	// 			{
	// 				while (p_dictionnaire_seq != NULL) //pour chaque sequence
	// 				{
	// 					pos_max = -1;
	// 					score_max = -100;
	// 					//pour chaque mot:
	// 					while (position <= 25)
	// 					{
	// 						n_sequence = p_dictionnaire_seq->numero_sequence;
	// 						for (cpt_mot=0; cpt_mot <= taille_motif; cpt_mot++)
	// 						{
	// 							p_mot->mot[cpt_mot] = p_dictionnaire_seq->sequence[position];
	// 							if (cpt_mot == taille_motif)
	// 							{
	// 								p_mot->mot[taille_motif] = '\0';
	// 								calcul_score(&p_mot, &matrice_PSSM, n_sequence, &p_dictionnaire_seq, taille_motif);
	// 								afficher_PSSM( &matrice_PSSM, taille_motif);
	// 							}
	// 							position++;
	// 						}
	// 						position = position -5;
	// 						//si le score de ce mot est supérieur au score max
	// 						if (p_mot->score_mot > score_max)
	// 						{
	// 							score_max = p_mot->score_mot;
	// 							pos_max = position;
	// 							strcpy(p_mot_selected->motif, p_mot->mot);
	// 						}
	// 						//p_mot->seq=n_sequence;
	// 						TPtr_Mot_Ameliorer_PSSM p_nouv_mot = malloc(sizeof(TMot_Ameliorer_PSSM));
	// 						p_mot->next_mot = p_nouv_mot;
	// 						p_mot = p_nouv_mot;
	// 					}
	// 					//Post condition: le mot de la sequence courante le plus proche de la PSSM est identifié !
	// 					position = 0;
	// 					TPtr_Cell_Motif_PSSM p_mot_selected_suiv = malloc(sizeof(TCell_Motif_PSSM));
	// 					p_mot_selected->suiv_motif = p_mot_selected_suiv;
	// 					p_mot_selected = p_mot_selected_suiv;
	// 					score_max = -100;
	// 					pos_max = 0;
	// 					p_dictionnaire_seq = p_dictionnaire_seq->suiv_seq;
	// 				}
	// 				calcul_nouvelle_PSSM(&tete_mot_selected_calcul_PSSM, &matrice_PSSM_nouv, ptr_info->nb_sequences, &Ct, taille_motif);
	// 				printf("Ct: %s \n", Ct);
	// 				distance_PSSM = dist_PSSM(&matrice_PSSM, &matrice_PSSM_nouv, &distance_PSSM, taille_motif);
	// 				if (distance_PSSM > 0.8)
	// 				{
	// 					for (i=0; i<nb_ligne; i++) //ancienne_matrice= nouv_matrice.
	// 					{
	// 						for(j=0; j<taille_motif; j++)
	// 						{
	// 							(matrice_PSSM)[i][j] = (matrice_PSSM_nouv)[i][j];
	// 						}
	// 					}
	// 					for (i=0; i<nb_ligne; i++) //reinitialisation de matrice_PSSM_nouv à 0.
	// 					{
	// 						for(j=0; j<taille_motif; j++)
	// 						{
	// 							(matrice_PSSM_nouv)[i][j] = 0;
	// 						}
	// 					}
	// 					p_dictionnaire_seq = ptr_ensemble;
	// 					p_mot_selected = tete_mot_selected;
	// 					p_mot = tete_mot;
	// 				}
	// 			}while(distance_PSSM > 0.8);
	// 			//Post-condition: tous les mots de longueur l maximisant le score sont contenu dans la structure chainée p_mot_selected.
	// 			//Le motif consensus de cette structure a été calculé
	// 			//Il faut maintenant calculer la distance de Hamming entre le motif consensus et tous les mots de longueur l contenu dans la structure p_mot
	// 			p_mot_selected = tete_mot_selected;
  //
	// 			//RAFFINER - Version 1:
	// 			st1 = distanceHammingSt1(&Ct, &p_mot_selected, &p_st1, taille_motif);
	// 			int v_St1_Pos[9]; //Position dans la liste chainée
	// 			quick_sort_ST(&p_st1, v_St1_Pos, ptr_info->nb_sequences);
	// 			fichier_sortie_st(&p_st1, v_St1_Pos, &Ct, ptr_info->nb_sequences);
	// 			printf("ST1= %d \n", st1);
  //
	// 			// RAFFINER - Version 2:
	// 			p_mot_selected = tete_mot_selected;
	// 			//T' (p_st2_prim) <- mot mi de longueur l de Si, mi minimisant Dh(mi, Ct)
	// 			TPtr_Mot_Ameliorer_PSSM p_mot_st2_prim = tete_mot_pour_st2_prim;
	// 			do{
	// 				st2 = distanceHammingSt2(&Ct, &p_mot_selected, &p_st2, taille_motif);
	// 				p_dictionnaire_seq = ptr_ensemble;
	// 				p_mot_st2_prim = tete_mot_pour_st2_prim;
	// 				st2_prim = distanceHammingSt2_prim(&Ct, &p_dictionnaire_seq, &p_mot_st2_prim, &p_st2_prim, taille_motif);
	// 				printf("St2: %d, St2 prim%d ", st2, st2_prim);
	// 				p_st2_prim = tete_st2_prim;
	// 				if (st2_prim > st2)
	// 				{
	// 					//Avant de faire ça il faudrait free l'ancienne liste st2;
	// 					TPtr_Cell_Motif_PSSM new_tete_mot_selected = malloc(sizeof(TCell_Motif_PSSM));
	// 					p_mot_selected = new_tete_mot_selected;
	// 					// T<-T';
	// 					while (p_st2_prim != NULL)
	// 					{
	// 						if (p_st2_prim->mot != NULL)
	// 						{
	// 							if (p_st2_prim->distance_hamming < 2)
	// 							{
	// 								strcpy(p_mot_selected->motif, p_st2_prim->mot);
	// 								printf("Nouveau ST2: %s \n", p_mot_selected->motif);
	// 								TPtr_Cell_Motif_PSSM p_next_mot_selected = malloc(sizeof(TCell_Motif_PSSM));
	// 								p_mot_selected->suiv_motif = p_next_mot_selected;
	// 								p_mot_selected = p_next_mot_selected;
	// 							}
	// 						}
	// 						p_st2_prim = p_st2_prim->next_mot;
	// 					}
	// 					//Calcul du nouveau motif consensus à partir de la nouvelle liste T
	// 					p_mot_selected = new_tete_mot_selected;
	// 					calcul_nouvelle_PSSM(&p_mot_selected, &matrice_PSSM_nouv, ptr_info->nb_sequences, &Ct, taille_motif);
	// 					p_mot_selected = new_tete_mot_selected;
	// 					Ptr_st new_tete_st2_prim = malloc(sizeof(st2));
	// 					p_st2_prim = new_tete_st2_prim;
	// 				}else{
	// 					convergence = 1;
	// 				}
	// 			}while(convergence == 0);
	// 			int v_St2_Pos[9]; //Position dans la liste chainée
	// 			quick_sort_ST(&p_st2, v_St2_Pos, ptr_info->nb_sequences);
	// 			fichier_sortie_st(&p_st2, v_St2_Pos, &Ct, ptr_info->nb_sequences);
	// 		}
	// 		tete_liste_kmer4 = tete_liste_kmer4->suiv_kmer;
	// 	}
	// }
	return(0);
}
