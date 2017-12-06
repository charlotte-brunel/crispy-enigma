#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "fonctions.h"
//------------------------------------------------------------------------------------------------------------
// lecture du fichier contenant les variables souhaitées, stockage de ces variables pour utilisation
void importer_parametres(int* taille_motif, int* d, int* nb_fenetres, int* nb_masques)
{
	FILE* ptr_fichier; //creation d'un pointeur sur le fichier
  char nom_fichier[30];

	printf("Quel est le nom du fichier que vous voulez utiliser pour importer les paramètres du programme ?\n");
	scanf("%s", nom_fichier);
  ptr_fichier = fopen(nom_fichier, "r"); //ouverture du fichier

	if( ptr_fichier == NULL) { free(ptr_fichier); return; }// si le fichier est vide on sort de la fonction

  fscanf(ptr_fichier, "#longueur du motif à identifier (IE. taille du masque): %d\n", taille_motif);
  fscanf(ptr_fichier, "#nombre maximal de substitutions autorisées: %d\n", d);
  fscanf(ptr_fichier, "#nombre de fenêtres dans les masques utilisés: %d\n", nb_fenetres);
  fscanf(ptr_fichier, "#nombre de masques à générer: %d\n", nb_masques);
  fclose(ptr_fichier); //fermeture du fichier
}
//------------------------------------------------------------------------------------------------------------
// Récupère les séquences et les stocke dans une liste chainée
void importer_sequences_fasta( TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq )
{
  FILE* ptr_fichier_fasta;
  char nom_fichier_fasta[30];
  char contenu_ligne[128];
  char c; // utiliser pour lire le fichier caractère par caractère
  int cptr = 0; //compteur du nombre de seq dans le fichier

	// déclaration des pointeurs permettant de créer la liste chainée de séquences
  TPtr_dictionnaire_sequences p_new = NULL;
  TPtr_dictionnaire_sequences p = *adr_tete_dict_seq;
	TPtr_info_dictionnaire_sequences ptr_info = *adr_tete_info_dict_seq;
  ptr_info->nb_sequences = 0; //initialisation du nombre de seq contenues dans le fichier à 0

	printf("Quel est le nom du fichier que vous voulez utiliser pour importer les séquences ?\n");
  scanf("%s", nom_fichier_fasta);

  ptr_fichier_fasta = fopen(nom_fichier_fasta, "r");
  if( ptr_fichier_fasta == NULL) { free(ptr_fichier_fasta); return; }// si le fichier est vide on sort de la fonction

	do //cette boucle permet de compter les séquences
	{
		c = fgetc (ptr_fichier_fasta); //lecture du fichier caractère par caractère
		if (c == '>') cptr++;  //incrémentation pour chaque nom de séquences
	} while (c != EOF); // l'action est répétée jusqu'à ce que la fin du fichier soit atteinte

	ptr_fichier_fasta = fopen(nom_fichier_fasta, "r");
	do
	{
		fscanf(ptr_fichier_fasta,"%s\n", contenu_ligne); // lecture du fichier ligne par ligne
	  if ( contenu_ligne[0] == '>')
	  {
			if ((ptr_info->nb_sequences) > 0)
			{
				p_new = malloc ( sizeof(TDictionnaire_Sequences));
				p->suiv_seq = p_new;
				p = p_new;
				if ( ptr_info->nb_sequences  == cptr) { p->suiv_seq = NULL;}
			}
			ptr_info->nb_sequences += 1;
			strcpy(p->nom_seq, contenu_ligne);
			p->numero_sequence = ptr_info->nb_sequences;
		}else{
			if ( p->sequence == NULL){
				strcpy(p->sequence, contenu_ligne);
	    }else{
	      strcat(p->sequence, contenu_ligne);
	    }
	  }
	}	while(!feof(ptr_fichier_fasta));
  fclose(ptr_fichier_fasta);
}
//------------------------------------------------------------------------------------------------------------
void afficher_sequences(TPtr_info_dictionnaire_sequences* ptr_info, TPtr_dictionnaire_sequences* ptr_ensemble )
{
	FILE* ptr_fichier;
  TPtr_dictionnaire_sequences p = *ptr_ensemble;
	int cptr = 0;
	printf("fonction afficher_sequences\n");
	printf("%d \n", (*ptr_info)->nb_sequences);

	ptr_fichier = fopen ( "verif_dico_fasta.txt" ,"w");
	while ( cptr < ((*ptr_info)->nb_sequences))
	{
		fprintf(ptr_fichier, "%d\n", p->numero_sequence);
		fputs( p->nom_seq, ptr_fichier);
		fputs( "\n" , ptr_fichier);
		fputs( p->sequence, ptr_fichier);
		fputs( "\n", ptr_fichier);
		printf(" %s\n", p->nom_seq );
		printf("%s\n", p->sequence);
		p = p->suiv_seq;
		cptr ++;
	}
}
//------------------------------------------------------------------------------------------------------------
int random_number(int max_number, int zero_excluded)
{
	int randomNumber;
	if(zero_excluded == 0) //on peut tomber sur 0 aléatoirement
  {
		randomNumber = rand() % max_number;
	}else{
    randomNumber = rand() % max_number +1;
	}
	return(randomNumber);
}
//------------------------------------------------------------------------------------------------------------
void generation_masque(void* adr_masque)
{
  int i;
  int *p_masque = adr_masque;
  int nb_fenetres_ouverte = 0;

  while (nb_fenetres_ouverte != nb_fenetres)
  {
	nb_fenetres_ouverte=0;
    for (i=0; i <= (taille_motif-1); i++)
    {
      p_masque[i] = random_number(2,0);
      if (p_masque[i] == 1) {  nb_fenetres_ouverte ++;  }
    }
  }
}
// //------------------------------------------------------------------------------------------------------------
void parcours_masque( void* adr_masque, TPtr_dictionnaire_sequences* adr_tete_dict_seq, TPtr_Cellkmer* adr_tete_liste_kmer, TPtr_CellSequence* adr_tete_liste_sequence, TPtr_CellPos* adr_tete_liste_pos)
{
	TPtr_dictionnaire_sequences p_dictionnaire_sequence = *adr_tete_dict_seq;
  TPtr_Cellkmer p_kmer = *adr_tete_liste_kmer;
  TPtr_CellSequence p_sequence = *adr_tete_liste_sequence;
  TPtr_CellPos p_pos = *adr_tete_liste_pos;
	
  int *p_masque = adr_masque;
  int position = 0;
  int cpt_pos_masque, position_kmer;
  int pos_kmer = 0;
  char k_mer[nb_fenetres+1]; // le k_mer mesure la taille du nombre de fen�tre ouverte dans le masque

  while (p_dictionnaire_sequence != NULL)
  {
    while (position < 30)
    {
      for(cpt_pos_masque = 0; cpt_pos_masque< taille_motif; cpt_pos_masque++)
      { // on parcourt les nucleotides sous le masque
        if (p_masque[cpt_pos_masque] == 1)
        { //si la fen�tre du masque est ouverte
          k_mer[pos_kmer] = p_dictionnaire_sequence->sequence[position];

          if (pos_kmer == 0) { position_kmer = position; }
          pos_kmer ++;
          if (pos_kmer == nb_fenetres)
          {
            k_mer[pos_kmer] = '\0' ;
            pos_kmer = 0;
          }
        }
        position = position+1;
      }
      generation_kmer(position_kmer, k_mer, &p_kmer, &p_sequence, &p_dictionnaire_sequence, &p_pos);
    }
    position = 0;
    p_dictionnaire_sequence = p_dictionnaire_sequence->suiv_seq;
  }
  return;
}
// //------------------------------------------------------------------------------------------------------------
// void generation_kmer(int position_kmer, char* k_mer, TPtr_Cellkmer* adr_liste_kmer, TPtr_CellSequence* adr_liste_sequence, TPtr_dictionnaire_sequences* ptr_ensemble, TPtr_CellPos* adr_liste_pos)
// {
//   TPtr_Cellkmer temp_p_kmer = *adr_liste_kmer; //on cr�� des pointeurs temporaire pour parcourir les listes
//   TPtr_dictionnaire_sequences p_dictionnaire_sequence = *ptr_ensemble;
//
//   if (temp_p_kmer->suiv_kmer == NULL) //cas du premier element de la liste
//   {
//     TPtr_Cellkmer nouveau_kmer = (TCellkmer*)malloc(sizeof(TCellkmer));
//     strcpy(temp_p_kmer->kmer, k_mer);
//     nouveau_kmer->suiv_kmer = NULL;
//     TPtr_CellSequence nouvelle_sequence = (TCellSequence*)malloc(sizeof(TCellSequence));
//     temp_p_kmer->suiv_kmer = nouveau_kmer;
//     temp_p_kmer->tete_sequence = nouvelle_sequence;
//     nouvelle_sequence->sequence = p_dictionnaire_sequence->numero_sequence;
//     nouvelle_sequence->suiv_sequence = NULL;
//     TPtr_CellPos nouvelle_position = (TCellPos*)malloc(sizeof(TCellPos));
//     nouvelle_sequence->tete_pos = nouvelle_position;
//     nouvelle_position->position = position_kmer;
//     nouvelle_position->suiv_pos = NULL;
//     return;
//   }else{
//     while (temp_p_kmer->suiv_kmer != NULL)
//     {
//       if (strcmp(temp_p_kmer->kmer, k_mer) == 0)   //si le kmer a d�j� �t� trouv�
//       {
//         TPtr_CellSequence p_liste_sequence = temp_p_kmer->tete_sequence;
//         //premier element de la liste chainee de sequence:
//         while (p_liste_sequence->suiv_sequence != NULL)
//         {
//           if (p_dictionnaire_sequence->numero_sequence == p_liste_sequence->sequence) //Si le kmer a d�j� �t� trouv� dans cette s�quence
//           {
//             TPtr_CellPos p_liste_pos = p_liste_sequence->tete_pos;
//             while (p_liste_pos->suiv_pos != NULL){
//               p_liste_pos = p_liste_pos->suiv_pos;
//             }
//             TPtr_CellPos nouvelle_position = (TCellPos*)malloc(sizeof(TCellPos)); // On cr�� une nouvelle brique de p�sition
//             nouvelle_position->position = position_kmer;
//             nouvelle_position->suiv_pos = NULL;
//             p_liste_pos->suiv_pos = nouvelle_position;
//             return;
//           }
//           p_liste_sequence = p_liste_sequence->suiv_sequence;
//         }
//         if (p_dictionnaire_sequence->numero_sequence == p_liste_sequence->sequence)  //Si le kmer a d�j� �t� trouv� dans cette s�quence
//         {
//           TPtr_CellPos p_liste_pos = p_liste_sequence->tete_pos;
//           while (p_liste_pos->suiv_pos != NULL)
//           {
//             p_liste_pos = p_liste_pos->suiv_pos;
//           }
//           TPtr_CellPos nouvelle_position = (TCellPos*)malloc(sizeof(TCellPos)); // On cr�� une nouvelle brique de p�sition
//           nouvelle_position->position = position_kmer;
//           nouvelle_position->suiv_pos = NULL;
//           p_liste_pos->suiv_pos = nouvelle_position;
//           return;
//         }
//         TPtr_CellSequence nouvelle_sequence = (TCellSequence*)malloc(sizeof(TCellSequence));
//         nouvelle_sequence->sequence = p_dictionnaire_sequence->numero_sequence;
//         nouvelle_sequence->suiv_sequence = NULL;
//         p_liste_sequence->suiv_sequence = nouvelle_sequence;
//         TPtr_CellPos nouvelle_tete_pos = (TCellPos*)malloc(sizeof(TCellPos)); // on cr�� une nouvelle tete de liste de position
//         nouvelle_sequence->tete_pos = nouvelle_tete_pos;
//         nouvelle_tete_pos->position = position_kmer;
//         nouvelle_tete_pos->suiv_pos = NULL;
//         return;
//       }
//       temp_p_kmer = temp_p_kmer->suiv_kmer;
//     }
//    // Si ce kmer n'a pas encore �t� trouv�: ajout en fin de boucle
//     TPtr_Cellkmer nouveau_kmer = (TCellkmer*)malloc(sizeof(TCellkmer));
//     strcpy(temp_p_kmer->kmer, k_mer);
//     nouveau_kmer->suiv_kmer = NULL;
//     TPtr_CellSequence nouvelle_sequence = (TCellSequence*)malloc(sizeof(TCellSequence));
//     temp_p_kmer->suiv_kmer = nouveau_kmer;
//     temp_p_kmer->tete_sequence = nouvelle_sequence;
//     nouvelle_sequence->sequence = p_dictionnaire_sequence->numero_sequence;
//     nouvelle_sequence->suiv_sequence = NULL;
//     TPtr_CellPos nouvelle_position = (TCellPos*)malloc(sizeof(TCellPos));
//     nouvelle_sequence->tete_pos = nouvelle_position;
//     nouvelle_position->position = position_kmer;
//     nouvelle_position->suiv_pos = NULL;
//     return;
//   }
// }
//
// //------------------------------------------------------------------------------------------------------------
// //Procedure pour remplir le dictionnaire de kmer s�lectionn� pour le calcul de la PSSM:
// void recuperer_motif_kmer(TPtr_Cellkmer* adr_parcours_kmer, TPtr_Cellkmer_selectionne *adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif_PSSM, TPtr_CellSequence* adr_cell_sequence, TPtr_CellPos* adr_cell_pos, TPtr_dictionnaire_sequences* ptr_ensemble, int nb_sequence_kmer)
// {
//   TPtr_Cellkmer p_kmer = *adr_parcours_kmer;
//   TPtr_Cellkmer_selectionne p_kmer_selectionne = *adr_tete_kmer_selectionne;
//   TPtr_dictionnaire_sequences tete_dictionnaire_sequence = *ptr_ensemble;
//   TPtr_dictionnaire_sequences p_dictionnaire_sequence = tete_dictionnaire_sequence;
//   int longueur_sequence = 30;
//   int longueur_motif = taille_motif;
//   int i, sequence_actuelle;
//
//   if (p_kmer_selectionne->suiv_kmer_selectionne == NULL) //la liste est vide:
//   {
//     TPtr_Cellkmer_selectionne prochain_kmer = malloc(sizeof(TCellkmer_selectionne)); //Je créer la prochaine brique
//     TPtr_Cell_Motif_PSSM nouveau_motif = malloc(sizeof(TCell_Motif_PSSM));
//     strcpy(p_kmer_selectionne->kmer, p_kmer->kmer);//Je remplis la brique actuelle
//     p_kmer_selectionne->nb_sequence = nb_sequence_kmer;
//     p_kmer_selectionne->suiv_kmer_selectionne = prochain_kmer; //je fais le chainage entre la brique actuelle et la prochaine
//     p_kmer_selectionne->tete_motif_PSSM = nouveau_motif; //Je recupere la bonne tete
//     TPtr_CellSequence p_sequence = p_kmer->tete_sequence;
//     do
//     {
//       TPtr_CellPos p_pos = p_sequence->tete_pos; // je recupere la liste de position pour chaque sequence
//       sequence_actuelle = p_sequence->sequence;
//       while (p_dictionnaire_sequence->numero_sequence != sequence_actuelle)  //je me place dans la bonne sequence
//       {
//         p_dictionnaire_sequence = p_dictionnaire_sequence->suiv_seq;
//       }
//       if((p_pos->position+longueur_motif) < longueur_sequence)
//       {
//         for (i = (p_pos->position); i < (p_pos->position+longueur_motif); i++)
//         {
//           nouveau_motif->motif[i-(p_pos->position)] = p_dictionnaire_sequence->sequence[i];
//         }
//         nouveau_motif->motif[longueur_motif] = '\0';
//         TPtr_Cell_Motif_PSSM prochain_motif = (TCell_Motif_PSSM*)malloc(sizeof(TCell_Motif_PSSM));
//         nouveau_motif->suiv_motif = prochain_motif;
//         nouveau_motif = prochain_motif;
//       }
//       p_sequence = p_sequence->suiv_sequence;
//     }while(p_sequence != NULL);
//     p_dictionnaire_sequence = tete_dictionnaire_sequence;
//     return;
//   }
//   do //Cas o� on a plusieurs kmer selectionn�
//   {
//     p_kmer_selectionne = p_kmer_selectionne->suiv_kmer_selectionne;
//   }while (p_kmer_selectionne->suiv_kmer_selectionne != NULL);
//
//   TPtr_Cellkmer_selectionne prochain_kmer = (TCellkmer_selectionne*)malloc(sizeof(TCellkmer_selectionne));
//   TPtr_Cell_Motif_PSSM nouveau_motif = (TCell_Motif_PSSM*)malloc(sizeof(TCell_Motif_PSSM));
//   strcpy(p_kmer_selectionne->kmer, p_kmer->kmer);
//   p_kmer_selectionne->nb_sequence = nb_sequence_kmer;
//   p_kmer_selectionne->suiv_kmer_selectionne = prochain_kmer;
//   p_kmer_selectionne->tete_motif_PSSM = nouveau_motif;
//   TPtr_CellSequence p_sequence = p_kmer->tete_sequence;
//   do
//   {
//     TPtr_CellPos p_pos = p_sequence->tete_pos; // je recupere la liste de position pour chaque sequence
//     sequence_actuelle = p_sequence->sequence;
//     while (p_dictionnaire_sequence->numero_sequence != sequence_actuelle)
//     {
//       p_dictionnaire_sequence = p_dictionnaire_sequence->suiv_seq;
//     } //je me place dans la bonne sequence
//
//     if((p_pos->position+longueur_motif) < longueur_sequence)
//     {
//       for (i = (p_pos->position); i < (p_pos->position+longueur_motif); i++)
//       {
//         nouveau_motif->motif[i-(p_pos->position)] = p_dictionnaire_sequence->sequence[i];
//       }
//       nouveau_motif->motif[longueur_motif] = '\0';
//
//       TPtr_Cell_Motif_PSSM prochain_motif = malloc(sizeof(TCell_Motif_PSSM));
//       nouveau_motif->suiv_motif = prochain_motif;
//       nouveau_motif = prochain_motif;
//     }
//     p_sequence = p_sequence->suiv_sequence;
//   }while(p_sequence != NULL);
//
//   p_dictionnaire_sequence = tete_dictionnaire_sequence;
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// //Fonction pour trouver si le kmer est trouv� dans chaque sequence:
// void kmer_present_dans_chaque_sequence(int nb_sequence, TPtr_Cellkmer* adr_cell_kmer, TPtr_CellSequence *adr_cell_sequence, TPtr_CellPos *adr_cell_pos, TPtr_dictionnaire_sequences* ptr_ensemble, TPtr_Cellkmer_selectionne* adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_cell_motif_PSSM)
// {
//   TPtr_Cellkmer p_parcours_kmer = *adr_cell_kmer;
//   TPtr_Cellkmer_selectionne tete_kmer_selectionne = *adr_cell_kmer_selectionne;
//   TPtr_Cellkmer_selectionne p_kmer_selectionne = *adr_cell_kmer_selectionne;
//   p_kmer_selectionne->suiv_kmer_selectionne = NULL;
//   TPtr_Cell_Motif_PSSM p_motif_PSSM = *adr_cell_motif_PSSM;
//   TPtr_CellSequence tete_sequence = *adr_cell_sequence;
//   TPtr_CellPos tete_pos = *adr_cell_pos;
//   TPtr_dictionnaire_sequences tete_dictionnaire_sequence = *ptr_ensemble;
//   int nb_sequence_par_kmer = 0;
//
//   while (p_parcours_kmer != NULL)
//   {
//     TPtr_CellSequence p_parcours_sequence = p_parcours_kmer->tete_sequence;
//     while (p_parcours_sequence != NULL )
//     {
//       nb_sequence_par_kmer ++;
//       p_parcours_sequence = p_parcours_sequence->suiv_sequence;
//     }
//     p_parcours_kmer->nb_sequence = nb_sequence_par_kmer;
//     if (nb_sequence_par_kmer >= 7) //Pour l"instant pour pouvoir continuer
//     {
//       recuperer_motif_kmer(&p_parcours_kmer, &p_kmer_selectionne, &p_motif_PSSM, &tete_sequence, &tete_pos, &tete_dictionnaire_sequence, nb_sequence_par_kmer);
//       p_kmer_selectionne = tete_kmer_selectionne;
//     }
//     nb_sequence_par_kmer = 0;
//     p_parcours_kmer = p_parcours_kmer->suiv_kmer;
//   }
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// void affichage_dictionnaire_kmer(TPtr_Cellkmer* adr_tete_kmer, TPtr_CellSequence* adr_tete_sequence, TPtr_CellPos* adr_tete_pos)
// {
//   TPtr_Cellkmer p_kmer = *adr_tete_kmer;
//   FILE* fichier_dictionnaire = NULL;
//   fichier_dictionnaire = fopen("dictionnaire_kmer.txt", "w");
//   while (p_kmer != NULL)
//   {
//     TPtr_CellSequence p_sequence = p_kmer->tete_sequence;
//     fprintf(fichier_dictionnaire, "KMER: %s \n", p_kmer->kmer);
//     while (p_sequence != NULL)
//     {
//       TPtr_CellPos p_pos = p_sequence->tete_pos;
//       fprintf(fichier_dictionnaire, "Sequence: %d \n", p_sequence->sequence);
//       while (p_pos != NULL)
//       {
//         fprintf(fichier_dictionnaire, "Position: %d \n ", p_pos->position);
//         p_pos = p_pos->suiv_pos;
//       }
//       p_sequence = p_sequence->suiv_sequence;
//     }
//     fprintf(fichier_dictionnaire, "Nb de séquence dans lesquelles le kmer est présent: %d \n", p_kmer->nb_sequence);
//     fprintf(fichier_dictionnaire, "\n\n\n");
//     p_kmer = p_kmer->suiv_kmer;
//   }
//   fclose(fichier_dictionnaire);
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// void affichage_motif_selectionne(TPtr_Cellkmer_selectionne* adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif)
// {
//   TPtr_Cellkmer_selectionne p_kmer_selectionne = *adr_tete_kmer_selectionne;
//   FILE* fichier_kmer_selectionne = NULL;
//   fichier_kmer_selectionne = fopen("kmer_selectionne.txt", "w");
//   while (p_kmer_selectionne != NULL)
//   {
//     TPtr_Cell_Motif_PSSM p_motif = p_kmer_selectionne->tete_motif_PSSM;
//     fprintf(fichier_kmer_selectionne, "KMER: %s \n", p_kmer_selectionne->kmer);
//     fprintf(fichier_kmer_selectionne, "Present dans %d sequences \n", p_kmer_selectionne->nb_sequence);
//     while(p_motif != NULL)
//     {
//       fprintf(fichier_kmer_selectionne, "MOTIF: %s \n", p_motif->motif);
//       p_motif = p_motif->suiv_motif;
//     }
//     p_kmer_selectionne = p_kmer_selectionne->suiv_kmer_selectionne;
//   }
//   fclose(fichier_kmer_selectionne);
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// void calcul_PSSM(TPtr_Cellkmer_selectionne *adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM *adr_cell_motif_PSSM, double*** adr_matrice_PSSM, int taille_motif)
// {
//   TPtr_Cellkmer_selectionne p_kmer_selectionne = *adr_cell_kmer_selectionne;
//   double nb_sequence = p_kmer_selectionne->nb_sequence;
//   double add = 1/nb_sequence;
//   double** p_matrice_PSSM = *adr_matrice_PSSM;
//   int i;
//   //calcul de la PSSM à partir de la liste chainée de motifs
//   TPtr_Cell_Motif_PSSM p_motif = p_kmer_selectionne->tete_motif_PSSM;
//   do
//   { // on remplit la PSSM
//     for (i=0; i < taille_motif; i++)
//     {
//       switch( p_motif->motif[i] )
//       {
//         case 'a':
//             p_matrice_PSSM[0][i] = p_matrice_PSSM[0][i] + add;
//             break;
//         case 't':
//             p_matrice_PSSM[1][i] = p_matrice_PSSM[1][i] + add;
//             break;
//         case 'c':
//             p_matrice_PSSM[2][i] = p_matrice_PSSM[2][i] + add;
//             break;
//         case 'g':
//             p_matrice_PSSM[3][i] = p_matrice_PSSM[3][i] + add;
//             break;
//       }
//     }
//     p_motif = p_motif->suiv_motif;
//   }while (p_motif != NULL);
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// void calcul_nouvelle_PSSM(TPtr_Cell_Motif_PSSM *adr_cell_mot_selected, double*** adr_matrice_PSSM_nouv, double nb_sequence, char (*adr_Ct)[6], int taille_motif)
// {
// 	FILE* file_nouv_PSSM;
// 	file_nouv_PSSM= fopen("Nouvelle_PSSM_Motif_Trouve.txt", "a"); // Dans ce fichier on va �crire la PSSM pr�visionnelle
//
//   double add = 1/nb_sequence;
//   TPtr_Cell_Motif_PSSM p_mot = *adr_cell_mot_selected;
//   double** p_matrice_PSSM_nouv = *adr_matrice_PSSM_nouv;
//
// 	char to_print;
// 	double maximum= -1;
// 	int i, j, k;
//   //calcul de la PSSM � partir de la liste chain�e de motif
//   do
//   { // on remplit la PSSM
//       for (i=0; i<taille_motif; i++)
//       {
//           switch(p_mot->motif[i])
//           {
//           case 'a':
//               p_matrice_PSSM_nouv[0][i]= p_matrice_PSSM_nouv[0][i] + add;
//               break;
//           case 't':
//               p_matrice_PSSM_nouv[1][i]= p_matrice_PSSM_nouv[1][i] + add;
//               break;
//           case 'c':
//               p_matrice_PSSM_nouv[2][i]= p_matrice_PSSM_nouv[2][i] + add;
//               break;
//           case 'g':
//               p_matrice_PSSM_nouv[3][i]= p_matrice_PSSM_nouv[3][i] + add;
//               break;
//           }
//
//       }
//       p_mot = p_mot->suiv_motif;
//   }while (p_mot != NULL);
//
// 	//ecriture du motif consensus dans le fichier a partir de la matrice PSSM:
// 	fprintf(file_nouv_PSSM, "\n\nMotif Consensus: \n" );
// 	for (j=0; j<taille_motif; j++)
// 	{
// 		for(k=0; k<4; k++)
// 		{
// 			if (p_matrice_PSSM_nouv[k][j]> maximum){
// 				maximum=p_matrice_PSSM_nouv[k][j];
// 				switch (k)
// 				{
// 					case 0:
// 							to_print= 'a';
// 							break;
// 					case 1:
// 							to_print= 't';
// 							break;
// 					case 2:
// 							to_print= 'c';
// 							break;
// 					case 3:
// 							to_print= 'g';
// 							break;
// 				}
// 			}
// 		}
// 		*adr_Ct[j]=to_print;
// 		printf("MOTIF CONSENSUS: %c", *adr_Ct[j]);
// 		maximum=-1; //on remet maximum a -1 avant d'�valuer la nucl�otide majoritaire de la s�quence suivante !
// 		fprintf(file_nouv_PSSM, "%c", to_print);
// 	}
// 	*adr_Ct[taille_motif+1]='\0';
// 	fprintf(file_nouv_PSSM, "\n\n\n");
// 	fclose(file_nouv_PSSM);
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// void afficher_PSSM( double*** adr_matrice_PSSM, int taille_motif)
// {
//   FILE* ptr_fichier_PSSM;
//   double** p_matrice_PSSM = *adr_matrice_PSSM;
//   int a, t, c, g;
//   ptr_fichier_PSSM = fopen("PSSM_Motif_Trouve.txt", "w"); // Dans ce fichier on va �crire la PSSM pr�visionnelle
//
//   //ecriture de la matrice PSSM dans le fichier d'info.
//   fprintf(ptr_fichier_PSSM, "\n\nPSSM: \n");
//   fprintf(ptr_fichier_PSSM, "a  ");
//   for (a=0; a < taille_motif; a++)
//   {
//     fprintf(ptr_fichier_PSSM, "%.2f ", p_matrice_PSSM[0][a]);
//   }
//   fprintf(ptr_fichier_PSSM, "\nt  ");
//   for (t=0; t < taille_motif; t++)
//   {
//     fprintf(ptr_fichier_PSSM, "%.2f ", p_matrice_PSSM[1][t]);
//   }
//   fprintf(ptr_fichier_PSSM, "\nc  ");
//   for (c=0; c < taille_motif; c++)
//   {
//     fprintf(ptr_fichier_PSSM, "%.2f ", p_matrice_PSSM[2][c]);
//   }
//   fprintf(ptr_fichier_PSSM, "\ng  ");
//   for (g=0; g < taille_motif; g++)
//   {
//     fprintf(ptr_fichier_PSSM, "%.2f ", p_matrice_PSSM[3][g]);
//   }
//   fprintf(ptr_fichier_PSSM, "\n\n\n");
//   fclose(ptr_fichier_PSSM);
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// void calcul_score(TPtr_Mot_Ameliorer_PSSM* adr_mot, double*** adr_matrice_PSSM, int n_sequence, TPtr_dictionnaire_sequences* ptr_ensemble, int taille_motif)
// {
// 	TPtr_dictionnaire_sequences p_dictionnaire_seq = *ptr_ensemble;
// 	TPtr_Mot_Ameliorer_PSSM p_mot = *adr_mot;
// 	double** p_matrice_PSSM = *adr_matrice_PSSM;
// 	int position = 0;
// 	int i;
// 	double nb_a = 0;
// 	double nb_t = 0;
// 	double nb_c = 0;
// 	double nb_g = 0;
// 	double prob_mot_background = 0; //probabilite du mot selon le background P(M|B)
// 	double prob_mot_PSSM = 0; // probabilite du mot selon la PSSM P(M|PSSM)
// 	double score; // P(M|PSSM)/P(M|B)
// 	//CALCUL DU BACKGROUND:
// 	//On vérifie qu'on est dans la bonne séquence:
// 	if (p_dictionnaire_seq->numero_sequence !=  n_sequence)
// 	{
// 		p_dictionnaire_seq = p_dictionnaire_seq->suiv_seq;
// 	}
// 	while (position < 30)
// 	{
// 		switch(p_dictionnaire_seq->sequence[position])
// 		{
// 			case 'a':
// 				nb_a++;
// 				break;
// 			case 't':
// 				nb_t++;
// 				break;
// 			case 'c':
// 				nb_c++;
// 				break;
// 			case 'g':
// 				nb_g++;
// 		}
// 		position++;
// 	}
// 	nb_a = nb_a/30;
// 	nb_t = nb_t/30;
// 	nb_c = nb_c/30;
// 	nb_g = nb_g/30;
// 	position = 0;
//
// 	for (i=0; i <= taille_motif; i++)
// 	{
// 		switch(p_mot->mot[i])
// 		{
// 			case 'a':
// 				prob_mot_background = prob_mot_background + nb_a;
// 				prob_mot_PSSM = prob_mot_PSSM + p_matrice_PSSM[0][i];
// 				break;
// 			case 't':
// 				prob_mot_background = prob_mot_background + nb_t;
// 				prob_mot_PSSM = prob_mot_PSSM + p_matrice_PSSM[1][i];
// 				break;
// 			case 'c':
// 				prob_mot_background = prob_mot_background + nb_c;
// 				prob_mot_PSSM = prob_mot_PSSM + p_matrice_PSSM[2][i];
// 				break;
// 			case 'g':
// 				prob_mot_background = prob_mot_background + nb_g;
// 				prob_mot_PSSM = prob_mot_PSSM + p_matrice_PSSM[3][i];
// 				break;
// 		}
// 	}
// 	score = prob_mot_PSSM/prob_mot_background;
// 	p_mot->score_mot = score;
// 	return;
// }
// //------------------------------------------------------------------------------------------------------------
// double dist_PSSM(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv, double* distance_PSSM, int taille_motif)
// {
// 	double** p_ancienne_matrice_PSSM = *adr_matrice_PSSM;
// 	double** p_nouvelle_matrice_PSSM = *adr_matrice_PSSM_nouv;
// 	int i, j;
// 	double somme_distance = 0;
// 	for (i=0; i< 4; i++)
// 	{
// 		for (j=0; j<taille_motif; j++)
// 		{
// 			somme_distance = somme_distance + fabs(p_ancienne_matrice_PSSM[i][j] - p_nouvelle_matrice_PSSM[i][j]);
// 		}
// 	}
// 	*distance_PSSM = somme_distance;
// 	printf("Somme Distance = %f \n", somme_distance);
// 	return(somme_distance);
// }
// //------------------------------------------------------------------------------------------------------------
// int distanceHammingSt1(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st1, int taille_motif)
// {
// 	int i;
// 	int Dh = 0;
//   int st1 = 0;
// 	TPtr_Cell_Motif_PSSM p_mot_selected = *adr_mot_selected;
// 	Ptr_st p_st1 = *adr_st1;
//
// 	while (p_mot_selected != NULL)
// 	{
// 		for (i=0; i<taille_motif; i++)
// 		{
// 			if (*adr_Ct[i] != p_mot_selected->motif[i])
// 			{
// 				Dh++;
// 				printf("Dh: %d \n", Dh);
// 			}
// 		}
//     if (Dh > 1){  st1++;  }
//     strcpy(p_st1->mot, p_mot_selected->motif);
//     printf("p_st1->mot: %s \n", p_st1->mot);
//     p_st1->distance_hamming = Dh;
//     printf("p_st1->distance: %d \n", p_st1->distance_hamming);
//     Ptr_st p_st1_suiv = malloc(sizeof(st1));
//     p_st1->next_mot = p_st1_suiv;
//     p_st1 = p_st1_suiv;
// 		Dh = 0;
// 		p_mot_selected = p_mot_selected->suiv_motif;
// 	}
// 	return(st1-1); // Parce que le dernier mot NULL de la liste chainée est compté
// }
// //------------------------------------------------------------------------------------------------------------
// int distanceHammingSt2(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st2, int taille_motif)
// {
//   printf("DISTANCE HAMMING ST2");
//   int st2 = 0;
// 	int i;
// 	int Dh = 0;
// 	TPtr_Cell_Motif_PSSM p_mot_selected = *adr_mot_selected;
// 	Ptr_st p_st2 = *adr_st2;
// 	printf("IN ST2");
// 	while (p_mot_selected != NULL)
// 	{
// 		for (i=0; i<taille_motif; i++)
// 		{
// 			if (*adr_Ct[i] != p_mot_selected->motif[i])
// 			{
// 				Dh++;
// 			}
// 		}
//     strcpy(p_st2->mot, p_mot_selected->motif);
//     printf("p_st2->mot: %s \n", p_st2->mot);
//     p_st2->distance_hamming = Dh;
//     printf("p_st2->distance: %d \n", p_st2->distance_hamming);
//     Ptr_st p_st2_suiv = malloc(sizeof(st2));
//     p_st2->next_mot = p_st2_suiv;
//     p_st2 = p_st2_suiv;
// 		if(Dh < 2)
// 		{
//       st2++;
// 		}
// 		Dh = 0;
// 		p_mot_selected = p_mot_selected->suiv_motif;
// 	}
// 	return(st2);
// }
// //------------------------------------------------------------------------------------------------------------
// int distanceHammingSt2_prim(char (*adr_Ct)[6], TPtr_dictionnaire_sequences* ptr_ensemble, TPtr_Mot_Ameliorer_PSSM *adr_mot, Ptr_st* adr_st2_prim, int taille_motif)
// {
// 	printf("IN ST2 PRIM");
// 	int i;
//   int st2_prim = 0;
//   int position = 0;
// 	int Dh_min = 8;
// 	int Dh = 0;
// 	TPtr_dictionnaire_sequences p_dictionnaire_seq = *ptr_ensemble;
// 	Ptr_st p_st2_prim = *adr_st2_prim;
//   TPtr_Mot_Ameliorer_PSSM p_mot = *adr_mot;
//
// 	while (p_dictionnaire_seq != NULL)
// 	{
// 		while (position <= 25)
// 		{
// 			for (i=0; i<taille_motif; i++)
// 			{
//         p_mot->mot[i] = p_dictionnaire_seq->sequence[position];
// 				if (*adr_Ct[i] != p_mot->mot[i])
// 				{
// 					Dh++;
// 				}
//         if(i == 4)
// 				{
//           p_mot->mot[taille_motif] = '\0';
//         }
//         position++;
// 			}
//       position = position-4;
// 			if (Dh <= Dh_min)
// 			{
//         Dh_min = Dh;
// 				strcpy(p_st2_prim->mot, p_mot->mot);
// 				p_st2_prim->distance_hamming = Dh;
// 			}
// 			Dh = 0;
//       TPtr_Mot_Ameliorer_PSSM p_nouv_mot = malloc(sizeof(TMot_Ameliorer_PSSM));
//       p_mot->next_mot = p_nouv_mot;
//       p_mot = p_nouv_mot;
//     }
//     position = 0;
//     Dh_min = 8;
//     p_dictionnaire_seq = p_dictionnaire_seq->suiv_seq;
// 		printf("p_st2->mot (minimum): %s \n", p_st2_prim->mot);
// 		printf("p_st2->distance (minimum): %d \n", p_st2_prim->distance_hamming);
// 		Ptr_st p_st2_prim_suiv = malloc(sizeof(st));
// 		p_st2_prim->next_mot = p_st2_prim_suiv;
// 		p_st2_prim = p_st2_prim_suiv;
// 	}
//   p_st2_prim = *adr_st2_prim;
//   while (p_st2_prim != NULL)
//   {
//     if(p_st2_prim->distance_hamming < 2)
//     {
//       st2_prim++;
//     }
//     p_st2_prim = p_st2_prim->next_mot;
//   }
// 	return(st2_prim-1); //On met -1 car la boucle précedente prend en compte la case NULL de la fin de liste chainé
// }
// //------------------------------------------------------------------------------------------------------------
// void trier(int* v_St1_Dh, int* v_St1_Pos, int g, int d)
// {
//   int indice_pivot;
//
//   if (g < d)
//   {
//     separer(v_St1_Dh, v_St1_Pos, g, d, &indice_pivot);
//     indice_pivot --;
//     separer(v_St1_Dh, v_St1_Pos, g, d, &indice_pivot);
//     trier(v_St1_Dh, v_St1_Pos, indice_pivot+1, d);
//   }
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// void separer(int* v_St1_Dh, int* v_St1_Pos, int g, int d, int* adr_indice_pivot)
// {
//   //PRECONDITIONS
//   // g < d
//   int bas, haut; //indices de position dans le vecteur
//   int comp, pivot;  //comp pour comparateur, c'est variables sont des entiers car on trie un vecteur d'entiers
//   int comp_pos, pivot_pos;
//
//   bas = g;
//   haut = d;
//   pivot = v_St1_Dh[bas];
//   pivot_pos = v_St1_Pos[bas];
//   comp = v_St1_Dh[haut];
//   comp_pos = v_St1_Pos[haut];
//
//   while (bas < haut)
//   {
//     if ( comp > pivot)
//     {
//       v_St1_Dh[haut] = comp;
//       v_St1_Pos[haut] = comp_pos;
//       haut --;
//       comp = v_St1_Dh[haut];
//       comp_pos = v_St1_Pos[haut];
//     }else{
//       v_St1_Dh[bas] = comp;
//       v_St1_Pos[bas] = comp_pos;
//       bas ++;
//       comp = v_St1_Dh[bas];
//       comp_pos = v_St1_Pos[bas];
//     }
//   }
//   v_St1_Dh[bas] = pivot;
//   v_St1_Pos[bas] = pivot_pos;
//   *adr_indice_pivot = bas;
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// void quick_sort_ST(Ptr_st* adr_st1, int* v_St1_Pos, int nb_sequences)
// {
//   Ptr_st p_st1 = *adr_st1;
//   int i;
//   int v_St1_Dh[9]; //Distance de Hamming de chaque occurrence
//   int borne_gauche = 0;
//   int borne_droite = 9;
//   for (i=0; i<nb_sequences; i++)
//   {
//     v_St1_Dh[i] = p_st1->distance_hamming;
//     v_St1_Pos[i] = i;
//     printf("Vecteur Distance Hamming: %d \n", v_St1_Dh[i]);
//     printf("Vecteur Pos: %d \n", v_St1_Pos[i]);
//     p_st1 = p_st1->next_mot;
//   }
//   printf("\n");
//
//   trier(v_St1_Dh, v_St1_Pos, borne_gauche, borne_droite);
//
//   for(i=0; i<10; i++)
//   {
//     printf("Vecteur Distance Hamming après tri: %d \n", v_St1_Dh[i]);
//     printf("Vecteur Pos après tri: %d \n", v_St1_Pos[i]);
//   }
//  printf("\n");
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
// void fichier_sortie_st(Ptr_st* adr_st1, int* v_St1_Pos, char (*adr_Ct)[6], int nb_sequences)
// {
//   Ptr_st p_st1 = *adr_st1;
//   int i, position, cpt;
//   FILE* fichier_sortie_st2;
//   fichier_sortie_st2 = fopen("score_ST2.txt", "w");
//
//   for (i=0; i<nb_sequences; i++)
//   {
//     p_st1 = *adr_st1;
//     position = v_St1_Pos[i];
//     for (cpt=0; cpt < position; cpt++)
//     {
//       p_st1 = p_st1->next_mot;
//     }
//     fprintf(fichier_sortie_st2, "OCCURRENCE: %s --> SCORE: %d \n", p_st1->mot, p_st1->distance_hamming);
//   }
//   fprintf(fichier_sortie_st2, "\n\nMOTIF CONSENSUS: \n");
//   for (i=0; i<taille_motif; i++)
//   {
//     fprintf(fichier_sortie_st2, "%c", *adr_Ct[i]);
//   }
//   return;
// }
// //------------------------------------------------------------------------------------------------------------
