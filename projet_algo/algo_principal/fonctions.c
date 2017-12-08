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
void importer_parametres(char* nom_fichier, int* taille_motif, int* d, int* nb_fenetres, int* nb_masques)
{
	FILE* ptr_fichier; //creation d'un pointeur sur le fichier
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
void importer_sequences_fasta(char* nom_fichier_fasta, TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq )
{
  FILE* ptr_fichier_fasta;
  char contenu_ligne[128];
  char c; // utiliser pour lire le fichier caractère par caractère
  int cptr = 0; //compteur du nombre de seq dans le fichier

	// déclaration des pointeurs permettant de créer la liste chainée de séquences
  TPtr_dictionnaire_sequences p_new = NULL;
  TPtr_dictionnaire_sequences p = *adr_tete_dict_seq;
	TPtr_info_dictionnaire_sequences ptr_info = *adr_tete_info_dict_seq;
  ptr_info->nb_sequences = 0; //initialisation du nombre de seq contenues dans le fichier à 0

  ptr_fichier_fasta = fopen(nom_fichier_fasta, "r");
  if( ptr_fichier_fasta == NULL) { free(ptr_fichier_fasta); return; }// si le fichier est vide on sort de la fonction

	do //cette boucle permet de compter les séquences
	{
		c = fgetc (ptr_fichier_fasta); //lecture du fichier caractère par caractère
		if (c == '>') cptr++;  //incrémentation pour chaque nom de séquences
	} while (c != EOF); // l'action est répétée jusqu'à ce que la fin du fichier soit atteinte
	fclose(ptr_fichier_fasta);

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
void afficher_sequences(TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq )
{
	FILE* ptr_fichier;
	TPtr_dictionnaire_sequences p = *adr_tete_dict_seq;
	TPtr_info_dictionnaire_sequences ptr_info = *adr_tete_info_dict_seq;
	int cptr = 0;
	// printf("fonction afficher_sequences\n");
	// printf("%d \n", ptr_info->nb_sequences);

	ptr_fichier = fopen ( "verif_dico_fasta.txt" ,"w");
	while ( cptr < (ptr_info->nb_sequences))
	{
		// fprintf(ptr_fichier, "%d\n", p->numero_sequence);
		fputs( p->nom_seq, ptr_fichier);
		fputs( "\n" , ptr_fichier);
		fputs( p->sequence, ptr_fichier);
		fputs( "\n", ptr_fichier);
		// printf(" %s\n", p->nom_seq );
		// printf("%s\n", p->sequence);
		p = p->suiv_seq;
		cptr ++;
	}
	fclose(ptr_fichier);
}
//------------------------------------------------------------------------------------------------------------
void liberation_dictionnaire_sequence( TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq )
{
	TPtr_info_dictionnaire_sequences p_info = *adr_tete_info_dict_seq;
	TPtr_dictionnaire_sequences p_seq = *adr_tete_dict_seq;
	TPtr_dictionnaire_sequences p_prec_seq = NULL;

	while ( p_seq->suiv_seq != NULL)
	{
		p_prec_seq = p_seq;
		p_seq = p_seq->suiv_seq;
		free (p_prec_seq);
	}
	free(p_info);
}
//------------------------------------------------------------------------------------------------------------
void creation_PSSM(double*** adr_matrice_PSSM)
{
	int i, j;
  *adr_matrice_PSSM = malloc(4* sizeof(double*));

  for (i=0; i<4; i++) //initialisation de la matrice à 0:
  {
		(*adr_matrice_PSSM)[i] = malloc(taille_motif * sizeof(double)); //On alloue des tableaux de 'taille2' variables.
    for(j=0; j < taille_motif; j++)
    {
      (*adr_matrice_PSSM)[i][j] = 0;
    }
		printf("\n" );
  }
}
//------------------------------------------------------------------------------------------------------------
void afficher_PSSM( double** matrice_PSSM)
{
  FILE* ptr_fichier_PSSM;
  int i, j;
  ptr_fichier_PSSM = fopen("PSSM_Motif_Trouve.txt", "w");

  fprintf(ptr_fichier_PSSM, "\n-----------PSSM------------- \n");
	for (i=0;  i<4; i++){
		if(i == 0){ fprintf(ptr_fichier_PSSM, "A  ");}
		if(i == 1){ fprintf(ptr_fichier_PSSM, "\nT  ");}
		if(i == 2){ fprintf(ptr_fichier_PSSM, "\nC  ");}
		if(i == 3){ fprintf(ptr_fichier_PSSM, "\nG  ");}
		for(j=0; j<taille_motif; j++){
		  fprintf(ptr_fichier_PSSM, "%.2f ", matrice_PSSM[i][j]);
		}
	}
  fprintf(ptr_fichier_PSSM, "\n----------------------------\n");
  fclose(ptr_fichier_PSSM);
}
//------------------------------------------------------------------------------------------------------------
void liberation_PSSM(double*** adr_matrice_PSSM)
{
	int i;
	double** p_PSSM = *adr_matrice_PSSM;

  for (i = 0; i<4; i++)
	{
		free(p_PSSM[i]);
	}
  free(p_PSSM);
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
int* generation_masque(int* masque)
{
  int i;
  // int p_masque = *adr_masque;
  int nb_fenetres_ouvertes;

  while (nb_fenetres_ouvertes != nb_fenetres)
  {
		nb_fenetres_ouvertes=0;
    for (i=0; i <= (taille_motif-1); i++)
    {
      (masque)[i] = random_number(2,0);
      if ((masque)[i] == 1) {  nb_fenetres_ouvertes ++;  }
    }
  }
	return(masque);
}
//------------------------------------------------------------------------------------------------------------
void generation_dictionnaire_kmer(int position_kmer, char* k_mer, char* motif, TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
{
	int i = 1;
	int j = 1;
	int is_kmer = 0; //est ce que la brique de ce kmer existe déjà?
	int is_seq = 0; //est ce que la brique de cette séquence existe déjà?

	TPtr_dictionnaire_sequences p_dictionnaire_sequence = tete_dict_seq;
  TPtr_Cellkmer p_kmer = tete_info_dict_kmer->tete_liste_kmer; //on créer des pointeurs temporaires pour parcourir les listes
	TPtr_Cellkmer p_new_kmer = p_kmer;
	TPtr_CellSequence p_Kseq ;
	TPtr_CellSequence p_new_Kseq ;
	TPtr_CellPos p_Kpos ;
	TPtr_CellPos p_new_Kpos ;

  if (tete_info_dict_kmer->nb_kmer == 0) //cas du premier element du dictionnaire de kmer
  {
		//construction bloc Kmer
		p_kmer = malloc(sizeof(TCellkmer));
		tete_info_dict_kmer->tete_liste_kmer = p_kmer;
		p_kmer->kmer = (char*) calloc(nb_fenetres+1, sizeof(char));
    strcpy(p_kmer->kmer, k_mer);
		tete_info_dict_kmer->nb_kmer = 1;
    p_kmer->suiv_kmer = NULL;
 		//construction du bloc de séquence
		p_Kseq = malloc(sizeof(TCellSequence));
		p_kmer->tete_sequence = p_Kseq; //connection bloc kmer et bloc seq
    p_Kseq->num_sequence = p_dictionnaire_sequence->numero_sequence;
    p_Kseq->suiv_sequence = NULL;
		p_kmer->nb_sequence = 1;
 		//construction du bloc de position
		p_Kpos = malloc(sizeof(TCellPos));
    p_Kseq->tete_pos = p_Kpos; //connection bloc seq et bloc position
    p_Kpos->position = position_kmer;
		p_Kpos->motif = (char*) calloc((taille_motif+1), sizeof(char));
		strcpy(p_Kpos->motif, motif);
    p_Kpos->suiv_pos = NULL;
    return;
  }
	else
	{
	  while( i <= tete_info_dict_kmer->nb_kmer ) //on parcourt les blocs de kmer
	  {
	    if ( strncmp(p_kmer->kmer, k_mer, nb_fenetres) == 0)   //si le kmer a déjà été trouvé ie. le bloc existe donc il existe au moins un bloc de seq et un bloc de pos
	    { //on est positionné sur le bon bloc Kmer
				is_kmer = 1;

				p_Kseq = p_kmer->tete_sequence;
	      while ( j <= p_kmer->nb_sequence ) //on parcourt les blocs de séquence
	      {
	        if (p_dictionnaire_sequence->numero_sequence == p_Kseq->num_sequence) //Si le kmer a déjà été trouvé dans la séquence et il y a plusieur bloc de seq
	        {	//on est positionné sur le bon bloc seq
						is_seq = 1;

	          p_Kpos = p_Kseq->tete_pos;
	          while (p_Kpos != NULL)
						{  //on parcourt les blocs de position jusqu'aiu bout de la liste chainée
	            p_Kpos = p_Kpos->suiv_pos;
	          }

	          p_new_Kpos = (TCellPos*)malloc(sizeof(TCellPos)); // On créer une nouvelle brique de position
						p_new_Kpos->position = position_kmer;
						p_new_Kpos->motif = (char*) calloc((taille_motif+1), sizeof(char));
						strcpy(p_new_Kpos->motif, motif);
						p_new_Kpos->suiv_pos = NULL;
						p_Kpos = p_new_Kpos; //chainage
						// return;
	        }
					if ( p_Kseq->suiv_sequence == NULL){ break; } else { p_Kseq = p_Kseq->suiv_sequence;  j++;}
	      }

				if (is_seq == 0) // pas d'occurrence de ce kmer dans la séquence
				{
					p_new_Kseq = (TCellSequence*)malloc(sizeof(TCellSequence)); // On créer une nouvelle brique de séquence
					p_Kseq->suiv_sequence = p_new_Kseq; //chainage
					p_new_Kseq->num_sequence = p_dictionnaire_sequence->numero_sequence;
					p_new_Kseq->suiv_sequence = NULL;
					p_kmer->nb_sequence = p_kmer->nb_sequence + 1;

					p_new_Kpos = (TCellPos*)malloc(sizeof(TCellPos)); // On créer une nouvelle brique de position
					p_new_Kseq->tete_pos = p_new_Kpos; //chainage
					p_new_Kpos->position = position_kmer;
					p_new_Kpos->motif = (char*) calloc((taille_motif+1), sizeof(char));
					strcpy(p_new_Kpos->motif, motif);
					p_new_Kpos->suiv_pos = NULL;
				}
	    }
			if (p_kmer->suiv_kmer == NULL){ break; } else { p_kmer = p_kmer->suiv_kmer;  i++;}
	  }

		if (is_kmer == 0)
		{
			// Si ce kmer n'a pas encore été trouvé ie.creation du bloc
			p_new_kmer = (TCellkmer*)malloc(sizeof(TCellkmer)); // On créer une nouvelle brique de kmer
			p_kmer->suiv_kmer = p_new_kmer; //chainage
			tete_info_dict_kmer->nb_kmer = tete_info_dict_kmer->nb_kmer + 1;
			p_new_kmer->kmer = (char*) calloc(nb_fenetres+1, sizeof(char));
			strcpy(p_new_kmer->kmer, k_mer);
			p_new_kmer->suiv_kmer = NULL;

			p_new_Kseq = (TCellSequence*)malloc(sizeof(TCellSequence)); // On créer une nouvelle brique de séquence
			p_new_kmer->tete_sequence = p_new_Kseq; //chainage
			p_new_Kseq->num_sequence = p_dictionnaire_sequence->numero_sequence;
			p_new_Kseq->suiv_sequence = NULL;
			p_new_kmer->nb_sequence = p_new_kmer->nb_sequence + 1;

			p_new_Kpos = (TCellPos*)malloc(sizeof(TCellPos)); // On créer une nouvelle brique de position
			p_new_Kseq->tete_pos = p_new_Kpos; //chainage
			p_new_Kpos->position = position_kmer;
			p_new_Kpos->motif = (char*) calloc(taille_motif+1, sizeof(char));
			strcpy(p_new_Kpos->motif, motif);
			p_new_Kpos->suiv_pos = NULL;
		}
	}
}
//------------------------------------------------------------------------------------------------------------
void parcours_masque_sur_seq( int* masque, TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
{
	TPtr_dictionnaire_sequences p_dictionnaire_sequence = tete_dict_seq;

  int position = 0;
  int i, position_kmer;
  int pos_kmer = 0; // correspond au nombre de fenetre ouvertes ie. la taille du kmer
  char* k_mer = (char*) calloc((nb_fenetres+1),sizeof(char));  // le k_mer mesure la taille du nombre de fenêtres ouvertes dans le masque,  !!! +1 pour le caractère de fin de chaine
	char* motif = (char*) calloc((taille_motif+1),sizeof(char));

  while (p_dictionnaire_sequence != NULL) //parcourt du dictionnaire de séquence
  {
    while (position < 25) //inferieur à longueur de la séquence
    {
      for(i = 0; i < taille_motif; i++) // on parcourt les nucléotides sous le masque
      {
        if ( masque[i] == 1) //si la fenêtre du masque est ouverte
        {
          k_mer[pos_kmer] = p_dictionnaire_sequence->sequence[position];
          if (pos_kmer == 0) {
						position_kmer = position;
					}
          pos_kmer ++;
          if (pos_kmer == nb_fenetres)
          {
            k_mer[pos_kmer] = '\0' ;
            pos_kmer = 0;
          }
        }
				position ++;
				motif[i] = p_dictionnaire_sequence->sequence[position]; // on récupère le motif correspondant au kmer
      }
			motif[taille_motif] = '\0'; // caractère de fin de chaîne
      generation_dictionnaire_kmer(position_kmer, k_mer, motif, p_dictionnaire_sequence , tete_info_dict_kmer);
			creation_liste_motifs(motif, tete_info_dict_kmer);
			position = position - taille_motif +1;
    }
    position = 0;
    p_dictionnaire_sequence = p_dictionnaire_sequence->suiv_seq;
  }
	free(k_mer);
	free(motif);
}
//------------------------------------------------------------------------------------------------------------
void affichage_dictionnaire_kmer(TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
{
  TPtr_Cellkmer p_kmer = tete_info_dict_kmer->tete_liste_kmer;
  TPtr_CellSequence p_sequence = p_kmer->tete_sequence;
	TPtr_CellPos p_pos = p_sequence->tete_pos;

  FILE* fichier_dictionnaire;
  fichier_dictionnaire = fopen("dictionnaire_kmer.txt", "w");

	while( p_kmer != NULL ) //on parcourt les blocs de kmer
  {
    fprintf(fichier_dictionnaire, "---KMER--- %s \n", p_kmer->kmer);
		fprintf(fichier_dictionnaire, "Nb de séquences: %2.f \n", p_kmer->nb_sequence);

		p_sequence = p_kmer->tete_sequence;
    while (p_sequence != NULL) // on parcourt la liste de séquences par kmer
    {
			fprintf(fichier_dictionnaire, "-----------\n");
      fprintf(fichier_dictionnaire, "Sequence_num:-> %d \n", p_sequence->num_sequence);

			p_pos = p_sequence->tete_pos;
      while (p_pos != NULL) //on parcout la liste de positions par séquence
      {
        fprintf(fichier_dictionnaire, "Positions:  -> %d \n", p_pos->position);
				fprintf(fichier_dictionnaire, "            -> %s \n", p_pos->motif);

        p_pos = p_pos->suiv_pos;
      }
      p_sequence = p_sequence->suiv_sequence;
    }

    fprintf(fichier_dictionnaire, "-------------------------------------------\n\n");
		p_kmer = p_kmer->suiv_kmer;
  }
  fclose(fichier_dictionnaire);
}
// //------------------------------------------------------------------------------------------------------------
void liberation_dictionnaire_kmer(TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
{
	TPtr_Cellkmer p_kmer = tete_info_dict_kmer->tete_liste_kmer;
	TPtr_Cellkmer p_prec_kmer = NULL;
	TPtr_CellSequence p_seq = p_kmer->tete_sequence;
	TPtr_CellSequence p_prec_seq = NULL;
	TPtr_CellPos p_pos = p_seq->tete_pos;
	TPtr_CellPos p_prec_pos = NULL;

	while( p_kmer != NULL ) //on parcourt les blocs de kmer
  {
		p_seq = p_kmer->tete_sequence;
    while (p_seq != NULL) // on parcourt la liste de séquences par kmer
    {
			p_pos = p_seq->tete_pos;
      while (p_pos != NULL) //on parcout la liste de positions par séquence
      {
				p_prec_pos = p_pos;
        p_pos = p_pos->suiv_pos;
				free(p_prec_pos->motif);
				free(p_prec_pos);
      }
			p_prec_seq = p_seq;
      p_seq = p_seq->suiv_sequence;
			free(p_prec_seq);
    }
		p_prec_kmer = p_kmer;
		p_kmer = p_kmer->suiv_kmer;
		free(p_prec_kmer->kmer);
		free(p_prec_kmer);
  }
}
//------------------------------------------------------------------------------------------------------------
void creation_liste_motifs(char* motif ,TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
{
	TPtr_Cell_Motif p = tete_info_dict_kmer->tete_liste_motif;
	TPtr_Cell_Motif p_new = NULL;

	if (tete_info_dict_kmer->tete_liste_motif == NULL)
	{
		tete_info_dict_kmer->tete_liste_motif = (TCell_Motif*)malloc(sizeof(TCell_Motif));
		p = tete_info_dict_kmer->tete_liste_motif;
	}
	while ( p != NULL)
	{
		p_new = (TCell_Motif*)malloc(sizeof(TCell_Motif));
		p_new->suiv_motif = NULL;
		p_new->motif = (char*) calloc(taille_motif+1, sizeof(char));
		stcpy(p_new->motif, motif);
		p->suiv_motif = p_new;

		p = p_new;
	}
}
//------------------------------------------------------------------------------------------------------------
void afficher_liste_motifs(TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
{
	TPtr_Cell_Motif p = tete_info_dict_kmer->tete_liste_motif;
	TPtr_Cell_Motif p_prec = NULL;

	FILE* fichier_motifs;
	fichier_motifs = fopen("dictionnaire_motifs.txt", "w");

	fprintf(fichier_motifs, "MOTIFS POUR CALCUL DE LA PSSM\n" );
	while( p != NULL)
	{
		fprintf(fichier_motifs, "%s\n", p->motif);
		p_prec = p;
		p = p->suiv_motifs;
	}
	fclose(fichier_motifs);
}
//------------------------------------------------------------------------------------------------------------
void liberation_liste_motifs(TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
{
	TPtr_Cell_Motif p = tete_info_dict_kmer->tete_liste_motif;
	TPtr_Cell_Motif p_prec = NULL;

	while( p != NULL)
	{
		p_prec = p;
		p = p->suiv_motif;
		free(p_prec);
	}
}
//------------------------------------------------------------------------------------------------------------
// void calcul_PSSM(TPtr_Cellkmer_selectionne *adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM *adr_cell_motif_PSSM, double*** adr_matrice_PSSM)
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
// }
//------------------------------------------------------------------------------------------------------------
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
// }
// //------------------------------------------------------------------------------------------------------------
