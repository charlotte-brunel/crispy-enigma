#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include "fonctions.h"
//==================================================================================================
//==================================================================================================
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
//==================================================================================================
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
//==================================================================================================
void affichage_dictionnaire_sequences(TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq )
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
//==================================================================================================
void liberation_dictionnaire_sequences( TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq )
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
//==================================================================================================
//==================================================================================================
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
//==================================================================================================
void initialisation_masque(int** adr_masque)
{
	int i;
  int* masque = *adr_masque;

	for (i=0; i<taille_motif; i++)
	{
		masque[i] = 0;
	}
}
//==================================================================================================
void affichage_masque(int* masque)
{
	int i;
	for (i=0; i<taille_motif; i++)
	{
		printf("%d", masque[i]);
	}
	printf("\n");
}
//==================================================================================================
void generation_masque(int** adr_masque)
{
  int i;
  int* masque = *adr_masque;
  int nb_fenetres_ouvertes;

  while (nb_fenetres_ouvertes != nb_fenetres)
  {
		nb_fenetres_ouvertes=0;
    for (i=0; i < taille_motif; i++)
    {
      (masque)[i] = random_number(2,0);
      if ((masque)[i] == 1) {  nb_fenetres_ouvertes ++;  }
    }
  }
}
//==================================================================================================
void generation_dictionnaire_kmers(int position_kmer, char* k_mer, char* motif, TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
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
//==================================================================================================
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
      generation_dictionnaire_kmers(position_kmer, k_mer, motif, p_dictionnaire_sequence , tete_info_dict_kmer);
			position = position - taille_motif +1;
    }
    position = 0;
    p_dictionnaire_sequence = p_dictionnaire_sequence->suiv_seq;
  }
	free(k_mer);
	free(motif);
}
//==================================================================================================
void affichage_dictionnaire_kmers(TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
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
//==================================================================================================
void liberation_dictionnaire_kmers(TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
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
		free(p_prec_kmer->motif_consensus);
		free(p_prec_kmer->kmer);
		liberation_PSSM(&(p_prec_kmer->PSSM_consensus));
		free(p_prec_kmer);
  }
}
//==================================================================================================
//==================================================================================================
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
  }
}
//==================================================================================================
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
//==================================================================================================
double calcul_distance_PSSMs(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv)
{
	double** p_ancienne_matrice_PSSM = *adr_matrice_PSSM;
	double** p_nouvelle_matrice_PSSM = *adr_matrice_PSSM_nouv;
	int i, j;
	double distance = 0;

	for (i=0; i< 4; i++)
	{
		for (j=0; j<taille_motif; j++)
		{
			distance = distance + fabs(p_ancienne_matrice_PSSM[i][j] - p_nouvelle_matrice_PSSM[i][j]);
		}
	}
	return(distance);
}
//==================================================================================================
void copie_PSSM(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv)
{
	double** p_matrice_PSSM = *adr_matrice_PSSM;
	double** p_nouvelle_matrice_PSSM = *adr_matrice_PSSM_nouv;
	int i, j;

	for (i=0; i< 4; i++)
	{
		for (j=0; j<taille_motif; j++)
		{
			p_matrice_PSSM[i][j] = p_nouvelle_matrice_PSSM[i][j];
		}
	}
}
//==================================================================================================
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
//==================================================================================================
void calcul_PSSM(TPtr_Cellkmer p_kmer , double*** adr_matrice_PSSM )
{
	TPtr_CellSequence p_seq ;
	TPtr_CellPos p_pos ;

  double** p_matrice_PSSM = *adr_matrice_PSSM;
  int i;

	p_seq = p_kmer->tete_sequence;
  while (p_seq != NULL) // on parcourt la liste de séquences par kmer
  {
		p_pos = p_seq->tete_pos;
		for (i=0; i < taille_motif; i++) //on parcourt chaque nucléotide du motif
    {
      switch(p_pos->motif[i] )
      {
        case 'a':
            p_matrice_PSSM[0][i] = p_matrice_PSSM[0][i] + (1/p_kmer->nb_sequence);
            break;
        case 't':
            p_matrice_PSSM[1][i] = p_matrice_PSSM[1][i] + (1/p_kmer->nb_sequence);
            break;
        case 'c':
            p_matrice_PSSM[2][i] = p_matrice_PSSM[2][i] + (1/p_kmer->nb_sequence);
            break;
        case 'g':
            p_matrice_PSSM[3][i] = p_matrice_PSSM[3][i] + (1/p_kmer->nb_sequence);
            break;
      }
    }
  p_seq = p_seq->suiv_sequence;
  }
}
//==================================================================================================
void calcul_PSSM_amelioree(TPtr_Cellkmer p_kmer , double*** adr_matrice_PSSM )
{
	TPtr_CellSequence p_seq ;
	TPtr_CellPos p_pos ;

	double** p_matrice_PSSM = *adr_matrice_PSSM;
	int i;

	p_seq = p_kmer->tete_sequence;
	while (p_seq != NULL) // on parcourt la liste de séquences par kmer
	{
		p_pos = p_seq->tete_pos_max;
		for (i=0; i < taille_motif; i++) //on parcourt chaque nucléotide du motif
		{
			switch(p_pos->motif[i] )
			{
				case 'a':
						p_matrice_PSSM[0][i] = p_matrice_PSSM[0][i] + (1/p_kmer->nb_sequence);
						break;
				case 't':
						p_matrice_PSSM[1][i] = p_matrice_PSSM[1][i] + (1/p_kmer->nb_sequence);
						break;
				case 'c':
						p_matrice_PSSM[2][i] = p_matrice_PSSM[2][i] + (1/p_kmer->nb_sequence);
						break;
				case 'g':
						p_matrice_PSSM[3][i] = p_matrice_PSSM[3][i] + (1/p_kmer->nb_sequence);
						break;
			}
		}
	p_seq = p_seq->suiv_sequence;
	}
}
//==================================================================================================
char* identification_motif_consensus(double** matrice_PSSM)
{
	int i, j;
	int maximum = -1;
	char* mot_consensus = (char*) calloc(taille_motif+1, sizeof(char));

	for(j=0; j<taille_motif; j++)
	{
		for (i=0;  i<4; i++)
		{
			if (matrice_PSSM[i][j] > maximum)
			{
				maximum = matrice_PSSM[i][j];
				switch (i)
				{
					case 0:
							mot_consensus[i] = 'a';
							break;
					case 1:
							mot_consensus[i] = 't';
							break;
					case 2:
							mot_consensus[i] = 'c';
							break;
					case 3:
							mot_consensus[i] = 'g';
							break;
				}
			}
		}
		maximum = -1; //on remet maximum a -1 avant d'�valuer la nucl�otide majoritaire de la s�quence suivante !
	}
	mot_consensus[taille_motif+1] = '\0';
	return(mot_consensus);
}
//==================================================================================================
void calculer_score(TPtr_CellPos* adr_p_pos, double** matrice_PSSM, TPtr_info_dictionnaire_sequences tete_info_dict_seq)
{
	TPtr_info_dictionnaire_sequences p_info_dict_seq = tete_info_dict_seq;
	TPtr_dictionnaire_sequences p_dictionnaire_seq = p_info_dict_seq->tete_dictionnaire_seq;
	TPtr_CellPos p_pos = *adr_p_pos;

	int i;
	int position = 0;
	double background[4] = {0};
	double Proba_motif_background = 0; //probabilite du mot selon le background P(M|B)
	double Proba_motif_PSSM = 0; // probabilite du mot selon la PSSM P(M|PSSM)
	double score; // P(M|PSSM)/P(M|B)

	//CALCUL DU BACKGROUND
	while (p_dictionnaire_seq != NULL)
	{
		while (position < 30) // idéalement < longueur séquence
		{
			switch( p_dictionnaire_seq->sequence[position] )
			{
				case 'a':
					background[0] = background[0] + (1/p_info_dict_seq->nb_sequences);
					break;
				case 't':
					background[1] = background[1] + (1/p_info_dict_seq->nb_sequences);
					break;
				case 'c':
					background[2] = background[2] + (1/p_info_dict_seq->nb_sequences);
					break;
				case 'g':
					background[3] = background[3] + (1/p_info_dict_seq->nb_sequences);
			}
			position++;
		}
		p_dictionnaire_seq = p_dictionnaire_seq->suiv_seq;
	}

	for (i=0; i <= taille_motif; i++)
	{
		switch(p_pos->motif[i])
		{
			case 'a':
				Proba_motif_background = Proba_motif_background * background[0];
				Proba_motif_PSSM = Proba_motif_PSSM * matrice_PSSM[0][i];
				break;
			case 't':
				Proba_motif_background = Proba_motif_background * background[1];
				Proba_motif_PSSM = Proba_motif_PSSM * matrice_PSSM[1][i];
				break;
			case 'c':
				Proba_motif_background = Proba_motif_background * background[2];
				Proba_motif_PSSM = Proba_motif_PSSM * matrice_PSSM[2][i];
				break;
			case 'g':
				Proba_motif_background = Proba_motif_background * background[3];
				Proba_motif_PSSM = Proba_motif_PSSM * matrice_PSSM[3][i];
				break;
		}
	}
	score = Proba_motif_PSSM/Proba_motif_background;
	p_pos->score = score;
}
//==================================================================================================
void recherche_motifs_maximisant_scores(TPtr_Cellkmer p_kmer, double** matrice_PSSM, TPtr_info_dictionnaire_sequences tete_info_dict_seq)
{
	double score_max = INT_MIN ;
  int pos_max = -1;

	TPtr_CellSequence p_seq = p_kmer->tete_sequence;
	TPtr_CellPos p_pos = p_seq->tete_pos;

	p_seq = p_kmer->tete_sequence;
  while (p_seq != NULL) // on parcourt la liste de séquences par kmer
  {
		pos_max = -1;
		score_max = INT_MIN;
		p_pos = p_seq->tete_pos;
    while (p_pos != NULL) //on parcout la liste de positions par séquence
    {
			calculer_score(&p_pos, matrice_PSSM, tete_info_dict_seq);
			if (p_pos->score > score_max)
			{
				score_max = p_pos->score;
				pos_max = p_pos->position;
				p_seq->tete_pos_max = p_pos;
			}
      p_pos = p_pos->suiv_pos;
    }
    p_seq = p_seq->suiv_sequence;
  }
	//Post condition: le mot de la sequence courante le plus proche de la PSSM est identifié !
}
//==================================================================================================
int distance_Hamming_St1(TPtr_Cellkmer p_kmer)
{
	TPtr_CellSequence p_seq = p_kmer->tete_sequence;
	TPtr_CellPos p_pos = p_seq->tete_pos;

	int i;
	int distance_hamming = 0;
  int st1 = 0;

	p_seq = p_kmer->tete_sequence;
	while (p_seq != NULL) // on parcourt la liste de séquences par kmer
	{
		p_pos = p_seq->tete_pos_max;
		for (i=0; i<taille_motif; i++)
		{
			if ( p_pos->motif[i] != p_kmer->motif_consensus[i] )
			{
				distance_hamming ++;
			}
		}
		if ( distance_hamming > d) { st1 ++;}

		p_seq = p_seq->suiv_sequence;
		distance_hamming = 0;
	}
	return(st1);
}
//==================================================================================================
int distance_Hamming_St2_T(TPtr_Cellkmer p_kmer)
{
	TPtr_CellSequence p_seq = p_kmer->tete_sequence;
	TPtr_CellPos p_pos = p_seq->tete_pos;

	int i;
	int distance_hamming = 0;
  int st2 = 0;

	p_seq = p_kmer->tete_sequence;
	while (p_seq != NULL) // on parcourt la liste de séquences par kmer
	{
		p_pos = p_seq->tete_pos_max;
		for (i=0; i<taille_motif; i++)
		{
			if ( p_pos->motif[i] != p_kmer->motif_consensus[i] )
			{
				distance_hamming ++;
			}
		}
		if ( distance_hamming <= d) { st2 ++;}
		p_seq = p_seq->suiv_sequence;
		distance_hamming = 0;
	}
	return(st2);
}
//==================================================================================================
int distance_Hamming_St2_T_prim(TPtr_Cellkmer p_kmer)
{
	TPtr_CellSequence p_seq = p_kmer->tete_sequence;
	TPtr_CellPos p_pos = p_seq->tete_pos;

	int i;
	int distance_hamming = 0;
  int st2_prim = 0;

	p_seq = p_kmer->tete_sequence;
	while (p_seq != NULL) // on parcourt la liste de séquences par kmer
	{
		p_pos = p_seq->tete_pos_dH_min;
		for (i=0; i<taille_motif; i++)
		{
			if ( p_pos->motif[i] != p_kmer->motif_consensus[i] )
			{
				distance_hamming ++;
			}
		}
		if ( distance_hamming <= d) { st2_prim ++;}
		p_seq = p_seq->suiv_sequence;
		distance_hamming = 0;
	}
	return(st2_prim);
}
//==================================================================================================
void recherche_motifs_minimisant_dHamming(TPtr_Cellkmer p_kmer)
{
	TPtr_CellSequence p_seq = p_kmer->tete_sequence;
	TPtr_CellPos p_pos = p_seq->tete_pos;

	int i;
	int dH_min = INT_MAX;
	int distance_hamming = 0;

	p_seq = p_kmer->tete_sequence;
	while (p_seq != NULL) // on parcourt la liste de séquences par kmer
	{
		dH_min = INT_MAX;

		p_pos = p_seq->tete_pos;
		p_seq->tete_pos_dH_min = p_pos;
		while (p_pos != NULL) //on parcout la liste de positions par séquence
    {
			for (i=0; i<taille_motif; i++)
			{
				if ( p_pos->motif[i] != p_kmer->motif_consensus[i] )
				{
					distance_hamming ++;
				}
			}
			if ( distance_hamming < dH_min)
			{
				dH_min = distance_hamming;
				p_seq->tete_pos_dH_min = p_pos;
			}
      p_pos = p_pos->suiv_pos;
    }
		p_seq = p_seq->suiv_sequence;
		distance_hamming = 0;
	}
}
//==================================================================================================
void raffiner_version1(TPtr_Cellkmer p_kmer)
{
	p_kmer->score_St1 = distance_Hamming_St1( p_kmer );
}
//==================================================================================================
void egalisation_T(TPtr_Cellkmer p_kmer)
{
	TPtr_CellSequence p_seq = p_kmer->tete_sequence;

	while (p_seq != NULL) // on parcourt la liste de séquences par kmer
	{
		p_seq->tete_pos_max = p_seq->tete_pos_max;

		p_seq = p_seq->suiv_sequence;
	}
}
//==================================================================================================
int verification_convergence_T(TPtr_Cellkmer p_kmer)
{
	TPtr_CellSequence p_seq = p_kmer->tete_sequence;
	int cpt = 0;

	while (p_seq != NULL) // on parcourt la liste de séquences par kmer
	{
		if (p_seq->tete_pos_max == p_seq->tete_pos_max) { cpt ++;}

		p_seq = p_seq->suiv_sequence;
	}
	if (cpt == p_kmer->nb_sequence){ return(1); }else{ return(0); }
}
//==================================================================================================
void raffiner_version2(TPtr_Cellkmer p_kmer)
{
	int score_St2 = 0;
	int score_St2_prim = 0;
	int convergence_T = false;
	double** PSSM;
	char* motif_consensus_prim = (char*) calloc((taille_motif+1), sizeof(char));

	do {
		score_St2 = distance_Hamming_St2_T( p_kmer );
		recherche_motifs_minimisant_dHamming( p_kmer); // définition de l'ensemble T'

		score_St2_prim = distance_Hamming_St2_T_prim(p_kmer);

		if (score_St2_prim > score_St2)
		{
			egalisation_T(p_kmer); // on égalise T à T'
			creation_PSSM(&PSSM);
			calcul_PSSM_amelioree(p_kmer, &PSSM);
			motif_consensus_prim = identification_motif_consensus(PSSM);
			p_kmer->motif_consensus = motif_consensus_prim;
			p_kmer->PSSM_consensus = PSSM;
		}
	} while( verification_convergence_T(p_kmer) == 0);

}
//==================================================================================================
void trier_ST1(TPtr_Cellkmer* (*adr_tableau_kmer_QS), int g, int d)
{
  int indice_pivot;

  if (g < d)
  {
    separer_ST1(adr_tableau_kmer_QS, g, d, &indice_pivot);
    trier_ST1(adr_tableau_kmer_QS, g, indice_pivot-1);
    trier_ST1(adr_tableau_kmer_QS, indice_pivot+1, d);
  }
}
//==================================================================================================
void separer_ST1(TPtr_Cellkmer* (*adr_tableau_kmer_QS), int g, int d, int* adr_indice_pivot)
{
  //PRECONDITIONS
  // g < d
  int bas, haut; //indices de position dans le vecteur
  int comp, pivot;  //comp pour comparateur, c'est variables sont des entiers car on trie un vecteur d'entiers
  int comp_pos, pivot_pos;
	TPtr_Cellkmer* p_tableau_kmer_QS = *adr_tableau_kmer_QS;

  bas = g;
  haut = d;
  pivot = p_tableau_kmer_QS[bas]->score_St1;
  pivot_pos = p_tableau_kmer_QS[bas]->score_St1;
  comp = p_tableau_kmer_QS[haut]->score_St1;
  comp_pos = p_tableau_kmer_QS[haut]->score_St1;

  while (bas < haut)
  {
    if ( comp > pivot)
    {
      p_tableau_kmer_QS[haut]->score_St1 = comp;
      p_tableau_kmer_QS[haut]->score_St1 = comp_pos;
      haut --;
      comp = p_tableau_kmer_QS[haut]->score_St1;
      comp_pos = p_tableau_kmer_QS[haut]->score_St1;
    }else{
      p_tableau_kmer_QS[bas]->score_St1 = comp;
      p_tableau_kmer_QS[bas]->score_St1 = comp_pos;
      bas ++;
      comp = p_tableau_kmer_QS[bas]->score_St1;
      comp_pos = p_tableau_kmer_QS[bas]->score_St1;
    }
  }
  p_tableau_kmer_QS[bas]->score_St1 = pivot;
  p_tableau_kmer_QS[bas]->score_St1 = pivot_pos;
  *adr_indice_pivot = bas;
}
//==================================================================================================
void trier_ST2(TPtr_Cellkmer* (*adr_tableau_kmer_QS), int g, int d)
{
  int indice_pivot;

  if (g < d)
  {
    separer_ST2(adr_tableau_kmer_QS, g, d, &indice_pivot);
    trier_ST2(adr_tableau_kmer_QS, g, indice_pivot-1);
    trier_ST2(adr_tableau_kmer_QS, indice_pivot+1, d);
  }
}
//==================================================================================================
void separer_ST2(TPtr_Cellkmer* (*adr_tableau_kmer_QS), int g, int d, int* adr_indice_pivot)
{
  //PRECONDITIONS
  // g < d
  int bas, haut; //indices de position dans le vecteur
  int comp, pivot;  //comp pour comparateur, c'est variables sont des entiers car on trie un vecteur d'entiers
  int comp_pos, pivot_pos;
	TPtr_Cellkmer* p_tableau_kmer_QS = *adr_tableau_kmer_QS;

  bas = g;
  haut = d;
  pivot = p_tableau_kmer_QS[bas]->score_St2;
  pivot_pos = p_tableau_kmer_QS[bas]->score_St2;
  comp = p_tableau_kmer_QS[haut]->score_St2;
  comp_pos = p_tableau_kmer_QS[haut]->score_St2;

  while (bas < haut)
  {
    if ( comp > pivot)
    {
      p_tableau_kmer_QS[haut]->score_St2 = comp;
      p_tableau_kmer_QS[haut]->score_St2 = comp_pos;
      haut --;
      comp = p_tableau_kmer_QS[haut]->score_St2;
      comp_pos = p_tableau_kmer_QS[haut]->score_St2;
    }else{
      p_tableau_kmer_QS[bas]->score_St2 = comp;
      p_tableau_kmer_QS[bas]->score_St2 = comp_pos;
      bas ++;
      comp = p_tableau_kmer_QS[bas]->score_St2;
      comp_pos = p_tableau_kmer_QS[bas]->score_St2;
    }
  }
  p_tableau_kmer_QS[bas]->score_St2 = pivot;
  p_tableau_kmer_QS[bas]->score_St2 = pivot_pos;
  *adr_indice_pivot = bas;
}
//==================================================================================================
void quick_sort_ST(TPtr_Cellkmer* (*adr_tableau_kmer_QS), TPtr_info_dictionnaire_kmer tete_info_dict_kmer)
{
	TPtr_Cellkmer* p_tableau_kmer_QS = *adr_tableau_kmer_QS;
	TPtr_Cellkmer p_kmer = tete_info_dict_kmer->tete_liste_kmer;

	int i = 0;
	int borne_gauche = 0;
	int borne_droite = tete_info_dict_kmer->nb_kmer - 1;

//on remplit le vecteur, chaque case contient le pointeur vers le bloc kmer correspondant
	while ( p_kmer != NULL)
	{
		if (i == 0) // cas de l'enregistrement du pointeur vers le premier bloc Kmer, donc on ne doit pas avancer
		{
			p_tableau_kmer_QS[i] = tete_info_dict_kmer->tete_liste_kmer;
			i++;
		}
		else
		{
			p_tableau_kmer_QS[i] = p_kmer->suiv_kmer;
			i++;
			p_kmer = p_kmer->suiv_kmer;
		}
	}

  trier_ST1(adr_tableau_kmer_QS, borne_gauche, borne_droite);
	trier_ST2(adr_tableau_kmer_QS, borne_gauche, borne_droite);
}
//==================================================================================================
//==================================================================================================
void generation_fichier_resultats(TPtr_Cellkmer *tableau_kmer_QS, int numero_essais, TPtr_info_dictionnaire_kmer tete_info_dict_kmer, TPtr_info_dictionnaire_sequences tete_info_dict_seq)
{

	TPtr_Cellkmer p_kmer = tete_info_dict_kmer->tete_liste_kmer;
  int i, j, iter_kmer;
	char nom_fichier[50];
	snprintf( nom_fichier, 50, "fichier_resultatse_%d.txt", numero_essais);

  FILE* fichier_resultats;
  fichier_resultats = fopen(nom_fichier, "w");

	for ( iter_kmer=0; iter_kmer < tete_info_dict_kmer->nb_kmer; iter_kmer++ )
	{
		if (p_kmer->nb_sequence == tete_info_dict_seq->nb_sequences )
		{
			fprintf(fichier_resultats, "\n=============================================== \n");

			fprintf(fichier_resultats, "Kmer: %s      Motif consensus: %s", tableau_kmer_QS[iter_kmer]->kmer, tableau_kmer_QS[iter_kmer]->motif_consensus );
			fprintf(fichier_resultats, " ST1: %d                  ST2: %d", tableau_kmer_QS[iter_kmer]->score_St1, tableau_kmer_QS[iter_kmer]->score_St2);

			fprintf(fichier_resultats, "\n-----------PSSM------------- \n");
			for (i=0;  i<4; i++){
				if(i == 0){ fprintf(fichier_resultats, "A  ");}
				if(i == 1){ fprintf(fichier_resultats, "\nT  ");}
				if(i == 2){ fprintf(fichier_resultats, "\nC  ");}
				if(i == 3){ fprintf(fichier_resultats, "\nG  ");}
				for(j=0; j<taille_motif; j++){
					fprintf(fichier_resultats, "%.2f ", tableau_kmer_QS[iter_kmer]->PSSM_consensus[i][j]);
				}
			}
			fprintf(fichier_resultats, "\n=============================================== \n");
		}
	}
	fclose(fichier_resultats);
}
//==================================================================================================
