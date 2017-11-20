#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <time.h>
#include "fonctions.h"

//------------------------------------------------------------------------------------------------------------
// lecture du fichier contenant les variables souhaitées, stockage de ces variables pour utilisation
void importer_parametres(int* l, int* d, int* k, int* nb_masques)
{
	FILE* ptr_fichier;
  char nom_fichier[30];

	printf("Quel est le nom du fichier que vous voulez utiliser pour importer les paramètres du programme ?\n");
	scanf("%s", nom_fichier);
  ptr_fichier = fopen(nom_fichier, "r");

	if( ptr_fichier == NULL) { free(ptr_fichier); return; }// si le fichier est vide on sort

  fscanf(ptr_fichier, "#longueur du motif à identifier (IE. taille du masque): %d\n", l);
  fscanf(ptr_fichier, "#nombre maximal de substitutions autorisées: %d\n", d);
  fscanf(ptr_fichier, "#nombre de fenêtres dans les masques utilisés: %d\n", k);
  fscanf(ptr_fichier, "#nombre de masques à générer: %d\n", nb_masques);

  fclose(ptr_fichier);
}
//---------------------------------------------------------------------------------------------------------
// Récupère les séquences et les stocke dans une liste chainée
void importer_sequences_fasta( TInfo_ensemble_sequences* ptr_info, TEnsemble_Sequences* ptr_ensemble )
{
  FILE* ptr_fichier_fasta;
  char nom_fichier_fasta[30];
  char contenu_ligne[128];
  char c;
  int cptr = 0;

  TPtr_ensemble_sequences p_new;
  TPtr_ensemble_sequences p = ptr_ensemble;
  ptr_info->nb_seq = 0;

	printf("Quel est le nom du fichier que vous voulez utiliser pour importer les séquences ?\n");
  scanf("%s", nom_fichier_fasta);

  ptr_fichier_fasta = fopen(nom_fichier_fasta, "r");
  if( ptr_fichier_fasta == NULL) { free(ptr_fichier_fasta); return; }// si le fichier est vide on sort

	do //cette boucle permet de compter les séquences
	{
		c = fgetc (ptr_fichier_fasta); //lecture du fichier caractère par caractère
		if (c == '>') cptr++;  //incrémentation pour chaque nom de séquences
	} while (c != EOF);

	ptr_fichier_fasta = fopen(nom_fichier_fasta, "r");
	do
	{
		fscanf(ptr_fichier_fasta,"%s\n", contenu_ligne);
	  if ( contenu_ligne[0] == '>')
	  {
			if ((ptr_info->nb_seq) > 0)
			{
				p_new = malloc ( sizeof(TEnsemble_Sequences));
				p->suiv_seq = p_new;
				p = p_new;
				if ( ptr_info->nb_seq  == cptr) { p->suiv_seq = NULL;}
			}
			ptr_info->nb_seq += 1;
			strcpy(p->nom_seq, contenu_ligne);
			// printf(" DOOM %s\n", p->nom_seq);
		}else{
			if ( p->seq == NULL){
				strcpy(p->seq, contenu_ligne);
	    }else{
	      strcat(p->seq, contenu_ligne);
	    }
	  }
	}	while(!feof(ptr_fichier_fasta));
  fclose(ptr_fichier_fasta);
}
//--------------------------------------------------------------------------------------------------

void afficher_sequences(TInfo_ensemble_sequences* ptr_info, TEnsemble_Sequences* ptr_ensemble )
{
	FILE * ptr_fichier;
  TPtr_ensemble_sequences p = ptr_ensemble;
	TPtr_ensemble_sequences p_prec = NULL;
	int cptr = 0;
	printf("fonction afficher_sequences\n");
	printf("%d \n", ptr_info->nb_seq);

	ptr_fichier = fopen ( "verif_dico_fasta.txt" ,"w");
	while ( cptr < (ptr_info->nb_seq))
	{
		fputs( p->nom_seq, ptr_fichier);
		fputs( "\n" ,ptr_fichier);
		fputs( p->seq,ptr_fichier);
		fputs( "\n" ,ptr_fichier);
		printf("nom_seq: %s\n", p->nom_seq );
		printf("%s\n", p->seq);
		p_prec = p;
		p = p->suiv_seq;
		cptr ++;
	}
}
//------------------------------------------------------------------------------------------------
//genere un chiffre random
int random_number(int max_number, int zero_excluded)
{
	int randomNumber;

	if(zero_excluded==0) //on peut tomber sur 0 aléatoirement
	{
		randomNumber= rand() % max_number;
	}else{
    randomNumber= rand() % max_number +1;
	}
	return(randomNumber);
}
//-------------------------------------------------------------------------------------------------
//genere un masque avec un nombre de fenetre ouverte en paramètre:
void generation_masque(int l, int* masque[l], int k)
{
  int i,j;
  int nb_fenetre_ouverte = 1;
  // while (nb_fenetre_ouverte < k)
	// {
    for (j=0; j<= (l -1); j++)
		{
			// printf("x -");
      masque[j] = random_number(2,0);
      // if (masque[j] == 1)
			// {
      // 	nb_fenetre_ouverte ++;
      // }
    }
		// for(i=0; i< l; i++)
		// {
		// 	printf("%d", masque[i]);
		// }
  // }
}
