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
void generation_masque(int l, int* masque, int k)
{
  int i,j;
  int nb_fenetre_ouverte = 1;
  while (nb_fenetre_ouverte < k)
	{
    for (j=0; j<= (l -1); j++)
		{
      masque[j] = random_number(2,0);
      if (masque[j] == 1)
			{
    		nb_fenetre_ouverte ++;
			}
    }
  }
	//test du masque
	for(i=0; i< l; i++)
	{
		printf("%d\n", masque[i]);
	}
}

//---------------------------------------------------------------------------------------------------------------------
void parcours_masque(int longueur_masque, int* adr_masque[longueur_masque], int nb_fenetre, int nb_sequence,
										ptr_struct_seq* adr_tete_struct_sequence, TPtr_Cellkmer* adr_tete_liste_kmer,
										TPtr_CellSequence* adr_tete_liste_sequence, TPtr_CellPos* adr_tete_liste_pos)
{
  ptr_struct_seq p_generation_seq = *adr_tete_struct_sequence;
  TPtr_Cellkmer p_kmer = *adr_tete_liste_kmer;
  TPtr_CellSequence p_sequence = *adr_tete_liste_sequence;
  TPtr_CellPos p_pos = *adr_tete_liste_pos;
	
  int position = 0;
  int cpt_pos_masque, position_kmer;
  int pos_kmer = 0;
  char k_mer[nb_fenetre + 1]; // le k_mer mesure la taille du nombre de fenêtre ouverte dans le masque

  while (p_generation_seq != NULL)
	{
    while (position < 30)
		{
      for(cpt_pos_masque = 0; cpt_pos_masque< longueur_masque; cpt_pos_masque++)  // on parcourt les nucleotides sous le masque
			{
        if (adr_masque[cpt_pos_masque] == 1) //si la fenêtre du masque est ouverte
				{
          k_mer[pos_kmer] = p_generation_seq->sequence[position];

          if (pos_kme == 0)
					{
            position_kmer = position;
          }
          pos_kmer ++;
          if (pos_kmer == nb_fenetre)
					{
            k_mer[pos_kmer] = '\0' ;
            pos_kmer = 0;
          }

        }
        position = position +1;
      }
    	generation_kmer(position_kmer, k_mer, &p_kmer, &p_sequence, &p_generation_seq, &p_pos);
    }
    position = 0;
    p_generation_seq = p_generation_seq->next_sequence;
  }
  return;

}


//----------------------------------------------------------------------------------------------------------------
void generation_kmer(int position_kmer, char* k_mer, TPtr_Cellkmer* adr_liste_kmer,
						TPtr_CellSequence* adr_liste_sequence, ptr_struct_seq* adr_liste_generation_sequence,
						TPtr_CellPos* adr_liste_pos)
{
	TPtr_Cellkmer temp_p_kmer = *adr_liste_kmer; //on créé des pointeurs temporaire pour parcourir les listes
	ptr_struct_seq p_generation_seq = *adr_liste_generation_sequence;

	if (temp_p_kmer->suiv_kmer == NULL) //cas du premier element de la liste
	{
	  TPtr_Cellkmer nouveau_kmer = (TCellkmer*) malloc(sizeof(TCellkmer));
	  strcpy(temp_p_kmer->kmer, k_mer);
	  nouveau_kmer->suiv_kmer = NULL;

	  TPtr_CellSequence nouvelle_sequence = (TCellSequence*) malloc(sizeof(TCellSequence));
	  temp_p_kmer->suiv_kmer = nouveau_kmer;
	  temp_p_kmer->tete_sequence = nouvelle_sequence;
	  nouvelle_sequence->sequence = p_generation_seq->numero_sequence;
	  nouvelle_sequence->suiv_sequence = NULL;

	  TPtr_CellPos nouvelle_position = (TCellPos*) malloc(sizeof(TCellPos));
	  nouvelle_sequence->tete_pos = nouvelle_position;
	  nouvelle_position->position = position_kmer;
	  nouvelle_position->suiv_pos = NULL;
	  return;

	}else{

	  while (temp_p_kmer->suiv_kmer != NULL)
		{
      if (strcmp(temp_p_kmer->kmer, k_mer) == 0)  //si le kmer a déjà été trouvé
			{
        TPtr_CellSequence p_liste_sequence = temp_p_kmer->tete_sequence;

        //premier element de la liste chainee de sequence:
        while (p_liste_sequence->suiv_sequence != NULL)
				{
            if (p_generation_seq->numero_sequence == p_liste_sequence->sequence) //Si le kmer a déjà été trouvé dans cette séquence
						{
              TPtr_CellPos p_liste_pos = p_liste_sequence->tete_pos;
              while (p_liste_pos->suiv_pos != NULL)
							{
              	p_liste_pos = p_liste_pos->suiv_pos;
              }
              TPtr_CellPos nouvelle_position = (TCellPos*) malloc(sizeof(TCellPos)); // On créé une nouvelle brique de pôsition
              nouvelle_position->position = position_kmer;
              nouvelle_position->suiv_pos = NULL;
              p_liste_pos->suiv_pos = nouvelle_position;
              return;

            }
            p_liste_sequence = p_liste_sequence->suiv_sequence;
        }
        if (p_generation_seq->numero_sequence == p_liste_sequence->sequence) //Si le kmer a déjà été trouvé dans cette séquence
				{
          TPtr_CellPos p_liste_pos = p_liste_sequence->tete_pos;
          while (p_liste_pos->suiv_pos != NULL)
					{
            p_liste_pos=p_liste_pos->suiv_pos;
          }
          TPtr_CellPos nouvelle_position = (TCellPos*) malloc(sizeof(TCellPos)); // On créé une nouvelle brique de pôsition
          nouvelle_position->position = position_kmer;
          nouvelle_position->suiv_pos = NULL;
          p_liste_pos->suiv_pos = nouvelle_position;
          return;

        }
        TPtr_CellSequence nouvelle_sequence = (TCellSequence*) malloc(sizeof(TCellSequence));
        nouvelle_sequence->sequence = p_generation_seq->numero_sequence;
        nouvelle_sequence->suiv_sequence = NULL;
        p_liste_sequence->suiv_sequence = nouvelle_sequence;

        TPtr_CellPos nouvelle_tete_pos = (TCellPos*) malloc(sizeof(TCellPos)); // on créé une nouvelle tete de liste de position
        nouvelle_sequence->tete_pos = nouvelle_tete_pos;
        nouvelle_tete_pos->position = position_kmer;
        nouvelle_tete_pos->suiv_pos = NULL;
        return;

      }
      temp_p_kmer=temp_p_kmer->suiv_kmer;
	  }

	 // Si ce kmer n'a pas encore été trouvé: ajout en fin de boucle
	  TPtr_Cellkmer nouveau_kmer = (TCellkmer*) malloc(sizeof(TCellkmer));
	  strcpy(temp_p_kmer->kmer, k_mer);
	  nouveau_kmer->suiv_kmer = NULL;

	  TPtr_CellSequence nouvelle_sequence =(TCellSequence*) malloc(sizeof(TCellSequence));
	  temp_p_kmer->suiv_kmer = nouveau_kmer;
	  temp_p_kmer->tete_sequence = nouvelle_sequence;
	  nouvelle_sequence->sequence = p_generation_seq->numero_sequence;
	  nouvelle_sequence->suiv_sequence = NULL;

	  TPtr_CellPos nouvelle_position = (TCellPos*) malloc(sizeof(TCellPos));
	  nouvelle_sequence->tete_pos = nouvelle_position;
	  nouvelle_position->position = position_kmer;
	  nouvelle_position->suiv_pos = NULL;
	  return;

	}
}
