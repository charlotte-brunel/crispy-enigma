#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "fonctions_algo_test.h"

//Fonction insert motif à enlever
void insert_motif(ptr_struct_seq* adr_tete_structure, int nb_sequence, ptr_liste_motif* tete_liste_motif)
{
  char subst = ' ';
  int i,k; // compteur
  ptr_struct_seq p = *adr_tete_structure;
  ptr_liste_motif p_motif = *tete_liste_motif;
  p_motif->next_motif = NULL;
  int position_subst = 0;
  int pos = 0; //position correspondant à la position où va être inséré le motif dans la séquence (entre 0 et 23)
  int nb_subst = 0; //nb_subst correspondant au nombre de substitutions dans le motif

  for (i=1; i <= nb_sequence; i++)
  {
    char motif[7] = "aaaaaa\0";
    pos = random_number(23, 0);
    nb_subst = random_number(3, 0);
    switch(nb_subst)
    {
      case 1:
        position_subst = random_number(5,0);
        random_nucleotide(&subst);
        motif[position_subst] = subst;
        break;
      case 2:
        position_subst = random_number(5,0);
        random_nucleotide(&subst);
        motif[position_subst] = subst;
        position_subst = random_number(5,0);
        random_nucleotide(&subst);
        motif[position_subst] = subst;
        break;
    }
    for (k=0; k < 6; k++)
    {
      p->sequence[pos+k] = motif[k]; // On remplace la séquence actuelle par une occurrence du motif
    }
    p = p->next_sequence;
    if (i != nb_sequence)
    {
      ptr_liste_motif p_motif_next = (liste_chaine_motif*)malloc(sizeof(liste_chaine_motif)); // On créer le prochain élément de la liste chainée de motifs
      strcpy(p_motif->motif_substitue, motif);
      p_motif->next_motif = p_motif_next;
      p_motif = p_motif_next;
    }else{
      strcpy(p_motif->motif_substitue, motif);
      p_motif->next_motif = NULL;
    }
  }
  return;
}

//---------------------------------------------------------------------------------------------------------------------
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
//---------------------------------------------------------------------------------------------------------------------
void random_nucleotide(void *adr_nucleo)
{
	char *p_nucleo = adr_nucleo;
	char nucleotide[4] = "atcg"; //on créer un array contenant les nucleotide
	int randomIndex = rand() % 4;
	*p_nucleo = nucleotide[randomIndex];
	return;
}
//---------------------------------------------------------------------------------------------------------------------
void generation_masque(int longueur_masque, void* adr_masque, int nb_fenetre)
{
  int i;
  int *p_masque = adr_masque;
  int nb_fenetre_ouverte=0;
  while (nb_fenetre_ouverte != nb_fenetre)
  {
  	nb_fenetre_ouverte = 0;
    for (i=0; i <= (longueur_masque-1); i++)
    {
      p_masque[i] = random_number(2,0);
      if (p_masque[i] == 1)
      {
        nb_fenetre_ouverte++;
      }
    }
  }
  return;
}
//---------------------------------------------------------------------------------------------------------------------
void generation_kmer(int position_kmer, char* k_mer, TPtr_Cellkmer* adr_liste_kmer, TPtr_CellSequence* adr_liste_sequence,
   ptr_struct_seq* adr_liste_generation_sequence, TPtr_CellPos* adr_liste_pos)
{
  TPtr_Cellkmer temp_p_kmer = *adr_liste_kmer; //on créer des pointeurs temporaire pour parcourir les listes
  ptr_struct_seq p_generation_seq = *adr_liste_generation_sequence;

  if (temp_p_kmer->suiv_kmer == NULL) //cas du premier element de la liste
  {
    TPtr_Cellkmer nouveau_kmer = (TCellkmer*)malloc(sizeof(TCellkmer));
    strcpy(temp_p_kmer->kmer, k_mer);
    nouveau_kmer->suiv_kmer = NULL;
    TPtr_CellSequence nouvelle_sequence=(TCellSequence*)malloc(sizeof(TCellSequence));
    temp_p_kmer->suiv_kmer = nouveau_kmer;
    temp_p_kmer->tete_sequence = nouvelle_sequence;
    nouvelle_sequence->sequence = p_generation_seq->numero_sequence;
    nouvelle_sequence->suiv_sequence = NULL;
    TPtr_CellPos nouvelle_position = (TCellPos*)malloc(sizeof(TCellPos));
    nouvelle_sequence->tete_pos = nouvelle_position;
    nouvelle_position->position = position_kmer;
    nouvelle_position->suiv_pos = NULL;
    return;
  }
  else
  {
    while (temp_p_kmer->suiv_kmer != NULL)
    {
      if (strcmp(temp_p_kmer->kmer, k_mer) == 0) //si le kmer a déjà été trouvé
      {
        TPtr_CellSequence p_liste_sequence = temp_p_kmer->tete_sequence;
        //premier element de la liste chainée de séquences:
        while (p_liste_sequence->suiv_sequence != NULL)
        {
          if (p_generation_seq->numero_sequence == p_liste_sequence->sequence) //Si le kmer a déjà été trouvé dans cette séquence
          {
            TPtr_CellPos p_liste_pos = p_liste_sequence->tete_pos;
            while (p_liste_pos->suiv_pos != NULL)
            {
              p_liste_pos=p_liste_pos->suiv_pos;
            }
            TPtr_CellPos nouvelle_position = (TCellPos*)malloc(sizeof(TCellPos)); // On créer une nouvelle brique de position
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
          TPtr_CellPos nouvelle_position = (TCellPos*)malloc(sizeof(TCellPos)); // On créer une nouvelle brique de position
          nouvelle_position->position = position_kmer;
          nouvelle_position->suiv_pos = NULL;
          p_liste_pos->suiv_pos = nouvelle_position;
          return;
        }
        TPtr_CellSequence nouvelle_sequence = (TCellSequence*)malloc(sizeof(TCellSequence));
        nouvelle_sequence->sequence = p_generation_seq->numero_sequence;
        nouvelle_sequence->suiv_sequence = NULL;
        p_liste_sequence->suiv_sequence = nouvelle_sequence;
        TPtr_CellPos nouvelle_tete_pos = (TCellPos*)malloc(sizeof(TCellPos)); // on créer une nouvelle tete de liste de position
        nouvelle_sequence->tete_pos = nouvelle_tete_pos;
        nouvelle_tete_pos->position = position_kmer;
        nouvelle_tete_pos->suiv_pos = NULL;
        return;
      }
      temp_p_kmer = temp_p_kmer->suiv_kmer;
    }

   // Si ce kmer n'a pas encore été trouvé: ajout en fin de boucle
    TPtr_Cellkmer nouveau_kmer = (TCellkmer*)malloc(sizeof(TCellkmer));
    strcpy(temp_p_kmer->kmer, k_mer);
    nouveau_kmer->suiv_kmer = NULL;
    TPtr_CellSequence nouvelle_sequence = (TCellSequence*)malloc(sizeof(TCellSequence));
    temp_p_kmer->suiv_kmer = nouveau_kmer;
    temp_p_kmer->tete_sequence = nouvelle_sequence;
    nouvelle_sequence->sequence = p_generation_seq->numero_sequence;
    nouvelle_sequence->suiv_sequence = NULL;
    TPtr_CellPos nouvelle_position = (TCellPos*)malloc(sizeof(TCellPos));
    nouvelle_sequence->tete_pos = nouvelle_position;
    nouvelle_position->position = position_kmer;
    nouvelle_position->suiv_pos = NULL;
    return;
  }
}


//---------------------------------------------------------------------------------------------------------------------
void parcours_masque(int longueur_masque, void* adr_masque, int nb_fenetre, int nb_sequence, ptr_struct_seq* adr_tete_struct_sequence,
   TPtr_Cellkmer* adr_tete_liste_kmer, TPtr_CellSequence* adr_tete_liste_sequence, TPtr_CellPos* adr_tete_liste_pos)
{
  ptr_struct_seq p_generation_seq = *adr_tete_struct_sequence;
  TPtr_Cellkmer p_kmer = *adr_tete_liste_kmer;
  TPtr_CellSequence p_sequence = *adr_tete_liste_sequence;
  TPtr_CellPos p_pos = *adr_tete_liste_pos;
  int *p_masque = adr_masque;
  int position = 0;
  int cpt_pos_masque, position_kmer;
  int pos_kmer = 0;
  char k_mer[nb_fenetre+1]; // le k_mer mesure la taille du nombre de fenêtres ouvertes dans le masque

  while (p_generation_seq != NULL)
  {
    while (position < 30)
    {
      for(cpt_pos_masque = 0; cpt_pos_masque < longueur_masque; cpt_pos_masque++) // on parcourt les nucleotides sous le masque
      {
        if (p_masque[cpt_pos_masque] == 1){ //si la fenêtre du masque est ouverte
          k_mer[pos_kmer] = p_generation_seq->sequence[position];
          if (pos_kmer == 0)
          {
            k_mer[pos_kmer] = '\0' ;
            pos_kmer = 0;
            position_kmer = position;
          }
          pos_kmer++;
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

//---------------------------------------------------------------------------------------------------------------------
//Procédure permettant de remplir le dictionnaire de kmers sélectionnés pour le calcul de la PSSM:
void recuperer_motif_kmer(TPtr_Cellkmer* adr_parcours_kmer, TPtr_Cellkmer_selectionne *adr_tete_kmer_selectionne,
  TPtr_Cell_Motif_PSSM* adr_tete_motif_PSSM, TPtr_CellSequence* adr_cell_sequence, TPtr_CellPos* adr_cell_pos,
   ptr_struct_seq* adr_generation_sequence, int nb_sequence_kmer)
{
  TPtr_Cellkmer p_kmer = *adr_parcours_kmer;
  printf("pkmer: %s \n", p_kmer->kmer);
  TPtr_Cellkmer_selectionne p_kmer_selectionne = *adr_tete_kmer_selectionne;
  ptr_struct_seq tete_generation_seq = *adr_generation_sequence;
  ptr_struct_seq p_generation_seq = tete_generation_seq;
  int longueur_sequence = 30;
  int longueur_motif = 5;
  int i, sequence_actuelle;

  if (p_kmer_selectionne->suiv_kmer_selectionne == NULL) //la liste est vide:
  {
    printf("premier tour de boucle \n");
    TPtr_Cellkmer_selectionne prochain_kmer = malloc(sizeof(TCellkmer_selectionne)); //création de la prochaine brique
    TPtr_Cell_Motif_PSSM nouveau_motif = malloc(sizeof(TCell_Motif_PSSM));
    strcpy(p_kmer_selectionne->kmer, p_kmer->kmer);//remplissage la brique actuelle
    printf("KMER: %s \n", p_kmer_selectionne->kmer);
    p_kmer_selectionne->nb_sequence = nb_sequence_kmer;
    printf("nb sequence: %d \n", p_kmer_selectionne->nb_sequence);
    p_kmer_selectionne->suiv_kmer_selectionne = prochain_kmer; // chainage entre la brique actuelle et la prochaine
    p_kmer_selectionne->tete_motif_PSSM = nouveau_motif; //récupération de la bonne tete
    TPtr_CellSequence p_sequence = p_kmer->tete_sequence;
    do
    {
      printf("parcours des sequences \n");
      TPtr_CellPos p_pos = p_sequence->tete_pos; // récupération de la liste de position pour chaque sequence
      sequence_actuelle = p_sequence->sequence;
      printf("p_generation_seq: %d , sequence actuelle: %d \n", p_generation_seq->numero_sequence, sequence_actuelle);
      while (p_generation_seq->numero_sequence != sequence_actuelle)  //on se positionne dans la bonne sequence
      {
        printf("p_generation_seq: %d \n", p_generation_seq->numero_sequence);
        p_generation_seq=p_generation_seq->next_sequence;
      }
      printf("numero sequence : %d, sequence: %s \n", p_generation_seq->numero_sequence, p_generation_seq->sequence);
      if((p_pos->position + longueur_motif)< longueur_sequence)
      {
        for (i=(p_pos->position); i<(p_pos->position + longueur_motif); i++)
        {
          nouveau_motif->motif[i-(p_pos->position)]=p_generation_seq->sequence[i];
        }
        nouveau_motif->motif[longueur_motif]= '\0';

        printf("Motif: %s \n", nouveau_motif->motif);
        TPtr_Cell_Motif_PSSM prochain_motif = (TCell_Motif_PSSM*)malloc(sizeof(TCell_Motif_PSSM));
        nouveau_motif->suiv_motif = prochain_motif;
        nouveau_motif = prochain_motif;
      }
      p_sequence = p_sequence->suiv_sequence;
    } while (p_sequence != NULL);
    p_generation_seq = tete_generation_seq;
    return;
  }
  printf("I have several kmer !!");
  do //Cas où on a plusieurs kmers selectionnés
  {
    p_kmer_selectionne = p_kmer_selectionne->suiv_kmer_selectionne;
    printf("KMER: %s \n", p_kmer_selectionne->kmer);
  } while (p_kmer_selectionne->suiv_kmer_selectionne != NULL);

  TPtr_Cellkmer_selectionne prochain_kmer = (TCellkmer_selectionne*)malloc(sizeof(TCellkmer_selectionne));
  TPtr_Cell_Motif_PSSM nouveau_motif = (TCell_Motif_PSSM*)malloc(sizeof(TCell_Motif_PSSM));
  strcpy(p_kmer_selectionne->kmer, p_kmer->kmer);
  p_kmer_selectionne->nb_sequence = nb_sequence_kmer;
  p_kmer_selectionne->suiv_kmer_selectionne = prochain_kmer;
  p_kmer_selectionne->tete_motif_PSSM = nouveau_motif;
  TPtr_CellSequence p_sequence = p_kmer->tete_sequence;
  do
  {
    TPtr_CellPos p_pos = p_sequence->tete_pos; // On récupère la liste de position pour chaque séquence
    sequence_actuelle = p_sequence->sequence;
    while (p_generation_seq->numero_sequence != sequence_actuelle)
    {
      p_generation_seq = p_generation_seq->next_sequence;
    } //On se place dans la bonne séquence
    printf("position motif: %d \n", p_pos->position);
    if((p_pos->position + longueur_motif) < longueur_sequence)
    {
      for (i = (p_pos->position); i < (p_pos->position + longueur_motif); i++)
      {
        nouveau_motif->motif[i-(p_pos->position)]=p_generation_seq->sequence[i];
      }
      nouveau_motif->motif[longueur_motif] = '\0';

      printf("Motif: %s \n", nouveau_motif->motif);

      TPtr_Cell_Motif_PSSM prochain_motif = malloc(sizeof(TCell_Motif_PSSM));
      nouveau_motif->suiv_motif = prochain_motif;
      nouveau_motif = prochain_motif;
    }
    p_sequence = p_sequence->suiv_sequence;
  } while (p_sequence != NULL);
  p_generation_seq = tete_generation_seq;
  return;
}

//---------------------------------------------------------------------------------------------------------------------
//Fonction vérifiant si le kmer donné est trouvé dans chaque séquence:
void kmer_present_dans_chaque_sequence(int nb_sequence, TPtr_Cellkmer* adr_cell_kmer, TPtr_CellSequence *adr_cell_sequence,
  TPtr_CellPos *adr_cell_pos, ptr_struct_seq* adr_cell_generation_sequence, TPtr_Cellkmer_selectionne* adr_cell_kmer_selectionne,
  TPtr_Cell_Motif_PSSM* adr_cell_motif_PSSM)
{
  TPtr_Cellkmer p_parcours_kmer = *adr_cell_kmer;
  TPtr_Cellkmer_selectionne tete_kmer_selectionne = *adr_cell_kmer_selectionne;
  TPtr_Cellkmer_selectionne p_kmer_selectionne = *adr_cell_kmer_selectionne;
  p_kmer_selectionne->suiv_kmer_selectionne = NULL;
  TPtr_Cell_Motif_PSSM p_motif_PSSM = *adr_cell_motif_PSSM;
  TPtr_CellSequence tete_sequence = *adr_cell_sequence;
  TPtr_CellPos tete_pos = *adr_cell_pos;
  ptr_struct_seq tete_generation_sequence = *adr_cell_generation_sequence;
  int nb_sequence_par_kmer = 0;

  while (p_parcours_kmer != NULL)
  {
    TPtr_CellSequence p_parcours_sequence = p_parcours_kmer->tete_sequence;
    while (p_parcours_sequence != NULL )
    {
      nb_sequence_par_kmer++;
      p_parcours_sequence = p_parcours_sequence->suiv_sequence;
    }

    printf("nb_sequence_par_kmer: %d \n", nb_sequence_par_kmer);
    p_parcours_kmer->nb_sequence = nb_sequence_par_kmer;
    if (nb_sequence_par_kmer >= 7)  //Pour l'instant pour pouvoir continuer
    {
      recuperer_motif_kmer(&p_parcours_kmer, &p_kmer_selectionne, &p_motif_PSSM, &tete_sequence, &tete_pos, &tete_generation_sequence, nb_sequence_par_kmer);
      p_kmer_selectionne = tete_kmer_selectionne;
    }
    nb_sequence_par_kmer = 0;
    p_parcours_kmer = p_parcours_kmer->suiv_kmer;
  }
  return;
}

//---------------------------------------------------------------------------------------------------------------------
void affichage_dictionnaire_kmer(TPtr_Cellkmer* adr_tete_kmer, TPtr_CellSequence* adr_tete_sequence, TPtr_CellPos* adr_tete_pos)
{
  TPtr_Cellkmer p_kmer = *adr_tete_kmer;
  FILE* fichier_dictionnaire = NULL;
  fichier_dictionnaire = fopen("dictionnaire_kmer.txt", "w");

  while (p_kmer != NULL)
  {
    TPtr_CellSequence p_sequence = p_kmer->tete_sequence;
    fprintf(fichier_dictionnaire, "KMER: %s \n", p_kmer->kmer);
    while (p_sequence != NULL)
    {
      TPtr_CellPos p_pos = p_sequence->tete_pos;
      fprintf(fichier_dictionnaire, "Sequence: %d \n", p_sequence->sequence);
      while (p_pos != NULL)
      {
        fprintf(fichier_dictionnaire, "Position: %d \n ", p_pos->position);
        p_pos = p_pos->suiv_pos;
      }
      p_sequence = p_sequence->suiv_sequence;
    }
    fprintf(fichier_dictionnaire, "\n\n\n");
    p_kmer = p_kmer->suiv_kmer;
  }
  fclose(fichier_dictionnaire);
  return;
}

//---------------------------------------------------------------------------------------------------------------------
void affichage_motif_selectionne(TPtr_Cellkmer_selectionne* adr_tete_kmer_selectionne,
  TPtr_Cell_Motif_PSSM* adr_tete_motif)
{
  printf("I m here");
  TPtr_Cellkmer_selectionne p_kmer_selectionne = *adr_tete_kmer_selectionne;
  FILE* fichier_kmer_selectionne = NULL;
  printf("KMER: %s \n", p_kmer_selectionne->kmer);
  fichier_kmer_selectionne = fopen("kmer_selectionne.txt", "w");

  while (p_kmer_selectionne != NULL)
  {
    TPtr_Cell_Motif_PSSM p_motif = p_kmer_selectionne->tete_motif_PSSM;
    fprintf(fichier_kmer_selectionne, "KMER: %s \n", p_kmer_selectionne->kmer);
    fprintf(fichier_kmer_selectionne, "Present dans %d sequences \n", p_kmer_selectionne->nb_sequence);
    while(p_motif != NULL)
    {
      fprintf(fichier_kmer_selectionne, "MOTIF: %s \n", p_motif->motif);
      p_motif = p_motif->suiv_motif;
    }
    p_kmer_selectionne = p_kmer_selectionne->suiv_kmer_selectionne;
  }
  fclose(fichier_kmer_selectionne);
  return;
}

//---------------------------------------------------------------------------------------------------------------------
void calcul_PSSM(TPtr_Cellkmer_selectionne *adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM *adr_cell_motif_PSSM, FILE** file_info, double*** adr_matrice_PSSM)
{
  printf("now i m here \n");
  *file_info = fopen("PSSM_Motif_Trouve.txt", "w"); // Dans ce fichier on va écrire la PSSM prévisionnelle
  TPtr_Cellkmer_selectionne p_kmer_selectionne = *adr_cell_kmer_selectionne;
  double** p_matrice_pssm = *adr_matrice_PSSM;
  double nb_sequence = p_kmer_selectionne->nb_sequence;
  int taille_motif = 5;
  int nb_ligne;
  printf("LET DO IIIIT");
/* ****************
 * MATRICE PSSM *
 *****************/

  p_matrice_pssm = malloc(4 * sizeof(*p_matrice_pssm));
  if (p_matrice_pssm == NULL){
    free(p_matrice_pssm);
  }

  for(nb_ligne=0 ; nb_ligne < 4 ; nb_ligne++)
  {
    p_matrice_pssm[nb_ligne] = malloc(taille_motif * sizeof(*(p_matrice_pssm[nb_ligne]))); //On alloue des tableaux de 'taille2' variables.
    if(p_matrice_pssm[nb_ligne] == NULL)
    {
      for(nb_ligne=0 ; nb_ligne < 4 ; nb_ligne++){
        free(p_matrice_pssm[nb_ligne]);
      }
    }
  }

  double add = 1/nb_sequence;
  char to_print;
  double maximum = -1;
  int i, j, k, a, t, c, g;
  //initialisation de la matrice à 0:
  for (j=0; j < 4; j++){
    for(k=0; k < taille_motif; k++){
      p_matrice_pssm[j][k]= 0;
    }
  }
  //calcul de la PSSM à partir de la liste chainée de motif
  TPtr_Cell_Motif_PSSM p_motif = p_kmer_selectionne->tete_motif_PSSM;
  do // on remplit la PSSM
  {
    for (i=0; i < 6; i++)
    {
      switch(p_motif->motif[i])
      {
        case 'a':
          p_matrice_pssm[0][i] = p_matrice_pssm[0][i] + add;
          break;
        case 't':
          p_matrice_pssm[1][i] = p_matrice_pssm[1][i] + add;
          break;
        case 'c':
          p_matrice_pssm[2][i] = p_matrice_pssm[2][i] + add;
          break;
        case 'g':
          p_matrice_pssm[3][i] = p_matrice_pssm[3][i] + add;
          break;
      }
    }
    p_motif = p_motif->suiv_motif;
  } while (p_motif != NULL);

  //écriture de la matrice PSSM dans le fichier d'info.
  fprintf(*file_info, "\n\nPSSM: \n");
  fprintf(*file_info, "a  ");
  for (a=0; a < taille_motif; a++){
    fprintf(*file_info, "%.2f ", p_matrice_pssm[0][a]);
  }
  fprintf(*file_info, "\nt  ");
  for (t=0; t < taille_motif; t++){
    fprintf(*file_info, "%.2f ", p_matrice_pssm[1][t]);
  }
  fprintf(*file_info, "\nc  ");
  for (c=0; c < taille_motif; c++){
    fprintf(*file_info, "%.2f ", p_matrice_pssm[2][c]);
  }
  fprintf(*file_info, "\ng  ");
  for (g=0; g < taille_motif; g++){
    fprintf(*file_info, "%.2f ", p_matrice_pssm[3][g]);
  }

  //écriture du motif consensus dans le fichier à partir de la matrice PSSM:
  fprintf(*file_info, "\n\nMotif Consensus: \n" );
  for (j=0; j < taille_motif; j++)
  {
    for(k=0; k < 4; k++)
    {
      if (p_matrice_pssm[k][j] > maximum)
      {
        maximum = p_matrice_pssm[k][j];
        switch (k)
        {
          case 0:
            to_print = 'a';
            break;
          case 1:
            to_print = 't';
            break;
          case 2:
            to_print = 'c';
            break;
          case 3:
            to_print = 'g';
            break;
        }
      }
    }
    maximum = -1; //on remet maximum à -1 avant d'évaluer la nucléotide majoritaire de la séquence suivante !
    fprintf(*file_info, "%c", to_print);
  }
  fprintf(*file_info, "\n\n\n");

  fclose(*file_info);
  return;
}