#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "fonctions_j.h"
//------------------------------------------------------------------------------------------------------------
//Fonction insert motif � enlever
void insert_motif(ptr_struct_seq* adr_tete_structure, int nb_sequence, ptr_liste_motif* tete_liste_motif)
{
  char subst= ' ';
  int i,k; // compteur
  ptr_struct_seq p= *adr_tete_structure;
  ptr_liste_motif p_motif= *tete_liste_motif;
  p_motif->next_motif= NULL;
  int position_subst=0;
  int pos=0; //position correspond � la position o� va �tre ins�r� le motif dans la s�quence (entre 0 et 23)
  int nb_subst=0; //nb_subst correspond au nombre de substitution dans le motif
  for (i=1; i<=nb_sequence; i++)
  {
    char motif[7]="aaaaaa\0";
    pos=random_number(23, 0);
    nb_subst=random_number(3, 0);
    switch(nb_subst)
    {
      case 1:
        position_subst= random_number(5,0);
        random_nucleotide(&subst);
        motif[position_subst]=subst;
        break;
      case 2:
        position_subst= random_number(5,0);
        random_nucleotide(&subst);
        motif[position_subst]=subst;
        position_subst= random_number(5,0);
        random_nucleotide(&subst);
        motif[position_subst]=subst;
        break;
    }
    for (k=0; k<(6); k++)
    {
      p->sequence[pos+k]= motif[k]; // On remplace la s�quence actuelle par une occurrence du motif
    }
    p=p->next_sequence;
    if (i != nb_sequence)
    {
      ptr_liste_motif p_motif_next= (liste_chaine_motif*)malloc(sizeof(liste_chaine_motif)); // On cr�e le prochain �l�ment de la liste chain�e de motif
      strcpy(p_motif->motif_substitue, motif);
      p_motif-> next_motif= p_motif_next;
      p_motif=p_motif_next;
    }else{
      strcpy(p_motif->motif_substitue, motif);
      p_motif-> next_motif= NULL;
    }
  }
  return;
}
//------------------------------------------------------------------------------------------------------------
 int random_number(int max_number, int zero_excluded)
 {
	int randomNumber;
	if(zero_excluded==0)
  { //on peut tomber sur 0 al�atoirement
		randomNumber= rand() % max_number;
	}else{
    randomNumber= rand() % max_number +1;
	}
	return(randomNumber);
}
//------------------------------------------------------------------------------------------------------------
void random_nucleotide(void *adr_nucleo)
{
	char *p_nucleo= adr_nucleo;
	char nucleotide[4] = "atcg"; //on cr�e un array contenant les nucleotide
	int randomIndex= rand() % 4;
	*p_nucleo= nucleotide[randomIndex];
	return;
	}
//------------------------------------------------------------------------------------------------------------
void generation_masque(int longueur_masque, void* adr_masque, int nb_fenetre)
{
  int i;
  int *p_masque= adr_masque;
  int nb_fenetre_ouverte=0;
  while (nb_fenetre_ouverte != nb_fenetre)
  {
  	nb_fenetre_ouverte=0;
    for (i=0; i<= (longueur_masque -1); i++)
    {
      p_masque[i]= random_number(2,0);
      if (p_masque[i]==1)
      {
        nb_fenetre_ouverte++;
      }
    }
  }
  return;
}
//------------------------------------------------------------------------------------------------------------
void generation_kmer(int position_kmer, char* k_mer, TPtr_Cellkmer* adr_liste_kmer, TPtr_CellSequence* adr_liste_sequence, ptr_struct_seq* adr_liste_generation_sequence, TPtr_CellPos* adr_liste_pos)
{
  TPtr_Cellkmer temp_p_kmer= *adr_liste_kmer; //on cr�� des pointeurs temporaire pour parcourir les listes
  ptr_struct_seq p_generation_seq= *adr_liste_generation_sequence;
  if (temp_p_kmer->suiv_kmer==NULL) //cas du premier element de la liste
  {
    TPtr_Cellkmer nouveau_kmer=(TCellkmer*)malloc(sizeof(TCellkmer));
    strcpy(temp_p_kmer->kmer, k_mer);
    nouveau_kmer->suiv_kmer= NULL;
    TPtr_CellSequence nouvelle_sequence=(TCellSequence*)malloc(sizeof(TCellSequence));
    temp_p_kmer->suiv_kmer= nouveau_kmer;
    temp_p_kmer->tete_sequence=nouvelle_sequence;
    nouvelle_sequence->sequence=p_generation_seq->numero_sequence;
    nouvelle_sequence->suiv_sequence= NULL;
    TPtr_CellPos nouvelle_position=(TCellPos*)malloc(sizeof(TCellPos));
    nouvelle_sequence->tete_pos= nouvelle_position;
    nouvelle_position->position=position_kmer;
    nouvelle_position->suiv_pos= NULL;
    return;
  }
  else
  {
    while (temp_p_kmer->suiv_kmer != NULL)
    {
      if (strcmp(temp_p_kmer->kmer, k_mer)==0)//si le kmer a d�j� �t� trouv�
      {
        TPtr_CellSequence p_liste_sequence=temp_p_kmer->tete_sequence;
        //premier element de la liste chainee de sequence:
        while (p_liste_sequence->suiv_sequence != NULL)
        {
          if (p_generation_seq->numero_sequence==p_liste_sequence->sequence)
          { //Si le kmer a d�j� �t� trouv� dans cette s�quence
            TPtr_CellPos p_liste_pos= p_liste_sequence->tete_pos;
            while (p_liste_pos->suiv_pos != NULL)
            {
              p_liste_pos=p_liste_pos->suiv_pos;
            }
            TPtr_CellPos nouvelle_position= (TCellPos*)malloc(sizeof(TCellPos)); // On cr�� une nouvelle brique de p�sition
            nouvelle_position->position=position_kmer;
            nouvelle_position->suiv_pos=NULL;
            p_liste_pos->suiv_pos= nouvelle_position;
            return;
          }
          p_liste_sequence=p_liste_sequence->suiv_sequence;
        }
        if (p_generation_seq->numero_sequence==p_liste_sequence->sequence)
        { //Si le kmer a d�j� �t� trouv� dans cette s�quence
          TPtr_CellPos p_liste_pos= p_liste_sequence->tete_pos;
          while (p_liste_pos->suiv_pos != NULL)
          {
            p_liste_pos=p_liste_pos->suiv_pos;
          }
          TPtr_CellPos nouvelle_position= (TCellPos*)malloc(sizeof(TCellPos)); // On cr�� une nouvelle brique de p�sition
          nouvelle_position->position=position_kmer;
          nouvelle_position->suiv_pos=NULL;
          p_liste_pos->suiv_pos= nouvelle_position;
          return;
        }
        TPtr_CellSequence nouvelle_sequence= (TCellSequence*)malloc(sizeof(TCellSequence));
        nouvelle_sequence->sequence= p_generation_seq->numero_sequence;
        nouvelle_sequence->suiv_sequence=NULL;
        p_liste_sequence->suiv_sequence= nouvelle_sequence;
        TPtr_CellPos nouvelle_tete_pos= (TCellPos*)malloc(sizeof(TCellPos)); // on cr�� une nouvelle tete de liste de position
        nouvelle_sequence->tete_pos= nouvelle_tete_pos;
        nouvelle_tete_pos->position=position_kmer;
        nouvelle_tete_pos->suiv_pos= NULL;
        return;
      }
      temp_p_kmer=temp_p_kmer->suiv_kmer;
    }
   // Si ce kmer n'a pas encore �t� trouv�: ajout en fin de boucle
    TPtr_Cellkmer nouveau_kmer=(TCellkmer*)malloc(sizeof(TCellkmer));
    strcpy(temp_p_kmer->kmer, k_mer);
    nouveau_kmer->suiv_kmer= NULL;
    TPtr_CellSequence nouvelle_sequence=(TCellSequence*)malloc(sizeof(TCellSequence));
    temp_p_kmer->suiv_kmer= nouveau_kmer;
    temp_p_kmer->tete_sequence=nouvelle_sequence;
    nouvelle_sequence->sequence=p_generation_seq->numero_sequence;
    nouvelle_sequence->suiv_sequence= NULL;
    TPtr_CellPos nouvelle_position=(TCellPos*)malloc(sizeof(TCellPos));
    nouvelle_sequence->tete_pos= nouvelle_position;
    nouvelle_position->position=position_kmer;
    nouvelle_position->suiv_pos= NULL;
    return;
  }
}
//------------------------------------------------------------------------------------------------------------
void parcours_masque(int longueur_masque, void* adr_masque, int nb_fenetre, int nb_sequence, ptr_struct_seq* adr_tete_struct_sequence, TPtr_Cellkmer* adr_tete_liste_kmer, TPtr_CellSequence* adr_tete_liste_sequence, TPtr_CellPos* adr_tete_liste_pos)
{
  ptr_struct_seq p_generation_seq = *adr_tete_struct_sequence;
  TPtr_Cellkmer p_kmer= *adr_tete_liste_kmer;
  TPtr_CellSequence p_sequence= *adr_tete_liste_sequence;
  TPtr_CellPos p_pos= *adr_tete_liste_pos;
  int *p_masque= adr_masque;
  int position=0;
  int cpt_pos_masque, position_kmer;
  int pos_kmer =0;
  char k_mer[nb_fenetre+1]; // le k_mer mesure la taille du nombre de fen�tre ouverte dans le masque
  while (p_generation_seq != NULL)
  {
    while (position<30)
    {
      for(cpt_pos_masque=0; cpt_pos_masque< longueur_masque; cpt_pos_masque++)// on parcourt les nucleotides sous le masque
      {
        if (p_masque[cpt_pos_masque]==1)//si la fen�tre du masque est ouverte
        {
          k_mer[pos_kmer]= p_generation_seq->sequence[position];
          if (pos_kmer==0)
          {
            position_kmer=position;
          }
          pos_kmer++;
          if (pos_kmer==nb_fenetre)
          {
            k_mer[pos_kmer]= '\0' ;
            pos_kmer=0;
          }
        }
        position= position +1;
      }
      generation_kmer(position_kmer, k_mer, &p_kmer, &p_sequence, &p_generation_seq, &p_pos);
    }
    position=0;
    p_generation_seq= p_generation_seq->next_sequence;
  }
  return;
}
//------------------------------------------------------------------------------------------------------------
//Procedure pour remplir le dictionnaire de kmer s�lectionn� pour le calcul de la PSSM:
void recuperer_motif_kmer(TPtr_Cellkmer* adr_parcours_kmer, TPtr_Cellkmer_selectionne *adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif_PSSM, TPtr_CellSequence* adr_cell_sequence, TPtr_CellPos* adr_cell_pos, ptr_struct_seq* adr_generation_sequence, int nb_sequence_kmer)
{
  TPtr_Cellkmer p_kmer= *adr_parcours_kmer;
  TPtr_Cellkmer_selectionne p_kmer_selectionne= *adr_tete_kmer_selectionne;
  ptr_struct_seq tete_generation_seq = *adr_generation_sequence;
  ptr_struct_seq p_generation_seq= tete_generation_seq;
  int longueur_sequence= 30;
  int longueur_motif=5;
  int i, sequence_actuelle;
  if (p_kmer_selectionne->suiv_kmer_selectionne == NULL) //la liste est vide:
  {
    TPtr_Cellkmer_selectionne prochain_kmer= malloc(sizeof(TCellkmer_selectionne)); //Je cr�e la prochaine brique
    TPtr_Cell_Motif_PSSM nouveau_motif= malloc(sizeof(TCell_Motif_PSSM));
    strcpy(p_kmer_selectionne->kmer, p_kmer->kmer);//Je remplis la brique actuelle
    p_kmer_selectionne->nb_sequence= nb_sequence_kmer;
    p_kmer_selectionne->suiv_kmer_selectionne=prochain_kmer; //je fais le chainage entre la brique actuelle et la prochaine
    p_kmer_selectionne->tete_motif_PSSM=nouveau_motif; //Je recupere la bonne tete
    TPtr_CellSequence p_sequence= p_kmer->tete_sequence;
    do
    {
      TPtr_CellPos p_pos= p_sequence->tete_pos; // je recupere la liste de position pour chaque sequence
      sequence_actuelle= p_sequence->sequence;
      while (p_generation_seq->numero_sequence != sequence_actuelle)
      { //je me place dans la bonne sequence
        p_generation_seq=p_generation_seq->next_sequence;
      }
      if((p_pos->position + longueur_motif)< longueur_sequence)
      {
        for (i=(p_pos->position); i<(p_pos->position + longueur_motif); i++)
        {
          nouveau_motif->motif[i-(p_pos->position)]=p_generation_seq->sequence[i];
        }
        nouveau_motif->motif[longueur_motif]= '\0';
        TPtr_Cell_Motif_PSSM prochain_motif= (TCell_Motif_PSSM*)malloc(sizeof(TCell_Motif_PSSM));
        nouveau_motif->suiv_motif=prochain_motif;
        nouveau_motif=prochain_motif;
      }
      p_sequence=p_sequence->suiv_sequence;
    }while(p_sequence != NULL);
    p_generation_seq= tete_generation_seq;
    return;
  }
  do //Cas o� on a plusieurs kmer selectionn�
  {
    p_kmer_selectionne=p_kmer_selectionne->suiv_kmer_selectionne;
  }while (p_kmer_selectionne->suiv_kmer_selectionne != NULL);

  TPtr_Cellkmer_selectionne prochain_kmer= (TCellkmer_selectionne*)malloc(sizeof(TCellkmer_selectionne));
  TPtr_Cell_Motif_PSSM nouveau_motif= (TCell_Motif_PSSM*)malloc(sizeof(TCell_Motif_PSSM));
  strcpy(p_kmer_selectionne->kmer, p_kmer->kmer);
  p_kmer_selectionne->nb_sequence= nb_sequence_kmer;
  p_kmer_selectionne->suiv_kmer_selectionne=prochain_kmer;
  p_kmer_selectionne->tete_motif_PSSM=nouveau_motif;
  TPtr_CellSequence p_sequence= p_kmer->tete_sequence;
  do
  {
    TPtr_CellPos p_pos= p_sequence->tete_pos; // je recupere la liste de position pour chaque sequence
    sequence_actuelle= p_sequence->sequence;
    while (p_generation_seq->numero_sequence != sequence_actuelle)
    {
      p_generation_seq=p_generation_seq->next_sequence;
    } //je me place dans la bonne sequence

    if((p_pos->position + longueur_motif)< longueur_sequence)
    {
      for (i=(p_pos->position); i<(p_pos->position + longueur_motif); i++)
      {
        nouveau_motif->motif[i-(p_pos->position)]=p_generation_seq->sequence[i];
      }
      nouveau_motif->motif[longueur_motif]= '\0';

      TPtr_Cell_Motif_PSSM prochain_motif= malloc(sizeof(TCell_Motif_PSSM));
      nouveau_motif->suiv_motif=prochain_motif;
      nouveau_motif=prochain_motif;
    }
    p_sequence=p_sequence->suiv_sequence;

  }while(p_sequence != NULL);
  p_generation_seq= tete_generation_seq;
  return;
}
//------------------------------------------------------------------------------------------------------------
//Fonction pour trouver si le kmer est trouv� dans chaque sequence:
void kmer_present_dans_chaque_sequence(int nb_sequence, TPtr_Cellkmer* adr_cell_kmer, TPtr_CellSequence *adr_cell_sequence, TPtr_CellPos *adr_cell_pos, ptr_struct_seq* adr_cell_generation_sequence, TPtr_Cellkmer_selectionne* adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_cell_motif_PSSM)
{
  TPtr_Cellkmer p_parcours_kmer= *adr_cell_kmer;
  TPtr_Cellkmer_selectionne tete_kmer_selectionne= *adr_cell_kmer_selectionne;
  TPtr_Cellkmer_selectionne p_kmer_selectionne= *adr_cell_kmer_selectionne;
  p_kmer_selectionne->suiv_kmer_selectionne= NULL;
  TPtr_Cell_Motif_PSSM p_motif_PSSM= *adr_cell_motif_PSSM;
  TPtr_CellSequence tete_sequence= *adr_cell_sequence;
  TPtr_CellPos tete_pos= *adr_cell_pos;
  ptr_struct_seq tete_generation_sequence= *adr_cell_generation_sequence;
  int nb_sequence_par_kmer=0;

  while (p_parcours_kmer != NULL)
  {
    TPtr_CellSequence p_parcours_sequence= p_parcours_kmer->tete_sequence;
    while (p_parcours_sequence != NULL )
    {
      nb_sequence_par_kmer++;
      p_parcours_sequence=p_parcours_sequence->suiv_sequence;
    }
    p_parcours_kmer->nb_sequence=nb_sequence_par_kmer;
    if (nb_sequence_par_kmer >= 7)
    { //Pour l"instant pour pouvoir continuer
      recuperer_motif_kmer(&p_parcours_kmer, &p_kmer_selectionne, &p_motif_PSSM, &tete_sequence, &tete_pos, &tete_generation_sequence, nb_sequence_par_kmer);
      p_kmer_selectionne= tete_kmer_selectionne;
    }
    nb_sequence_par_kmer=0;
    p_parcours_kmer=p_parcours_kmer->suiv_kmer;
  }
  return;
}
//------------------------------------------------------------------------------------------------------------
void affichage_dictionnaire_kmer(TPtr_Cellkmer* adr_tete_kmer, TPtr_CellSequence* adr_tete_sequence, TPtr_CellPos* adr_tete_pos)
{
  TPtr_Cellkmer p_kmer= *adr_tete_kmer;
  FILE* fichier_dictionnaire= NULL;
  fichier_dictionnaire= fopen("dictionnaire_kmer.txt", "w");
  while (p_kmer != NULL)
  {
    TPtr_CellSequence p_sequence=p_kmer->tete_sequence;
    fprintf(fichier_dictionnaire, "KMER: %s \n", p_kmer->kmer);
    while (p_sequence != NULL)
    {
      TPtr_CellPos p_pos= p_sequence->tete_pos;
      fprintf(fichier_dictionnaire, "Sequence: %d \n", p_sequence->sequence);
      while (p_pos != NULL)
      {
        fprintf(fichier_dictionnaire, "Position: %d \n ", p_pos->position);
        p_pos=p_pos->suiv_pos;
      }
      p_sequence=p_sequence->suiv_sequence;
    }
    fprintf(fichier_dictionnaire, "Nb de séquence dans lesquelles le kmer est présent: %d \n", p_kmer->nb_sequence);
    fprintf(fichier_dictionnaire, "\n\n\n");
    p_kmer=p_kmer->suiv_kmer;
  }
  fclose(fichier_dictionnaire);
  return;
}
//------------------------------------------------------------------------------------------------------------
void affichage_motif_selectionne(TPtr_Cellkmer_selectionne* adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif)
{
  TPtr_Cellkmer_selectionne p_kmer_selectionne= *adr_tete_kmer_selectionne;
  FILE* fichier_kmer_selectionne= NULL;
  fichier_kmer_selectionne= fopen("kmer_selectionne.txt", "w");
  while (p_kmer_selectionne != NULL)
  {
    TPtr_Cell_Motif_PSSM p_motif= p_kmer_selectionne->tete_motif_PSSM;
    fprintf(fichier_kmer_selectionne, "KMER: %s \n", p_kmer_selectionne->kmer);
    fprintf(fichier_kmer_selectionne, "Present dans %d sequences \n", p_kmer_selectionne->nb_sequence);
    while(p_motif != NULL)
    {
      fprintf(fichier_kmer_selectionne, "MOTIF: %s \n", p_motif->motif);
      p_motif=p_motif->suiv_motif;
    }
    p_kmer_selectionne=p_kmer_selectionne->suiv_kmer_selectionne;
  }
  fclose(fichier_kmer_selectionne);
  return;
}
//------------------------------------------------------------------------------------------------------------
void calcul_PSSM(TPtr_Cellkmer_selectionne *adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM *adr_cell_motif_PSSM, FILE** file_info, double*** adr_matrice_PSSM)
{
  *file_info= fopen("PSSM_Motif_Trouve.txt", "w"); // Dans ce fichier on va �crire la PSSM pr�visionnelle
  TPtr_Cellkmer_selectionne p_kmer_selectionne= *adr_cell_kmer_selectionne;
  double nb_sequence=p_kmer_selectionne->nb_sequence;
  int taille_motif=5;
  double add= 1/nb_sequence;
  char to_print;
  double maximum= -1;
  double** p_matrice_PSSM= *adr_matrice_PSSM;
  int i, j, a, t, c, g, k;
  //calcul de la PSSM � partir de la liste chain�e de motif
  TPtr_Cell_Motif_PSSM p_motif=p_kmer_selectionne->tete_motif_PSSM;
  do
  { // on remplit la PSSM
    for (i=0; i<6; i++)
    {
      switch(p_motif->motif[i])
      {
        case 'a':
            p_matrice_PSSM[0][i]= p_matrice_PSSM[0][i] + add;
            break;
        case 't':
            p_matrice_PSSM[1][i]= p_matrice_PSSM[1][i] + add;
            break;
        case 'c':
            p_matrice_PSSM[2][i]= p_matrice_PSSM[2][i] + add;
            break;
        case 'g':
            p_matrice_PSSM[3][i]= p_matrice_PSSM[3][i] + add;
            break;
      }
    }
    p_motif=p_motif->suiv_motif;
  }while (p_motif != NULL);

  //ecriture de la matrice PSSM dans le fichier d'info.
  fprintf(*file_info, "\n\nPSSM: \n");
  fprintf(*file_info, "a  ");
  for (a=0; a<taille_motif; a++)
  {
    fprintf(*file_info, "%.2f ", p_matrice_PSSM[0][a]);
  }
  fprintf(*file_info, "\nt  ");
  for (t=0; t<taille_motif; t++)
  {
    fprintf(*file_info, "%.2f ", p_matrice_PSSM[1][t]);
  }
  fprintf(*file_info, "\nc  ");
  for (c=0; c<taille_motif; c++)
  {
    fprintf(*file_info, "%.2f ", p_matrice_PSSM[2][c]);
  }
  fprintf(*file_info, "\ng  ");
  for (g=0; g<taille_motif; g++)
  {
    fprintf(*file_info, "%.2f ", p_matrice_PSSM[3][g]);
  }

  //ecriture du motif consensus dans le fichier a partir de la matrice PSSM:
  fprintf(*file_info, "\n\nMotif Consensus: \n" );
  for (j=0; j<taille_motif; j++)
  {
    for(k=0; k<4; k++)
    {
      if (p_matrice_PSSM[k][j]> maximum)
      {
        maximum=p_matrice_PSSM[k][j];
        switch (k)
        {
          case 0:
              to_print= 'a';
              break;
          case 1:
              to_print= 't';
              break;
          case 2:
              to_print= 'c';
              break;
          case 3:
              to_print= 'g';
              break;
        }
      }
    }
    maximum=-1; //on remet maximum a -1 avant d'�valuer la nucl�otide majoritaire de la s�quence suivante !
    fprintf(*file_info, "%c", to_print);
  }
  fprintf(*file_info, "\n\n\n");
  fclose(*file_info);
  return;
}
//------------------------------------------------------------------------------------------------------------
void calcul_nouvelle_PSSM(TPtr_Cell_Motif_PSSM *adr_cell_mot_selected, double*** adr_matrice_PSSM_nouv, double nb_sequence, char (*adr_Ct)[6])
{
  FILE* file_nouv_PSSM;
  file_nouv_PSSM= fopen("Nouvelle_PSSM_Motif_Trouve.txt", "a"); // Dans ce fichier on va �crire la PSSM pr�visionnelle
  int taille_motif=5;
  double add= 1/nb_sequence;
  TPtr_Cell_Motif_PSSM p_mot= *adr_cell_mot_selected;
  char to_print;
  double maximum= -1;
  double** p_matrice_PSSM_nouv= *adr_matrice_PSSM_nouv;
  int i, j, a, t, c, g, k;
  //calcul de la PSSM � partir de la liste chain�e de motif
  do
  { // on remplit la PSSM
    for (i=0; i<6; i++)
    {
      switch(p_mot->motif[i])
      {
      case 'a':
          p_matrice_PSSM_nouv[0][i]= p_matrice_PSSM_nouv[0][i] + add;
          break;
      case 't':
          p_matrice_PSSM_nouv[1][i]= p_matrice_PSSM_nouv[1][i] + add;
          break;
      case 'c':
          p_matrice_PSSM_nouv[2][i]= p_matrice_PSSM_nouv[2][i] + add;
          break;
      case 'g':
          p_matrice_PSSM_nouv[3][i]= p_matrice_PSSM_nouv[3][i] + add;
          break;
      }
    }
    p_mot=p_mot->suiv_motif;
  }while (p_mot != NULL);

  //ecriture de la matrice PSSM dans le fichier d'info.
  fprintf(file_nouv_PSSM, "\n\nPSSM: \n");
  fprintf(file_nouv_PSSM, "a  ");
  for (a=0; a<taille_motif; a++)
  {
    fprintf(file_nouv_PSSM, "%.2f ", p_matrice_PSSM_nouv[0][a]);
  }
  fprintf(file_nouv_PSSM, "\nt  ");
  for (t=0; t<taille_motif; t++)
  {
    fprintf(file_nouv_PSSM, "%.2f ", p_matrice_PSSM_nouv[1][t]);
  }
  fprintf(file_nouv_PSSM, "\nc  ");
  for (c=0; c<taille_motif; c++)
  {
    fprintf(file_nouv_PSSM, "%.2f ", p_matrice_PSSM_nouv[2][c]);
  }
  fprintf(file_nouv_PSSM, "\ng  ");
  for (g=0; g<taille_motif; g++)
  {
    fprintf(file_nouv_PSSM, "%.2f ", p_matrice_PSSM_nouv[3][g]);
  }
  //ecriture du motif consensus dans le fichier a partir de la matrice PSSM:
  fprintf(file_nouv_PSSM, "\n\nMotif Consensus: \n" );
  for (j=0; j<taille_motif; j++)
  {
    for(k=0; k<4; k++)
    {
      if (p_matrice_PSSM_nouv[k][j]> maximum){
        maximum=p_matrice_PSSM_nouv[k][j];
        switch (k)
        {
          case 0:
              to_print= 'a';
              break;
          case 1:
              to_print= 't';
              break;
          case 2:
              to_print= 'c';
              break;
          case 3:
              to_print= 'g';
              break;
        }
      }
    }
    *adr_Ct[j]=to_print;
    printf("MOTIF CONSENSUS: %c", *adr_Ct[j]);
    maximum=-1; //on remet maximum a -1 avant d'�valuer la nucl�otide majoritaire de la s�quence suivante !
    fprintf(file_nouv_PSSM, "%c", to_print);
  }
  *adr_Ct[taille_motif+1]='\0';
  fprintf(file_nouv_PSSM, "\n\n\n");
  fclose(file_nouv_PSSM);
  return;
}
//------------------------------------------------------------------------------------------------------------
void calcul_score(TPtr_Mot_Ameliorer_PSSM* adr_mot, double*** adr_matrice_PSSM, int n_sequence, ptr_struct_seq* adr_generation_sequence, int longueur_masque)
{
  ptr_struct_seq p_generation_seq= *adr_generation_sequence;
  TPtr_Mot_Ameliorer_PSSM p_mot= *adr_mot;
  double** p_matrice_PSSM= *adr_matrice_PSSM;
  int position=0;
  int i;
  double nb_a=0;
  double nb_t=0;
  double nb_c=0;
  double nb_g=0;
  double prob_mot_background=0; //probabilite du mot selon le background P(M|B)
  double prob_mot_PSSM=0; // probabilite du mot selon la PSSM P(M|PSSM)
  double score; // P(M|PSSM)/P(M|B)
  //CALCUL DU BACKGROUND:
  //On vérifie qu'on est dans la bonne séquence:
  if (p_generation_seq->numero_sequence !=  n_sequence)
  {
    p_generation_seq=p_generation_seq->next_sequence;
  }
  while (position<30)
  {
    switch(p_generation_seq->sequence[position])
    {
      case 'a':
      nb_a++;
      break;
      case 't':
      nb_t++;
      break;
      case 'c':
      nb_c++;
      break;
      case 'g':
      nb_g++;
    }
    position++;
  }
  nb_a= nb_a/30;
  nb_t= nb_t/30;
  nb_c= nb_c/30;
  nb_g= nb_g/30;
  position=0;

  for (i=0; i<= longueur_masque; i++)
  {
    switch(p_mot->mot[i])
    {
      case 'a':
        prob_mot_background= prob_mot_background + nb_a;
        prob_mot_PSSM= prob_mot_PSSM + p_matrice_PSSM[0][i];
        break;
      case 't':
        prob_mot_background= prob_mot_background + nb_t;
        prob_mot_PSSM= prob_mot_PSSM + p_matrice_PSSM[1][i];
        break;
      case 'c':
        prob_mot_background= prob_mot_background + nb_c;
        prob_mot_PSSM= prob_mot_PSSM + p_matrice_PSSM[2][i];
        break;
      case 'g':
        prob_mot_background= prob_mot_background + nb_g;
        prob_mot_PSSM= prob_mot_PSSM + p_matrice_PSSM[3][i];
        break;
    }
  }
  score= prob_mot_PSSM/prob_mot_background;
  p_mot->score_mot= score;
  return;
}
//------------------------------------------------------------------------------------------------------------
double dist_PSSM(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv, double* distance_PSSM)
{
	double** p_ancienne_matrice_PSSM= *adr_matrice_PSSM;
	double** p_nouvelle_matrice_PSSM= *adr_matrice_PSSM_nouv;
	int i, j;
	double somme_distance=0;
	for (i=0; i< 4; i++)
	{
		for (j=0; j<6; j++)
		{
			somme_distance= somme_distance + fabs(p_ancienne_matrice_PSSM[i][j] - p_nouvelle_matrice_PSSM[i][j]);
		}
	}
	*distance_PSSM= somme_distance;
	printf("Somme Distance= %f \n", somme_distance);
	return(somme_distance);
}
//------------------------------------------------------------------------------------------------------------
int distanceHammingSt1(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st1)
{
	int i;
	int Dh=0;
  int st1=0;
	TPtr_Cell_Motif_PSSM p_mot_selected= *adr_mot_selected;
	Ptr_st p_st1= *adr_st1;

	while (p_mot_selected !=NULL)
	{
		for (i=0; i<5; i++)
		{
			if (*adr_Ct[i] != p_mot_selected->motif[i])
			{
				Dh++;
				printf("Dh: %d \n", Dh);
			}
		}
    if (Dh>1)
		{
		  st1++;
		}
    strcpy(p_st1->mot, p_mot_selected->motif);
    printf("p_st1->mot: %s \n", p_st1->mot);
    p_st1->distance_hamming=Dh;
    printf("p_st1->distance: %d \n", p_st1->distance_hamming);
    Ptr_st p_st1_suiv= malloc(sizeof(st1));
    p_st1->next_mot= p_st1_suiv;
    p_st1=p_st1_suiv;
		Dh=0;
		p_mot_selected=p_mot_selected->suiv_motif;
	}
	return(st1-1); // Parce que le dernier mot NULL de la liste chainée est compté
}
//------------------------------------------------------------------------------------------------------------
int distanceHammingSt2(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st2)
{
  printf("DISTANCE HAMMING ST2");
  int st2=0;
	int i;
	int Dh=0;
	TPtr_Cell_Motif_PSSM p_mot_selected= *adr_mot_selected;
	Ptr_st p_st2= *adr_st2;
	printf("IN ST2");
	while (p_mot_selected !=NULL)
	{
		for (i=0; i<5; i++)
		{
			if (*adr_Ct[i] != p_mot_selected->motif[i])
			{
				Dh++;
			}
		}
    strcpy(p_st2->mot, p_mot_selected->motif);
    printf("p_st2->mot: %s \n", p_st2->mot);
    p_st2->distance_hamming=Dh;
    printf("p_st2->distance: %d \n", p_st2->distance_hamming);
    Ptr_st p_st2_suiv= malloc(sizeof(st2));
    p_st2->next_mot= p_st2_suiv;
    p_st2=p_st2_suiv;
		if (Dh<2)
		{
      st2++;
		}
		Dh=0;
		p_mot_selected=p_mot_selected->suiv_motif;
	}
	return(st2);
}
//------------------------------------------------------------------------------------------------------------
int distanceHammingSt2_prim(char (*adr_Ct)[6], ptr_struct_seq* adr_generation_sequence, TPtr_Mot_Ameliorer_PSSM *adr_mot, Ptr_st* adr_st2_prim)
{
	printf("IN ST2 PRIM");
	int i;
  int st2_prim=0;
  int position=0;
	int Dh_min=8;
	int Dh=0;
	ptr_struct_seq p_generation_seq= *adr_generation_sequence;
	Ptr_st p_st2_prim= *adr_st2_prim;
  TPtr_Mot_Ameliorer_PSSM p_mot= *adr_mot;

	while (p_generation_seq !=NULL)
	{
		while (position<=25)
		{
			for (i=0; i<5; i++)
			{
        p_mot->mot[i]=p_generation_seq->sequence[position];
				if (*adr_Ct[i] != p_mot->mot[i])
				{
					Dh++;
				}
        if (i==4){
          p_mot->mot[5]= '\0';
        }
        position++;
			}
      position= position-4;
			if (Dh<=Dh_min)
			{
        Dh_min=Dh;
				strcpy(p_st2_prim->mot, p_mot->mot);
				p_st2_prim->distance_hamming=Dh;
			}
			Dh=0;
      TPtr_Mot_Ameliorer_PSSM p_nouv_mot=malloc(sizeof(TMot_Ameliorer_PSSM));
      p_mot->next_mot= p_nouv_mot;
      p_mot=p_nouv_mot;
    }
    position=0;
    Dh_min=8;
    p_generation_seq=p_generation_seq->next_sequence;
		printf("p_st2->mot (minimum): %s \n", p_st2_prim->mot);
		printf("p_st2->distance (minimum): %d \n", p_st2_prim->distance_hamming);
		Ptr_st p_st2_prim_suiv= malloc(sizeof(st));
		p_st2_prim->next_mot= p_st2_prim_suiv;
		p_st2_prim=p_st2_prim_suiv;
	}
  p_st2_prim= *adr_st2_prim;
  while (p_st2_prim != NULL)
  {
    if(p_st2_prim->distance_hamming < 2)
    {
      st2_prim++;
    }
    p_st2_prim=p_st2_prim->next_mot;
  }
	return(st2_prim-1); //On met -1 car la boucle précedente prend en compte la case NULL de la fin de liste chainé
}
//------------------------------------------------------------------------------------------------------------
void trier(int* v_St1_Dh, int* v_St1_Pos, int g, int d)
{
  int indice_pivot;
  if (g < d)
  {
    separer(v_St1_Dh, v_St1_Pos, g, d, &indice_pivot);
    indice_pivot --;
    separer(v_St1_Dh, v_St1_Pos, g, d, &indice_pivot);
    trier(v_St1_Dh, v_St1_Pos, indice_pivot+1, d);
  }
  return;
}
//--------------------------------------------------------------------------------------------------
void separer(int* v_St1_Dh, int* v_St1_Pos, int g, int d, int* adr_indice_pivot)
{
  //PRECONDITIONS
  // g < d
  int bas, haut; //indices de position dans le vecteur
  int comp, pivot;  //comp pour comparateur, c'est variables sont des entiers car on trie un vecteur d'entiers
  int comp_pos, pivot_pos;

  bas = g;
  haut = d;
  pivot = v_St1_Dh[bas];
  pivot_pos= v_St1_Pos[bas];
  comp = v_St1_Dh[haut];
  comp_pos=v_St1_Pos[haut];

  while (bas < haut)
  {
    if ( comp > pivot)
    {
      v_St1_Dh[haut] = comp;
      v_St1_Pos[haut]= comp_pos;
      haut --;
      comp = v_St1_Dh[haut];
      comp_pos=v_St1_Pos[haut];
    }else{
      v_St1_Dh[bas] = comp;
      v_St1_Pos[bas]= comp_pos;
      bas ++;
      comp = v_St1_Dh[bas];
      comp_pos=v_St1_Pos[bas];
    }
  }
  v_St1_Dh[bas] = pivot;
  v_St1_Pos[bas]= pivot_pos;
  *adr_indice_pivot = bas;
  return;
}
//------------------------------------------------------------------------------------------------------------
void quick_sort_ST(Ptr_st* adr_st1, int* v_St1_Pos)
{
  Ptr_st p_st1= *adr_st1;
  int i;
  int v_St1_Dh[9]; //Distance de Hamming de chaque occurrence
  int borne_gauche = 0;
  int borne_droite = 9;
  for (i=0; i<10; i++)
  {
    v_St1_Dh[i]= p_st1->distance_hamming;
    v_St1_Pos[i]= i;
    printf("Vecteur Distance Hamming: %d \n", v_St1_Dh[i]);
    printf("Vecteur Pos: %d \n", v_St1_Pos[i]);
    p_st1=p_st1->next_mot;
  }
  printf("\n");

  trier(v_St1_Dh, v_St1_Pos, borne_gauche, borne_droite);

  for(i=0; i<10; i++)
  {
    printf("Vecteur Distance Hamming après tri: %d \n", v_St1_Dh[i]);
    printf("Vecteur Pos après tri: %d \n", v_St1_Pos[i]);
  }
  printf("\n");
  return;
}
//--------------------------------------------------------------------------------------------------
void fichier_sortie_st(Ptr_st* adr_st1, int* v_St1_Pos, char (*adr_Ct)[6])
{
  Ptr_st p_st1= *adr_st1;
  int i, position, cpt;
  FILE* fichier_sortie_st2;
  fichier_sortie_st2= fopen("score_ST2.txt", "w");

  for (i=0; i<10; i++)
  {
    p_st1=*adr_st1;
    position= v_St1_Pos[i];
    for (cpt=0; cpt<position; cpt++)
    {
      p_st1=p_st1->next_mot;
    }
    fprintf(fichier_sortie_st2, "OCCURRENCE: %s --> SCORE: %d \n", p_st1->mot, p_st1->distance_hamming);
  }
  fprintf(fichier_sortie_st2, "\n\nMOTIF CONSENSUS: \n");
  for (i=0; i<5; i++)
  {
    fprintf(fichier_sortie_st2, "%c", *adr_Ct[i]);
  }
  return;
}
