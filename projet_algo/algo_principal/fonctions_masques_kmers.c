
//------------------------------------------------------------------------------------------------------------
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
void generation_kmer(int position_kmer, char* k_mer, TPtr_Cellkmer* adr_liste_kmer, TPtr_CellSequence* adr_liste_sequence, TPtr_ensemble_sequences* ptr_ensemble, TPtr_CellPos* adr_liste_pos)
{
  TPtr_Cellkmer temp_p_kmer= *adr_liste_kmer; //on cr�� des pointeurs temporaire pour parcourir les listes
  TPtr_ensemble_sequences p_dictionnaire_sequence= *ptr_ensemble;

  if (temp_p_kmer->suiv_kmer==NULL) //cas du premier element de la liste
  {
    TPtr_Cellkmer nouveau_kmer=(TCellkmer*)malloc(sizeof(TCellkmer));
    strcpy(temp_p_kmer->kmer, k_mer);
    nouveau_kmer->suiv_kmer= NULL;
    TPtr_CellSequence nouvelle_sequence=(TCellSequence*)malloc(sizeof(TCellSequence));
    temp_p_kmer->suiv_kmer= nouveau_kmer;
    temp_p_kmer->tete_sequence=nouvelle_sequence;
    nouvelle_sequence->sequence=p_dictionnaire_sequence->numero_sequence;
    nouvelle_sequence->suiv_sequence= NULL;
    TPtr_CellPos nouvelle_position=(TCellPos*)malloc(sizeof(TCellPos));
    nouvelle_sequence->tete_pos= nouvelle_position;
    nouvelle_position->position=position_kmer;
    nouvelle_position->suiv_pos= NULL;
    return;
  }else{
    while (temp_p_kmer->suiv_kmer != NULL)
    {
      if (strcmp(temp_p_kmer->kmer, k_mer)==0)   //si le kmer a d�j� �t� trouv�
      {
        TPtr_CellSequence p_liste_sequence=temp_p_kmer->tete_sequence;
        //premier element de la liste chainee de sequence:
        while (p_liste_sequence->suiv_sequence != NULL)
        {
          if (p_dictionnaire_sequence->numero_sequence==p_liste_sequence->sequence) //Si le kmer a d�j� �t� trouv� dans cette s�quence
          {
            TPtr_CellPos p_liste_pos= p_liste_sequence->tete_pos;
            while (p_liste_pos->suiv_pos != NULL){
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
        if (p_dictionnaire_sequence->numero_sequence==p_liste_sequence->sequence)  //Si le kmer a d�j� �t� trouv� dans cette s�quence
        {
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
        nouvelle_sequence->sequence= p_dictionnaire_sequence->numero_sequence;
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
    nouvelle_sequence->sequence=p_dictionnaire_sequence->numero_sequence;
    nouvelle_sequence->suiv_sequence= NULL;
    TPtr_CellPos nouvelle_position=(TCellPos*)malloc(sizeof(TCellPos));
    nouvelle_sequence->tete_pos= nouvelle_position;
    nouvelle_position->position=position_kmer;
    nouvelle_position->suiv_pos= NULL;
    return;
  }
}
//------------------------------------------------------------------------------------------------------------
void parcours_masque(int longueur_masque, void* adr_masque, int nb_fenetre, int nb_sequence, TPtr_ensemble_sequences* ptr_ensemble, TPtr_Cellkmer* adr_tete_liste_kmer, TPtr_CellSequence* adr_tete_liste_sequence, TPtr_CellPos* adr_tete_liste_pos)
{
	TPtr_ensemble_sequences p_dictionnaire_sequence= *ptr_ensemble;
  TPtr_Cellkmer p_kmer= *adr_tete_liste_kmer;
  TPtr_CellSequence p_sequence= *adr_tete_liste_sequence;
  TPtr_CellPos p_pos= *adr_tete_liste_pos;
  int *p_masque= adr_masque;
  int position=0;
  int cpt_pos_masque, position_kmer;
  int pos_kmer =0;
  char k_mer[nb_fenetre+1]; // le k_mer mesure la taille du nombre de fen�tre ouverte dans le masque

  while (p_dictionnaire_sequence != NULL)
  {
    while (position<30)
    {
      for(cpt_pos_masque=0; cpt_pos_masque< longueur_masque; cpt_pos_masque++)
      { // on parcourt les nucleotides sous le masque
        if (p_masque[cpt_pos_masque]==1)
        { //si la fen�tre du masque est ouverte
          k_mer[pos_kmer]= p_dictionnaire_sequence->sequence[position];

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
      generation_kmer(position_kmer, k_mer, &p_kmer, &p_sequence, &p_dictionnaire_sequence, &p_pos);
    }
    position=0;
    p_dictionnaire_sequence= p_dictionnaire_sequence->suiv_seq;
  }
  return;
}
//------------------------------------------------------------------------------------------------------------
//Procedure pour remplir le dictionnaire de kmer s�lectionn� pour le calcul de la PSSM:
void recuperer_motif_kmer(TPtr_Cellkmer* adr_parcours_kmer, TPtr_Cellkmer_selectionne *adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif_PSSM, TPtr_CellSequence* adr_cell_sequence, TPtr_CellPos* adr_cell_pos, TPtr_ensemble_sequences* ptr_ensemble, int nb_sequence_kmer)
{
  TPtr_Cellkmer p_kmer= *adr_parcours_kmer;
  TPtr_Cellkmer_selectionne p_kmer_selectionne= *adr_tete_kmer_selectionne;
  TPtr_ensemble_sequences tete_dictionnaire_sequence = *ptr_ensemble;
  TPtr_ensemble_sequences p_dictionnaire_sequence= tete_dictionnaire_sequence;
  int longueur_sequence= 30;
  int longueur_motif=5;
  int i, sequence_actuelle;

  if (p_kmer_selectionne->suiv_kmer_selectionne == NULL) //la liste est vide:
  {
    TPtr_Cellkmer_selectionne prochain_kmer= malloc(sizeof(TCellkmer_selectionne)); //Je créer la prochaine brique
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
      while (p_dictionnaire_sequence->numero_sequence != sequence_actuelle)  //je me place dans la bonne sequence
      {
        p_dictionnaire_sequence=p_dictionnaire_sequence->suiv_seq;
      }
      if((p_pos->position + longueur_motif)< longueur_sequence)
      {
        for (i=(p_pos->position); i<(p_pos->position + longueur_motif); i++)
        {
          nouveau_motif->motif[i-(p_pos->position)]=p_dictionnaire_sequence->sequence[i];
        }
        nouveau_motif->motif[longueur_motif]= '\0';

        TPtr_Cell_Motif_PSSM prochain_motif= (TCell_Motif_PSSM*)malloc(sizeof(TCell_Motif_PSSM));
        nouveau_motif->suiv_motif=prochain_motif;
        nouveau_motif=prochain_motif;
      }
      p_sequence=p_sequence->suiv_sequence;
    }while(p_sequence != NULL);
    p_dictionnaire_sequence= tete_dictionnaire_sequence;
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
      while (p_dictionnaire_sequence->numero_sequence != sequence_actuelle)
      {
        p_dictionnaire_sequence=p_dictionnaire_sequence->suiv_seq;
      } //je me place dans la bonne sequence

      if((p_pos->position + longueur_motif)< longueur_sequence)
      {
        for (i=(p_pos->position); i<(p_pos->position + longueur_motif); i++)
        {
          nouveau_motif->motif[i-(p_pos->position)]=p_dictionnaire_sequence->sequence[i];
        }
        nouveau_motif->motif[longueur_motif]= '\0';

        TPtr_Cell_Motif_PSSM prochain_motif= malloc(sizeof(TCell_Motif_PSSM));
        nouveau_motif->suiv_motif=prochain_motif;
        nouveau_motif=prochain_motif;
      }
      p_sequence=p_sequence->suiv_sequence;
    }while(p_sequence != NULL);

  p_dictionnaire_sequence= tete_dictionnaire_sequence;
  return;
}
//------------------------------------------------------------------------------------------------------------
//Fonction pour trouver si le kmer est trouv� dans chaque sequence:
void kmer_present_dans_chaque_sequence(int nb_sequence, TPtr_Cellkmer* adr_cell_kmer, TPtr_CellSequence *adr_cell_sequence, TPtr_CellPos *adr_cell_pos, TPtr_ensemble_sequences* ptr_ensemble, TPtr_Cellkmer_selectionne* adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_cell_motif_PSSM)
{
  TPtr_Cellkmer p_parcours_kmer= *adr_cell_kmer;
  TPtr_Cellkmer_selectionne tete_kmer_selectionne= *adr_cell_kmer_selectionne;
  TPtr_Cellkmer_selectionne p_kmer_selectionne= *adr_cell_kmer_selectionne;
  p_kmer_selectionne->suiv_kmer_selectionne= NULL;
  TPtr_Cell_Motif_PSSM p_motif_PSSM= *adr_cell_motif_PSSM;
  TPtr_CellSequence tete_sequence= *adr_cell_sequence;
  TPtr_CellPos tete_pos= *adr_cell_pos;
  TPtr_ensemble_sequences tete_dictionnaire_sequence = *ptr_ensemble;
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
    if (nb_sequence_par_kmer >= 7) //Pour l"instant pour pouvoir continuer
    {
      recuperer_motif_kmer(&p_parcours_kmer, &p_kmer_selectionne, &p_motif_PSSM, &tete_sequence, &tete_pos, &tete_dictionnaire_sequence, nb_sequence_par_kmer);
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
