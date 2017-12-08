#ifndef HEADER_FONCTIONS
#define HEADER_FONCTIONS
#include <stdio.h>

  #define TAILLE_MAX_SEQ 800
/*************************************************************************************************************
 * * *                                       VARIABLES GLOBALES                                          * * *
 *************************************************************************************************************/
  extern bool convergence;
  extern int taille_motif;
  extern int d;
  extern int nb_fenetres;
  extern int nb_masques;
/*************************************************************************************************************
 * * *                                         STRUCTURES                                                * * *
 *************************************************************************************************************/
  typedef struct TDictionnaire_Sequences
  {
    int numero_sequence;
    char nom_seq[120];
    char sequence[TAILLE_MAX_SEQ];
    struct TDictionnaire_Sequences* suiv_seq;
  } TDictionnaire_Sequences;
  typedef TDictionnaire_Sequences* TPtr_dictionnaire_sequences;

  typedef struct TInfo_dictionnaire_sequences
  {
    int nb_sequences;
    TPtr_dictionnaire_sequences tete_dictionnaire_seq; //tete pointant sur la liste chainée de sequences
  } TInfo_dictionnaire_sequences;
  typedef TInfo_dictionnaire_sequences* TPtr_info_dictionnaire_sequences;

/*******************************
 * STRUCTURE DICTIONNAIRE KMER *
 *******************************/


  //1er element: liste chainee des kmers
  typedef struct TCellkmer
  {
    struct TCellkmer* suiv_kmer; //Pointeur sur le kmer suivant
    struct TCellSequence* tete_sequence; //Pointer qui pointe sur une liste chain�e de s�quence o� le kmer est retrouv�
    double nb_sequence; // nombre de séquence dans lesquelles le kmer est présent
    char* kmer; //"TCC" par exemple
  } TCellkmer;
  typedef TCellkmer* TPtr_Cellkmer; //Pointeur sur Tcellkmer

  //deuxieme element du dictionnaire kmer: liste chainee de sequence possedant le kmer:
  typedef struct TCellSequence
  {
    int num_sequence; //correspond au numéro de la séquence dans le dictionnaire de séquences
    struct TCellSequence* suiv_sequence; // Pointeur qui pointe l'element suivant de la liste chainee de sequence
    struct TCellPos* tete_pos; // Pointeur qui pointe sur le premier element de la structure position qui repertorie toutes les positions o� le kmer a �t� trouv� dans une s�quence
  } TCellSequence;
  typedef TCellSequence* TPtr_CellSequence; // Pointeur sur TCellSequence

  // troisieme element du dictionnaire de kmer: liste chainee de position ou le kmer a ete trouve dans une sequence.
  typedef struct TCellPos
  {
    int position; // position du kmer dans la s�quence
    struct TCellPos* suiv_pos; //pointe sur la position suivante
  } TCellPos;
  typedef TCellPos* TPtr_CellPos;

/******************************
 * STRUCTURE KMER SELECTIONNE *
 ******************************/
  //1er element: liste chainee des kmers
  typedef struct TCellkmer_selectionne
  {
    struct TCellkmer_selectionne* suiv_kmer_selectionne;
    struct TCell_Motif_PSSM* tete_motif_PSSM;
    double nb_sequence;
    char* kmer;
  } TCellkmer_selectionne;
  typedef TCellkmer_selectionne* TPtr_Cellkmer_selectionne;

  //2eme element: liste chainee de motif pour lequels on va calculer la PSSM:
  typedef struct TCell_Motif_PSSM
  {
    struct TCell_Motif_PSSM* suiv_motif;
    char motif[];
  } TCell_Motif_PSSM;
  typedef TCell_Motif_PSSM* TPtr_Cell_Motif_PSSM;

  typedef struct TInfo_dictionnaire_kmer
  {
    int nb_kmer;
    TPtr_Cellkmer tete_liste_kmer;
    TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne;
  } TInfo_dictionnaire_kmer;
  typedef TInfo_dictionnaire_kmer* TPtr_info_dictionnaire_kmer;

/****************************************************
 * STRUCTURE CHAINEE DE MOTS POUR AMELIORER LA PSSM *
 ****************************************************/
  typedef struct TMot_Ameliorer_PSSM
  {
    double score_mot; //score du mot
    struct TMot_Ameliorer_PSSM * next_mot; //pointeur sur mot suivant
    char mot[]; //Mot de longueur motif
  } TMot_Ameliorer_PSSM;
  typedef TMot_Ameliorer_PSSM* TPtr_Mot_Ameliorer_PSSM;

/************************************************
 * STRUCTURE CHAINEE ST1/ST2 CONTENANT LES MOTS *
 ************************************************/
  typedef struct st
  {
    // position /s si meme motif pour plusieurs séquences
    int distance_hamming; //distance de hamming entre le mot et Ct
    struct st * next_mot; //pointeur sur mot suivant
    char mot[]; //Mot qui ont une distance de hamming supérieur à 2
  } st;
  typedef st* Ptr_st;

/*************************************************************************************************************
 * * *                                           FONCTIONS                                               * * *
 *************************************************************************************************************/
  void importer_parametres(char* nom_fichier, int* taille_motif, int* d, int* nb_fenetres, int* nb_masques);
  void importer_sequences_fasta(char* nom_fichier_fasta, TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);
  void afficher_sequences(TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);
  void liberation_dictionnaire_sequence( TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);

  int random_number(int max_number, int zero_excluded);

  int* generation_masque(int* masque);
  void parcours_masque( int* masque, TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void generation_dictionnaire_kmer(int position_kmer, char* k_mer, TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void affichage_dictionnaire_kmer(TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void liberation_dictionnaire_kmer(TPtr_info_dictionnaire_kmer tete_info_dict_kmer);

  // void recuperer_motif_kmer(TPtr_Cellkmer* adr_parcours_kmer, TPtr_Cellkmer_selectionne *adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif_PSSM, TPtr_CellSequence* adr_cell_sequence, TPtr_CellPos* adr_cell_pos, TPtr_dictionnaire_sequences* ptr_ensemble, int nb_sequence_kmer);
  // void kmer_present_dans_chaque_sequence(int nb_sequences, TPtr_info_dictionnaire_kmer tete_info_dict_kmer);

  // void affichage_motif_selectionne(TPtr_Cellkmer_selectionne* adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif);
  //
  void creation_PSSM(double*** adr_matrice_PSSM);
  // void calcul_PSSM(TPtr_Cellkmer_selectionne *adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM *adr_cell_motif_PSSM, double*** adr_matrice_PSSM, int taille_motif);
  // void calcul_nouvelle_PSSM(TPtr_Cell_Motif_PSSM *adr_cell_mot_selected, double*** adr_matrice_PSSM, double nb_sequence, char (*adr_Ct)[6], int taille_motif);
  void afficher_PSSM( double** matrice_PSSM);
  void liberation_PSSM(double*** adr_matrice_PSSM);
  //
  // void calcul_score(TPtr_Mot_Ameliorer_PSSM* adr_mot, double*** adr_matrice_PSSM, int n_sequence, TPtr_dictionnaire_sequences* ptr_ensemble, int taille_motif);
  // double dist_PSSM(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv, double* distance_PSSM, int taille_motif);
  //
  // int distanceHammingSt1(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st1, int taille_motif);
  // int distanceHammingSt2(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st2, int taille_motif);
  // int distanceHammingSt2_prim(char (*adr_Ct)[6], TPtr_dictionnaire_sequences* ptr_ensemble, TPtr_Mot_Ameliorer_PSSM *adr_mot, Ptr_st* adr_st2_prim, int taille_motif);
  //
  // void quick_sort_ST(Ptr_st* adr_st1, int* v_St1_Pos, int nb_sequences);
  // void trier(int* v_St1_Dh, int* v_St1_Pos, int g, int d);
  // void separer(int* v_St1_Dh, int* v_St1_Pos, int g, int d, int* adr_indice_pivot);
  //
  // void fichier_sortie_st(Ptr_st* adr_st1, int* v_St1_Pos, char (*adr_Ct)[6], int nb_sequences);

#endif
