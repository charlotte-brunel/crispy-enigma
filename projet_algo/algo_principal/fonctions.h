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
    int numero_sequence;		// // kmer_present_dans_chaque_sequence(tete_info_dict_seq->nb_sequences, &tete_liste_kmer2,, &tete_liste_pour_recup_motif, &tete_liste_kmer_selectionne, &tete_liste_motif_PSSM);
    // //
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
    char* motif; //on stocke le motif équivalent au kmer
  } TCellPos;
  typedef TCellPos* TPtr_CellPos;

  // quatrième élément contenant tous les motifs correspondants aux kmers trouvés
  typedef struct TCell_Motif
  {
    struct TCell_Motif* suiv_motif;
    char* motif;
  } TCell_Motif;
  typedef TCell_Motif* TPtr_Cell_Motif;

  // header contenant les tetes des elements précédents
  typedef struct TInfo_dictionnaire_kmer
  {
    int nb_kmer;
    TPtr_Cellkmer tete_liste_kmer;
    TPtr_Cellkmer_selectionne tete_liste_kmer_selectionne;
    TPtr_Cell_Motif tete_liste_motif;
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
  void affichage_dictionnaire_sequences(TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);
  void liberation_dictionnaire_sequences( TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);

  int random_number(int max_number, int zero_excluded);
  int* generation_masque(int* masque);
  void parcours_masque_sur_seq( int* masque, TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_kmer tete_info_dict_kmer);

  void generation_dictionnaire_kmers(int position_kmer, char* k_mer, char* motif, TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void affichage_dictionnaire_kmers(TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void liberation_dictionnaire_kmers(TPtr_info_dictionnaire_kmer tete_info_dict_kmer);

  void creation_liste_motifs( char* motif ,TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void afficher_liste_motifs(TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void liberation_liste_motifs(TPtr_info_dictionnaire_kmer tete_info_dict_kmer);

  void creation_PSSM(double*** adr_matrice_PSSM);
  void afficher_PSSM( double** matrice_PSSM);
  void liberation_PSSM(double*** adr_matrice_PSSM);
  // void calcul_PSSM(TPtr_Cellkmer_selectionne *adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM *adr_cell_motif_PSSM, double*** adr_matrice_PSSM, int taille_motif);
  // void calcul_nouvelle_PSSM(TPtr_Cell_Motif_PSSM *adr_cell_mot_selected, double*** adr_matrice_PSSM, double nb_sequence, char (*adr_Ct)[6], int taille_motif);
  // double calcul_distance_PSSMs(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv, double* distance_PSSM, int taille_motif);
  //
  // void calcul_score(TPtr_Mot_Ameliorer_PSSM* adr_mot, double*** adr_matrice_PSSM, int n_sequence, TPtr_dictionnaire_sequences* ptr_ensemble, int taille_motif);
  // int distanceHammingSt1(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st1, int taille_motif);
  // int distanceHammingSt2(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st2, int taille_motif);
  // int distanceHammingSt2_prim(char (*adr_Ct)[6], TPtr_dictionnaire_sequences* ptr_ensemble, TPtr_Mot_Ameliorer_PSSM *adr_mot, Ptr_st* adr_st2_prim, int taille_motif);
  //
  // void quick_sort_ST(Ptr_st* adr_st1, int* v_St1_Pos, int nb_sequences);
  // void trier(int* v_St1_Dh, int* v_St1_Pos, int g, int d);
  // void separer(int* v_St1_Dh, int* v_St1_Pos, int g, int d, int* adr_indice_pivot);
  //
  // void generation_fichier_résultats(Ptr_st* adr_st1, int* v_St1_Pos, char (*adr_Ct)[6], int nb_sequences);

#endif
