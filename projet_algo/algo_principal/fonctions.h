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
  extern int nb_fenetre;
  extern int nb_masques;
/*************************************************************************************************************
 * * *                                         STRUCTURES                                                * * *
 *************************************************************************************************************/
  typedef struct TEnsemble_Sequences TEnsemble_Sequences;
  struct TEnsemble_Sequences
  {
   int numero_sequence;
   char nom_seq[120];
   char sequence[TAILLE_MAX_SEQ];
   struct TEnsemble_Sequences* suiv_seq;
  };
  typedef TEnsemble_Sequences* TPtr_ensemble_sequences;

  typedef struct TInfo_ensemble_sequences TInfo_ensemble_sequences;
  struct TInfo_ensemble_sequences
  {
    double nb_sequences;
    TPtr_ensemble_sequences tete_ensemble_seq; //tete pointant sur la liste chainée de sequences
  };
  typedef TInfo_ensemble_sequences* TPtr_info_ensemble_sequences;

/*******************************
 * STRUCTURE DICTIONNAIRE KMER *
 *******************************/
  //1er element: liste chainee des kmers
  typedef struct TCellkmer TCellkmer;
  struct TCellkmer{
    struct TCellkmer* suiv_kmer; //Pointeur sur le kmer suivant
    struct TCellSequence* tete_sequence; //Pointer qui pointe sur une liste chain�e de s�quence o� le kmer est retrouv�
    int nb_sequence; // nombre de séquence dans lesquelles le kmer est présent
    char kmer[]; //"TCC" par exemple
  };
  typedef TCellkmer* TPtr_Cellkmer; //Pointeur sur Tcellkmer

  //deuxieme element du dictionnaire kmer: liste chainee de sequence possedant le kmer:
  typedef struct TCellSequence TCellSequence;
  struct TCellSequence{
    int sequence; //sequence 1 --> 1; sequence 2 --> 2
    struct TCellSequence* suiv_sequence; // Pointeur qui pointe l'element suivant de la liste chainee de sequence
    struct TCellPos* tete_pos; // Pointeur qui pointe sur le premier element de la structure position qui repertorie toutes les positions o� le kmer a �t� trouv� dans une s�quence
  };
  typedef TCellSequence* TPtr_CellSequence; // Pointeur sur TCellSequence

  // troisieme element du dictionnaire de kmer: liste chainee de position ou le kmer a ete trouve dans une sequence.
  typedef struct TCellPos TCellPos;
  struct TCellPos{
    int position; // position du kmer dans la s�quence
    struct TCellPos* suiv_pos; //pointe sur la position suivante
  };
  typedef TCellPos* TPtr_CellPos;

/******************************
 * STRUCTURE KMER SELECTIONNE *
 ******************************/
  //1er element: liste chainee des kmers
  typedef struct TCellkmer_selectionne TCellkmer_selectionne;
  struct TCellkmer_selectionne{
    struct TCellkmer_selectionne* suiv_kmer_selectionne;
    struct TCell_Motif_PSSM* tete_motif_PSSM;
    int nb_sequence;
    char kmer[];
  };
  typedef TCellkmer_selectionne* TPtr_Cellkmer_selectionne;

  //2eme element: liste chainee de motif pour lequels on va calculer la PSSM:
  typedef struct TCell_Motif_PSSM TCell_Motif_PSSM;
  struct TCell_Motif_PSSM{
    struct TCell_Motif_PSSM* suiv_motif;
    char motif[];
  };
  typedef TCell_Motif_PSSM* TPtr_Cell_Motif_PSSM;

/****************************
 * LISTE CHAINEE DES MOTIFS *
 ****************************/
  typedef struct liste_chaine_motif liste_chaine_motif; //chaque motif sera stock� dans cette liste pour effectuer le calcul de la PSSM
  struct liste_chaine_motif{
    char motif_substitue[6];
    struct liste_chaine_motif *next_motif; //pointeur sur le prochain num�ro de s�quence (ex: '>Seq2:')
  };
  typedef liste_chaine_motif* ptr_liste_motif;

/****************************************************
 * STRUCTURE CHAINEE DE MOTS POUR AMELIORER LA PSSM *
 ****************************************************/
  typedef struct TMot_Ameliorer_PSSM TMot_Ameliorer_PSSM;
  struct TMot_Ameliorer_PSSM
  {
    double score_mot; //score du mot
    struct TMot_Ameliorer_PSSM * next_mot; //pointeur sur mot suivant
    char mot[]; //Mot de longueur motif
  };
  typedef TMot_Ameliorer_PSSM* TPtr_Mot_Ameliorer_PSSM;

/************************************************
 * STRUCTURE CHAINEE ST1/ST2 CONTENANT LES MOTS *
 ************************************************/
  typedef struct st st;
  struct st
  {
   int distance_hamming; //distance de hamming entre le mot et Ct
   struct st * next_mot; //pointeur sur mot suivant
   char mot[]; //Mot qui ont une distance de hamming supérieur à 2
  };
  typedef st* Ptr_st;

/*************************************************************************************************************
 * * *                                           FONCTIONS                                               * * *
 *************************************************************************************************************/
  void importer_parametres(int* taille_motif, int* d, int* nb_fenetre, int* nb_masques);
  void importer_sequences_fasta(TPtr_info_ensemble_sequences* ptr_info, TPtr_ensemble_sequences* ptr_ensemble );
  void afficher_sequences(TPtr_info_ensemble_sequences* ptr_info, TPtr_ensemble_sequences* ptr_ensemble );

  int random_number(int max_number, int zero_excluded);
  void generation_masque(int longueur_masque, void* adr_masque, int nb_fenetre);

  void generation_kmer(int position_kmer, char* k_mer, TPtr_Cellkmer* adr_liste_kmer, TPtr_CellSequence* adr_liste_sequence, TPtr_ensemble_sequences* ptr_ensemble, TPtr_CellPos* adr_liste_pos);
  void parcours_masque(int longueur_masque, void* adr_masque, int nb_fenetre, int nb_sequence, TPtr_ensemble_sequences* ptr_ensemble, TPtr_Cellkmer* adr_tete_liste_kmer, TPtr_CellSequence* adr_tete_liste_sequence, TPtr_CellPos* adr_tete_liste_pos);

  void recuperer_motif_kmer(TPtr_Cellkmer* adr_parcours_kmer, TPtr_Cellkmer_selectionne *adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif_PSSM, TPtr_CellSequence* adr_cell_sequence, TPtr_CellPos* adr_cell_pos, TPtr_ensemble_sequences* ptr_ensemble, int nb_sequence_kmer);
  void kmer_present_dans_chaque_sequence(int nb_sequence, TPtr_Cellkmer* adr_cell_kmer, TPtr_CellSequence *adr_cell_sequence, TPtr_CellPos *adr_cell_pos, TPtr_ensemble_sequences* ptr_ensemble, TPtr_Cellkmer_selectionne* adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_cell_motif_PSSM);

  void affichage_dictionnaire_kmer(TPtr_Cellkmer* adr_tete_kmer, TPtr_CellSequence* adr_tete_sequence, TPtr_CellPos* adr_tete_pos);
  void affichage_motif_selectionne(TPtr_Cellkmer_selectionne* adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif);

  void calcul_PSSM(TPtr_Cellkmer_selectionne *adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM *adr_cell_motif_PSSM, double*** adr_matrice_PSSM, int taille_motif);
  void calcul_nouvelle_PSSM(TPtr_Cell_Motif_PSSM *adr_cell_mot_selected, double*** adr_matrice_PSSM, double nb_sequence, char (*adr_Ct)[6], int taille_motif);
  void afficher_PSSM( double*** adr_matrice_PSSM, int taille_motif);

  void calcul_score(TPtr_Mot_Ameliorer_PSSM* adr_mot, double*** adr_matrice_PSSM, int n_sequence, TPtr_ensemble_sequences* ptr_ensemble, int taille_motif);
  double dist_PSSM(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv, double* distance_PSSM, int taille_motif);

  int distanceHammingSt1(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st1, int taille_motif);
  int distanceHammingSt2(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st2, int taille_motif);
  int distanceHammingSt2_prim(char (*adr_Ct)[6], TPtr_ensemble_sequences* ptr_ensemble, TPtr_Mot_Ameliorer_PSSM *adr_mot, Ptr_st* adr_st2_prim, int taille_motif);

  void quick_sort_ST(Ptr_st* adr_st1, int* v_St1_Pos, int nb_sequences);
  void trier(int* v_St1_Dh, int* v_St1_Pos, int g, int d);
  void separer(int* v_St1_Dh, int* v_St1_Pos, int g, int d, int* adr_indice_pivot);

  void fichier_sortie_st(Ptr_st* adr_st1, int* v_St1_Pos, char (*adr_Ct)[6], int nb_sequences);

#endif
