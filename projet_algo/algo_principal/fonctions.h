#ifndef HEADER_FONCTIONS
#define HEADER_FONCTIONS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <time.h>

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

 /*********************************
  *   DICTIONNAIRE DE SEQUENCES   *
  *********************************/
  typedef struct TDictionnaire_Sequences
  {
    int numero_sequence;		// kmer_present_dans_chaque_sequence(tete_info_dict_seq->nb_sequences, &tete_liste_kmer2,, &tete_liste_pour_recup_motif, &tete_liste_kmer_selectionne, &tete_liste_motif_PSSM);
    char nom_seq[120];
    char sequence[TAILLE_MAX_SEQ];
    struct TDictionnaire_Sequences* suiv_seq;
  } TDictionnaire_Sequences;
  typedef TDictionnaire_Sequences* TPtr_dictionnaire_sequences;

  typedef struct TInfo_dictionnaire_sequences
  {
    double nb_sequences;
    TPtr_dictionnaire_sequences tete_dictionnaire_seq; //tete pointant sur la liste chainée de sequences
  } TInfo_dictionnaire_sequences;
  typedef TInfo_dictionnaire_sequences* TPtr_info_dictionnaire_sequences;

/*********************************
 *     DICTIONNAIRE DE KMERs     *
 *********************************/
  //premier element du dictionnaire de kmer:: liste chainee des kmers
  typedef struct TCellkmer
  {
    struct TCellkmer* suiv_kmer; //Pointeur sur le kmer suivant
    struct TCellSequence* tete_sequence; //Pointer qui pointe sur une liste chain�e de s�quence o� le kmer est retrouv�
    double nb_sequence; // nombre de séquence dans lesquelles le kmer est présent
    int score_St1;
    int score_St2;
    double*** PSSM_consensus; //PSSM associée au motif consensus
    char* motif_consensus;
    char* kmer; //"TCC" par exemple
  } TCellkmer;
  typedef TCellkmer* TPtr_Cellkmer; //Pointeur sur Tcellkmer

  //deuxieme element du dictionnaire kmer: liste chainee de sequence possedant le kmer:
  typedef struct TCellSequence
  {
    int num_sequence; //correspond au numéro de la séquence dans le dictionnaire de séquences
    struct TCellSequence* suiv_sequence; // Pointeur qui pointe l'element suivant de la liste chainee de sequence
    struct TCellPos* tete_pos; // Pointeur qui pointe sur le premier element de la structure position qui repertorie toutes les positions o� le kmer a �t� trouv� dans une s�quence
    struct TCellPos* tete_pos_max;
    struct TCellPos* tete_pos_dH_min;
  } TCellSequence;
  typedef TCellSequence* TPtr_CellSequence; // Pointeur sur TCellSequence

  // troisieme element du dictionnaire de kmer: liste chainee de position ou le kmer a ete trouve dans une sequence.
  typedef struct TCellPos
  {
    int position; // position du kmer dans la s�quence
    double score;
    struct TCellPos* suiv_pos; //pointe sur la position suivante
    char* motif; //on stocke le motif équivalent au kmer
  } TCellPos;
  typedef TCellPos* TPtr_CellPos;

  // header contenant les tetes des elements précédents
  typedef struct TInfo_dictionnaire_kmer
  {
    int nb_kmer;
    TPtr_Cellkmer tete_liste_kmer;
  } TInfo_dictionnaire_kmer;
  typedef TInfo_dictionnaire_kmer* TPtr_info_dictionnaire_kmer;

/*************************************************************************************************************
 * * *                                           FONCTIONS                                               * * *
 *************************************************************************************************************/
  void importer_parametres(char* nom_fichier, int* taille_motif, int* d, int* nb_fenetres, int* nb_masques);
  void importer_sequences_fasta(char* nom_fichier_fasta, TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);
  void affichage_dictionnaire_sequences(TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);
  void liberation_dictionnaire_sequences( TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);

  int random_number(int max_number, int zero_excluded);
  void initialisation_masque(int** adr_masque);
  void affichage_masque(int* masque);
  void generation_masque(int** masque);
  void generation_dictionnaire_kmers(int position_kmer, char* k_mer, char* motif, TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void parcours_masque_sur_seq( int* masque, TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void affichage_dictionnaire_kmers(TPtr_info_dictionnaire_kmer tete_info_dict_kmer);
  void liberation_dictionnaire_kmers(TPtr_info_dictionnaire_kmer tete_info_dict_kmer);

  void creation_PSSM(double*** adr_matrice_PSSM);
  void afficher_PSSM( double** matrice_PSSM);
  double calcul_distance_PSSMs(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv );
  void copie_PSSM(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv);
  void liberation_PSSM(double*** adr_matrice_PSSM);
  void calcul_PSSM(TPtr_Cellkmer p_kmer , double*** adr_matrice_PSSM );
  void calcul_PSSM_amelioree(TPtr_Cellkmer p_kmer , double*** adr_matrice_PSSM );

  char* identification_motif_consensus(double** matrice_PSSM);
  double calculer_score(char* motif_tmp, double** matrice_PSSM,  TPtr_dictionnaire_sequences tete_dict_seq,  TPtr_info_dictionnaire_sequences tete_info_dict_seq);
  void recherche_motifs_maximisant_scores(TPtr_Cellkmer p_kmer, double** matrice_PSSM,  TPtr_dictionnaire_sequences tete_dict_seq, TPtr_info_dictionnaire_sequences tete_info_dict_seq);
  int distance_Hamming_St1(TPtr_Cellkmer p_kmer);
  int distance_Hamming_St2(TPtr_Cellkmer p_kmer);
  int distance_Hamming_St2_T_prim(TPtr_Cellkmer p_kmer);
  void recherche_motifs_minimisant_dHamming(TPtr_Cellkmer p_kmer,TPtr_info_dictionnaire_sequences tete_info_dict_seq);
  void raffiner_version1(TPtr_Cellkmer p_kmer);
  void egalisation_T(TPtr_Cellkmer p_kmer);
  int verification_convergence_T(TPtr_Cellkmer p_kmer);
  void raffiner_version2(TPtr_Cellkmer p_kmer, TPtr_info_dictionnaire_sequences tete_info_dict_seq);

  void trier_ST1(TPtr_Cellkmer* (*adr_tableau_kmer_QS), int g, int d);
  void separer_ST1(TPtr_Cellkmer* (*adr_tableau_kmer_QS), int g, int d, int* adr_indice_pivot);
  void trier_ST2(TPtr_Cellkmer* (*adr_tableau_kmer_QS), int g, int d);
  void separer_ST2(TPtr_Cellkmer* (*adr_tableau_kmer_QS), int g, int d, int* adr_indice_pivot);
  void quick_sort_ST(TPtr_Cellkmer* (*adr_tableau_kmer_QS), TPtr_info_dictionnaire_kmer tete_info_dict_kmer);

  void generation_fichier_resultats(TPtr_Cellkmer *tableau_kmer_QS, int numero_essais, TPtr_info_dictionnaire_kmer tete_info_dict_kmer, TPtr_info_dictionnaire_sequences tete_info_dict_seq);

#endif
