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
    int numero_sequence; //numero de la séquence (1ère, 2ème ... N ème)
    char nom_seq[120]; // nom de la séquence dans le fichier FASTA ( "> Seq 1" )
    char sequence[TAILLE_MAX_SEQ]; // Séquence ("ataccgacagt")
    struct TDictionnaire_Sequences* suiv_seq; // Pointeur sur la prochaine séquence
  } TDictionnaire_Sequences;
  typedef TDictionnaire_Sequences* TPtr_dictionnaire_sequences; // Pointeur sur la structure TDictionnaire_Sequences

  typedef struct TInfo_dictionnaire_sequences
  {
    int nb_sequences; // nombre de séquence dans le fichier FASTA
    TPtr_dictionnaire_sequences tete_dictionnaire_seq; // tête pointant sur la liste chainée de sequences
  } TInfo_dictionnaire_sequences;
  typedef TInfo_dictionnaire_sequences* TPtr_info_dictionnaire_sequences; // Pointeur sur la structure TInfo_dictionnaire_sequences

/*******************************
 * STRUCTURE DICTIONNAIRE KMER *
 *******************************/
 
  //1er élément: liste chainée des kmers
  typedef struct TCellkmer
  {
    struct TCellkmer* suiv_kmer; //Pointeur sur le kmer suivant
    struct TCellSequence* tete_sequence; //Pointeur sur une liste chainée de séquences où le kmer est trouvé
    double nb_sequence; // nombre de séquences dans lesquelles le kmer est présent
    char kmer[]; //"TCC" par exemple
  } TCellkmer;
  typedef TCellkmer* TPtr_Cellkmer; //Pointeur sur Tcellkmer

  //deuxieme element du dictionnaire kmer: liste chainee de sequence possedant le kmer:
  typedef struct TCellSequence
  {
    int sequence; // Numéro de la séquence (sequence 1 --> 1; sequence 2 --> 2)
    struct TCellSequence* suiv_sequence; // Pointeur sur l'élément suivant de la liste chainée de sequence
    struct TCellPos* tete_pos; // Pointeur sur le premier élément de la structure position 
  } TCellSequence;
  typedef TCellSequence* TPtr_CellSequence; // Pointeur sur TCellSequence

  // troisième élément du dictionnaire de kmer: liste chainée de position où le kmer a été trouvé dans une sequence.
  typedef struct TCellPos
  {
    int position; // position du kmer dans la séquence
    struct TCellPos* suiv_pos; //pointe sur la position suivante
  } TCellPos;
  typedef TCellPos* TPtr_CellPos; // Pointeur sur TCellPos

/******************************
 * STRUCTURE KMER SELECTIONNE * : Permet de stocker les Kmers présents dans chaque séquence. 
 ******************************   Ils seront ensuite utilisé pour calculer la PSSM */
 
  //1er element: liste chainee des kmers
  typedef struct TCellkmer_selectionne
  {
    struct TCellkmer_selectionne* suiv_kmer_selectionne; // Pointeur sur le kmer suivant 
    struct TCell_Motif_PSSM* tete_motif_PSSM; // Pointeur sur la structure TCell_Motif_PSSM
    double nb_sequence; // Nombre de séquence où l'on trouve le kmers (défini en double pour effectuer des opérations dessus)
    char kmer[]; // Kmer
  } TCellkmer_selectionne;
  typedef TCellkmer_selectionne* TPtr_Cellkmer_selectionne; // Pointeur sur TCellkmer_selectionne

  //2eme element: liste chainee de motif pour lequels on va calculer la PSSM:
  typedef struct TCell_Motif_PSSM
  {
    struct TCell_Motif_PSSM* suiv_motif; // Pointeur sur le motif suivant
    char motif[]; // Motif (par exemple "aatca")
  } TCell_Motif_PSSM;
  typedef TCell_Motif_PSSM* TPtr_Cell_Motif_PSSM; // Pointeur que la structure TCell_Motif_PSSM

/****************************************************
 * STRUCTURE CHAINEE DE MOTS POUR AMELIORER LA PSSM * : Permet la sélection des mots maximisant le score
 ****************************************************/
  typedef struct TMot_Ameliorer_PSSM
  {
    double score_mot; //score du mot
    struct TMot_Ameliorer_PSSM * next_mot; //pointeur sur mot suivant
    char mot[]; //Mot de longueur motif
  } TMot_Ameliorer_PSSM;
  typedef TMot_Ameliorer_PSSM* TPtr_Mot_Ameliorer_PSSM; // Pointeur sur TMot_Ameliorer_PSSM

/************************************************
 * STRUCTURE CHAINEE ST1/ST2 CONTENANT LES MOTS * : Contient les score ST1 et ST2
 ************************************************/
  typedef struct st
  {
    // position /s si meme motif pour plusieurs séquences
    int distance_hamming; //distance de hamming entre le mot et Ct
    struct st * next_mot; //pointeur sur mot suivant
    char mot[]; //Mot qui ont une distance de hamming supérieur à 2
  } st;
  typedef st* Ptr_st; // Pointeur sur st

/*************************************************************************************************************
 * * *                                           FONCTIONS                                               * * *
 *************************************************************************************************************/
  void importer_parametres(int* taille_motif, int* d, int* nb_fenetres, int* nb_masques);
  void importer_sequences_fasta(TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);
  void afficher_sequences(TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);
  void liberation_dictionnaire_sequence( TPtr_info_dictionnaire_sequences* adr_tete_info_dict_seq, TPtr_dictionnaire_sequences* adr_tete_dict_seq);


/**********
RANDOM_NUMBER est une fonction qui sélectionne un nombre au hasard. Le paramètre zero_excluded permet d'indiquer si 0 fait parti du choix 
pour la sélection du nombre. max_number indique le seuil maximal. Par exemple random_number(5, 0) sélectionnera au hasard un chiffre entre 0 et 5.
random_number (12, 1) sélectionnera au hasard un nombre entre 1 et 12. Cette fonction est utilisée dans generation_masque.
**********/

  int random_number(int max_number, int zero_excluded);
  
/**********
GENERATION_MASQUE permet de génerer le masque qui va parcourir les séquences. Le masque est représenté par un vecteur de 0 et de 1:
__________
|0|1|0|1|0|  : Un zéro signifie que la fenêtre est fermée, le 1, que la fenêtre est ouverte.

On ne peut pas avoir un masque avec un nombre de fenêtre ouverte/fermée différent que celui spécifié par l'utilisateur
**********/
  void generation_masque(void* adr_masque);

/**********
PARCOURS_MASQUE: Dans cette fonction le masque parcourt les séquences et récupère les kmers, leurs positions, leurs séquences pour remplir 
le dictionnaire de kmers. Pour chaque séquence, le masque se déplace sur la séquence. On récupère le premier kmer, on se déplace avec le masque 
sur la séquence d'une longueur de masque, on récupère le deuxième kmer et ainsi de suite. Les informations sur les kmers sont envoyés à la 
fonction generation_kmer pour remplir le dictionnaire
**********/

  void parcours_masque( void* adr_masque, TPtr_dictionnaire_sequences* ptr_ensemble, TPtr_Cellkmer* adr_tete_liste_kmer, TPtr_CellSequence* adr_tete_liste_sequence, TPtr_CellPos* adr_tete_liste_pos);
/**********
GENERATION_KMER: cette fonction permet de remplir le dictionnaire de kmer. On utilise pour chaque structure deux pointeurs pour le chainage.
Un pointeur p qui pointe sur la brique actuelle et un p_prec qui pointe sur la brique précédente. On distingue 3 cas possibles:

1) Le dictionnaire est vide: On crée une brique kmer, une brique sequence et une brique pos. La brique kmer pointe sur la brique sequence qui pointe 
sur la brique pos. On remplit les trois briques avec les bonnes informations. Le pointeur p sur chaque brique suivante pointe sur NULL. p_prec=p
car on a qu'une seule brique par structure.
2) On ajoute un kmer qui a déjà été trouvé: on vérifie si ce kmer a été trouvé dans une séquence connue ou une nouvelle séquence. Si la séquence est
connue on créera juste une nouvelle brique position. Sinon il faudra créer une nouvelle brique séquence qui pointe sur une nouvelle brique position.
Pour le chainage on a p qui pointe sur la brique nouvellement créée, p->suiv = NULL, et p_prec pointe sur la brique précedemment créée. On fait 
pointer p_prec -> suiv = p et on déplace p_prec sur p.
3) On ajoute un nouveau kmer: en parcourant la liste de kmer on n'en a pas trouvé de similaire. On crée une brique kmer, une brique séquence et une
brique position. De nouveau, la brique kmer pointe sur la brique sequence qui pointe sur la brique pos. On remplit les trois briques avec les bonnes
informations. On réalise le chainage précedemment décrit avec p et p_prec.
**********/
  
  void generation_kmer(int position_kmer, char* k_mer, TPtr_Cellkmer* adr_liste_kmer, TPtr_CellSequence* adr_liste_sequence, TPtr_dictionnaire_sequences* ptr_ensemble, TPtr_CellPos* adr_liste_pos);

/**********
RECUPERER_MOTIF_KMER et KMER_PRESENT_DANS_CHAQUE_SEQUENCE: ces deux fonctions marchent de paires et permettent de récupérer les occurrences 
correspondant aux kmers présents dans chaque séquence. Dans un premier temps, la fonction k_mer_present_dans_chaque_sequence parcourt le dictionnaire
de kmer et compte le nombre de séquences où il est présent. Si le kmer est présent dans 80% des séquences on le sélectionne et on le sauvegarde dans
la structure TCellkmer_selectionne. La sauvegarde se déroule dans la fonction recuperer_motif_kmer: au même moment on parcourt les séquences pour
récupérer le motif correspondant au kmer (il sera stocké dans la liste chaînée TCell_Motif_PSSM).
**********/
  //
  // void recuperer_motif_kmer(TPtr_Cellkmer* adr_parcours_kmer, TPtr_Cellkmer_selectionne *adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif_PSSM, TPtr_CellSequence* adr_cell_sequence, TPtr_CellPos* adr_cell_pos, TPtr_dictionnaire_sequences* ptr_ensemble, int nb_sequence_kmer);
  // void kmer_present_dans_chaque_sequence(int nb_sequence, TPtr_Cellkmer* adr_cell_kmer, TPtr_CellSequence *adr_cell_sequence, TPtr_CellPos *adr_cell_pos, TPtr_dictionnaire_sequences* ptr_ensemble, TPtr_Cellkmer_selectionne* adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_cell_motif_PSSM);
  //
  
/**********
AFFICHAGE_DICTIONNAIRE_KMER et AFFICHE_MOTIF_SELECTIONNE: affiche le contenu du dictionnaire de kmer initial dans 'dictionnaire_kmer.txt'
et le contenu du dictionnaire des kmers sélectionnés (car présents dans la majorité des séquences) dans 'kmer_selectionne.txt'.  Ce sont
des fichiers sorties intermédiaires qui permettent de mieux visualiser le déroulement du programme.
**********/
  // void affichage_dictionnaire_kmer(TPtr_Cellkmer* adr_tete_kmer, TPtr_CellSequence* adr_tete_sequence, TPtr_CellPos* adr_tete_pos);
  // void affichage_motif_selectionne(TPtr_Cellkmer_selectionne* adr_tete_kmer_selectionne, TPtr_Cell_Motif_PSSM* adr_tete_motif);
  //
  void creation_PSSM(double*** adr_matrice_PSSM); 
  
/**********
CALCUL_PSSM et CALCUL_NOUVELLE_PSSM: la fonction calcul PSSM permet de calculer la PSSM et de déterminer le motif consensus à partir des kmers 
sélectionnés précédemment. On compare les motifs récupérés pour chaque occurrence du kmer et on calcule le pourcentage de chaque nucléotide à 
travers les motifs à chaque position. Le motif consensus est constitué par les nucléotides les plus présentes parmi les motifs à chaque position.
La fonction calcul_nouvelle_PSSM réalise les mêmes étapes mais ne possède pas les même paramètres d'entrées. Une possible amélioration du 
code serait de fusionner ces deux fonctions.
**********/

  // void calcul_PSSM(TPtr_Cellkmer_selectionne *adr_cell_kmer_selectionne, TPtr_Cell_Motif_PSSM *adr_cell_motif_PSSM, double*** adr_matrice_PSSM, int taille_motif);
  // void calcul_nouvelle_PSSM(TPtr_Cell_Motif_PSSM *adr_cell_mot_selected, double*** adr_matrice_PSSM, double nb_sequence, char (*adr_Ct)[6], int taille_motif);
  void afficher_PSSM( double** adr_matrice_PSSM);
  void liberation_PSSM(double*** adr_matrice_PSSM);
  //
  
/**********
CALCUL_SCORE: pour chaque mot de chaque séquence on calcule le score de ce mot. score= p(mot|PSSM)/p(mot|Background).
**********/

  // void calcul_score(TPtr_Mot_Ameliorer_PSSM* adr_mot, double*** adr_matrice_PSSM, int n_sequence, TPtr_dictionnaire_sequences* ptr_ensemble, int taille_motif);
/**********
DIST_PSSM: calcule la distance entre deux matrices PSSM. Cette distance est calculée de la façon suivante: pour chaque case des deux matrices on fait:
| (case[i][j]_matrice1) - (case[i][j]_matrice2) | où i représente le nombre de colonne et j le nombre de ligne. On somme le score de chaque case pour 
obtenir la distance totale entres les deux matrices.
**********/
  
  // double dist_PSSM(double*** adr_matrice_PSSM, double*** adr_matrice_PSSM_nouv, double* distance_PSSM, int taille_motif);
  //
/**********
DISTANCEHAMMINGST1, DISTANCEHAMMINGST2 et DISTANCEHAMMINGST2_PRIM: dans ces trois fonctions on calcule la distance de Hamming entre chaque mot et 
le motif consensus. Contrairement à distanceHammingSt1 et distanceHammingSt2, dans distanceHammingSt2_prim on ne garde que le mot de la séquence
qui minimise la distance de Hamming avec le motif consensus. Dans les autres fonctions la distance de Hamming était gardée et stockée pour tous les 
mots.
**********/
  // int distanceHammingSt1(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st1, int taille_motif);
  // int distanceHammingSt2(char (*adr_Ct)[6], TPtr_Cell_Motif_PSSM* adr_mot_selected, Ptr_st* adr_st2, int taille_motif);
  // int distanceHammingSt2_prim(char (*adr_Ct)[6], TPtr_dictionnaire_sequences* ptr_ensemble, TPtr_Mot_Ameliorer_PSSM *adr_mot, Ptr_st* adr_st2_prim, int taille_motif);
  //
/**********
QUICK_SORT, TRIER et SEPARER: les mots, leurs scores ST1/ST2 et leurs positions sont stockés dans les listes chaînées st. Pour trier ces listes chaînées
on procéde de la façon suivante: dans quick_sort deux vecteurs sont créés à partir des listes chaînées
-le premier vecteur contient le score de chaque mot (on a autant de scores et de mots que de séquences) 
-Le deuxième vecteut contient la position des scores dans la liste chaînée.
On trie le premier vecteur avec la méthode quick_sort mais à chaque déplacement on réalise le même dans le vecteur de position.
En sortie on aura un vecteur de score trié et le vecteur de position indiquera où chercher les scores triés dans la liste chainée.
**********/
  // void quick_sort_ST(Ptr_st* adr_st1, int* v_St1_Pos, int nb_sequences);
  // void trier(int* v_St1_Dh, int* v_St1_Pos, int g, int d);
  // void separer(int* v_St1_Dh, int* v_St1_Pos, int g, int d, int* adr_indice_pivot);
  //
/**********
FICHIER_SORTIE_ST: cette fonction permet d'afficher les motifs sélectionnés dans chaque séquence avec leurs positions et leurs scores ST1 et ST2.
Pour afficher les motifs dans le bon ordre de score on utilise le vecteur de position générer avec le quick-sort. Dans ce vecteur, on aura dans la
première case, la position du mot avec le score le plus petit. Pour récupérer ses informations il faudra alors parcourir la liste jusqu'à la bonne 
position.
**********/
  // void fichier_sortie_st(Ptr_st* adr_st1, int* v_St1_Pos, char (*adr_Ct)[6], int nb_sequences);
  

#endif
