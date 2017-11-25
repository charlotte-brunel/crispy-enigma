#define TAILLE_MAX_SEQ 800
/*******************************************************************************
 * * *                   VARIABLES GLOBALES                                * * *
 *******************************************************************************/

extern bool convergence;
extern int l;
extern int d;
extern int k;
extern int nb_masques;

/*******************************************************************************
 * * *                      STRUCTURES                                     * * *
 *******************************************************************************/
typedef struct TInfo_ensemble_sequences TInfo_ensemble_sequences;
struct TInfo_ensemble_sequences
{
  int nb_seq;
  struct TEnsemble_Sequences* tete_ensemble_seq;
};
typedef TInfo_ensemble_sequences* TPtr_info_ensemble_sequences;


typedef struct TEnsemble_Sequences TEnsemble_Sequences;
struct TEnsemble_Sequences
{
  char nom_seq[120];
  char seq[TAILLE_MAX_SEQ];
  struct TEnsemble_Sequences* suiv_seq;
};
typedef TEnsemble_Sequences* TPtr_ensemble_sequences;

/*
 *******************************
 * STRUCTURE DICTIONNAIRE KMER *
 *******************************
1er élément: liste chainée des kmers
*/
typedef struct TCellkmer TCellkmer;
struct TCellkmer{
    struct TCellkmer* suiv_kmer; //Pointer sur le kmer suivant
    struct TCellSequence* tete_sequence; //Pointer qui pointe sur une liste chainée de séquence où le kmer est retrouvé
    char kmer[]; //"TCC" par exemple

};
typedef TCellkmer* TPtr_Cellkmer; //Pointeur sur Tcellkmer

//deuxième élément du dictionnaire kmer: liste chainée de sequence possédant le kmer:

typedef struct TCellSequence TCellSequence;
struct TCellSequence{
    int sequence; //sequence 1 --> 1; sequence 2 --> 2
    struct TCellSequence* suiv_sequence; // Pointeur qui pointe l'élément suivant de la liste chainée de séquence
    struct TCellPos* tete_pos; // Pointeur qui pointe sur le premier élément de la structure position qui répertorie toutes les positions où le kmer a été trouvé dans une séquence
};
typedef TCellSequence* TPtr_CellSequence; // Pointeur sur TCellSéquence

// troisième élément du dictionnaire de kmer: liste chainée de position où le kmer a été trouvé dans une séquence.

typedef struct TCellPos TCellPos;
struct TCellPos{
    int position; // position du kmer dans la séquence
    struct TCellPos* suiv_pos; //pointe sur la position suivante
};
typedef TCellPos* TPtr_CellPos;

/*******************************************************************************
 * * *                      FONCTIONS                                      * * *
 *******************************************************************************/

void importer_parametres(int* l, int* d, int* k, int* nb_masques);
void importer_sequences_fasta(TInfo_ensemble_sequences* ptr_info, TEnsemble_Sequences* ptr_ensemble );
void afficher_sequences(TInfo_ensemble_sequences* ptr_info, TEnsemble_Sequences* ptr_ensemble  );
