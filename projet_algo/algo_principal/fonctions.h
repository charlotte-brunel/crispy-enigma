#define TAILLE_MAX_SEQ 800

//---STRUCTURES------------------------------------------------

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


//---FONCTIONS--------------------------

void importer_parametres(int* l, int* d, int* k, int* nb_masques);
void importer_sequences_fasta(TInfo_ensemble_sequences* ptr_info, TEnsemble_Sequences* ptr_ensemble );
void afficher_sequences(TInfo_ensemble_sequences* ptr_info, TEnsemble_Sequences* ptr_ensemble  );
int random_number(int max_number, int zero_excluded);
void generation_masque(int l, int* masque[l], int k);


//---VARIABLES GLOBALES---------------------------------

extern bool convergence;
extern int l;
extern int d;
extern int k;
extern int nb_masques;
//
