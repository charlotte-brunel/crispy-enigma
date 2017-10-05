typedef Struct{ //structure des colonnes
  int col;
  double val;
  struct TCellCol* suiv;
} TCellCol;
typedef TCellCol* TPtrCellCol;

typedef Struct{ //structure des lignes
  int lig;
  struct TCellLig* suiv;
  TPtrCellCol tete_cols;
} TCellLig;
typedef TCellLig* TPtrCellLig;

typedef Struct{ //structure des lignes
  int nb_max_ligs;
  int nb_max_cols;
  double val_def;
  TPtrCellLig tete_ligs;
} TCellMat;
typedef TCellMat* TMat;

// définition de l'initialisation de la fonction
void init(TMat p_mat, int nb_max_lignes, int nb_max_colonnes, double val_def, TPtrCellLig tete_ligs){
  p_mat->nb_max_ligs = nb_max_lignes;
  p_mat->nb_max_cols = nb_max_colonnes;
  p_mat->val_def = val_def;
  p_mat->tete_ligs = tete_ligs;
}

//appel de l'initialisation
init(p_mat,10,3,4,null);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void recherche_i(TPtrCellLig tete_ligs, int i, TPtrCellLig* adr_p_prec_i, TPtrCellLig* adr_p_i){
  TPtrCellLig p, prec;
  P = tete_ligs;
  prec = null;
  while(p != null) {
    if(p->lig == i) {
      *adr_p_i = p; *adr_p_prec_i= prec;
      return;
    }
    if(p->lig > i){
      *adr_p_i = null;; adr_p_prec_i = prec;
      return;
    }
  }
  *adr_p_i = null;
  if(tete_ligs == null) { *adr_p_prec_i = null; }
    else { *adr_p_prec_i = prec; }
}


void recherche_j(TPtrCellCol tete_Cols, int j, TPtrCellCol* adr_p_prec_j, TPtrCellCol* adr_p_j){
  TPtrCellCol p, prec;
  P = tete_Cols;
  prec = null;
  while(p != null) {
    if(p->Col == j) {
      *adr_p_j = p; *adr_p_prec_j= prec;
      return;
    }
    if(p->Col > j) {
      *adr_p_j = null;; adr_p_prec_j = prec;
      return;
    }
  }
  *adr_p_j = null;
  if(tete_Cols == null) { *adr_p_prec_j = null; }
    else { *adr_p_prec_j = prec; }
}


bool recherche_i_j(TMat p_mat, int i, int j, TPtrCellLig* adr_p_prec_i, TPtrCellLig* adr_p_i, TPtrCellCol* adr_p_prec_j, TPtrCellCol* adr_p_j) {

  recherche_i(p_mat, i, adr_p_prec_i, adr_p_i);
  if(*adr_p_i == null){ // la cellule i n'existe pas
    return(false);
  }
  recherche_j((*adr_p_i)->tete_cols, j, adr_p_prec_j, adr_p_j);
  if (*adr_p_j == null){
    retun(false);
  }
  return(true);
}


void suppression_ij(TMat p_mat, TPtrCellLig p_prec_i, TPtrCellLig* adr_p_i,TPtrCellCol p_prec_j, TPtrCellCol* adr_p_j){
  if(p_prec_j == null) { *adr_p_i->tete_cols = *adr_p_j->suiv; }
  else{ p_prec_j-suiv = *adr_p_j->suiv; }
  free( *adr_p_j);
  *adr_p_j = null;
  if(*adr_p_i->tete_cols == null) {
    if (p_prec_i == null){ p_mat->tete_ligs = *adr_p_i->suiv;}
    else{ p_prec_i->suiv = *adr_p_i->suiv;}
    free( *adr_p_i);
    *adr_p_i = null;
  }
}


void insertion(TMat p_mat, TPtrCellLig p_prec_i, TPtrCellLig* adr_p_i, TPtrCellCol p_prec_j, TPtrCellCol* adr_p_j, int i, int j, double val){
  if(*adr_p_i == null) {
    *adr_p_i =(TPtrCellLig)malloc(sizeof TCellLig);
    if(*adr_p_i == null) { printf("Pas assez de mémoire"); exit(1);}
    *adr_p_i->lig = i;
    if(p_prec_i == null) {
      *adr_p_i->suiv = p_mat->tete_ligs;
      p_mat->tete_ligs = *adr_p_i;
    }else{
      *adr_p_i->suiv = p_prec_i->suiv;
      p_prec_i->suiv = *adr_p_i;
    }
  }

  *adr_p_j =(TPtrCellCol)malloc(sizeof TCellCol);
  if(*adr_p_j == null) { printf("Pas assez de mémoire"); exit(1);}
  *adr_p_j->col =j;
  if(p_prec_j == null) {
    *adr_p_j->suiv = *adr_p_i->tete_cols;
    *adr_p_i->tete_cols = *adr_p_j;
  }else{
    *adr_p_j->suiv = p_prec_j->suiv;
    p_prec_j->suiv = *adr_p_j;
  }
}

// définition de la fonction Set()
void Set(TMat p_mat, int i, int j, double val){
  TPtrCellLig p_i; p_prec_i;
  TPtrCellCol p_j; p_prec_j;
  bool exist_struct;
  exist_struct = recherche_i_j(p_mat, i, j, &p_prec_i, &p_i, &p_prec_j, &p_j);
  switch(exist_stuct){
    true:switch(val == p_mat->val_def){
            true: { suppression_ij(p_mat, p_prec_i, &p_i, p_prec_j, &p_j); } break;
            false: { p_j->val = val; } break;
          }
    break;
    false:switch(val == p_mat->val_def){
            true: {} break;
            false: { insertion(p_mat, &p_prec_i, &p_i, &p_prec_j, &p_j, i, j, val);} break;
          }
    break;
  }
}

// définition de la fonction Get()
double Get(TMat p_mat, int i, int j){
  TPtrCellLig p_i; p_prec_i;
  TPtrCellCol p_j; p_prec_j;
  bool exist_struct;
  exist_struct = recherche_i_j(p_mat, i, j, &p_prec_i, &p_i, &p_prec_j, &p_j);
  if(exist_struct) { return(p_j->val); }
  return(p_mat->val_def);
}
