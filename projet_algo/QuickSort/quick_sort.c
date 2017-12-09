#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "quick_sort.h"

//--------------------------------------------------------------------------------------------------
void trier(int* v, int g, int d)
{
  int indice_pivot;

  if (g < d)
  {
    separer(v, g, d, &indice_pivot);
    indice_pivot --;
    trier(v, g, indice_pivot-1);
    trier(v, indice_pivot+1, d);
  }
  return;
}
//--------------------------------------------------------------------------------------------------
void separer(int* v, int g, int d, int* adr_indice_pivot)
{
  //PRECONDITIONS
  // g < d

  int bas, haut; //indices de position dans le vecteur
  int comp, pivot;  //comp pour comparateur, c'est variables sont des entiers car on trie un vecteur d'entiers

  bas = g;
  haut = d;
  pivot = v[bas];
  comp = v[haut];

  while (bas < haut)
  {
    if ( comp > pivot)
    {
      v[haut] = comp;
      haut --;
      comp = v[haut];
    }else{
      v[bas] = comp;
      bas ++;
      comp = v[bas];
    }
  }
  v[bas] = pivot;
  *adr_indice_pivot = bas;

  //POST-CONDITIONS
  return;
}

/*************************************************
***                 MAIN                       ***
**************************************************/
int main()
{
  int i;
  int borne_gauche = 0;
  int borne_droite = 4;
  int v[] = {2,6,4,1,9};

  for(i=0; i<5; i++)
  {
    printf("%d", v[i]);
  }
  printf("\n");

  trier(v, borne_gauche, borne_droite);

  for(i=0; i<5; i++)
  {
    printf("%d", v[i]);
  }
  printf("\n");

  return(0);
}
