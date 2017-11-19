//genere un chiffre random

int random_number(int max_number, int zero_excluded){
	int randomNumber;
	if(zero_excluded==0){ //on peut tomber sur 0 aléatoirement
		randomNumber= rand() % max_number;
	}else{
        randomNumber= rand() % max_number +1;
	}
	return(randomNumber);
	}
  
  //genere un masque avec un nombre de fenetre ouverte en paramètre:
  
  void generation_masque(int longueur_masque, int* adr_masque[longueur_masque], int nb_fenetre){
    int i;
    int nb_fenetre_ouverte=1;
    while (nb_fenetre_ouverte < nb_fenetre){
        for (i=0; i<= (longueur_masque -1); i++){
            adr_masque[i]= random_number(2,0);
            if (adr_masque[i]==1){
                nb_fenetre_ouverte++;
            }
        }
    }
    return;

}


/*  APPEL DE LA FONCTION GENERATION_MASQUE

int longueur_masque=5;
int adr_masque[longueur_masque];
int nb_fenetre=2;

generation_masque(longueur_masque, &adr_masque[longueur_masque], nb_fenetre);

*/
