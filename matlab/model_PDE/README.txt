Datoteke:
- params.m: nastavitev vrednosti parametrov in njihovo shranjevanje v datoteko params.mat, ki jo preberemo v simulacijski datoteki (po�eni pred prvo simulacijo);
- repressilator_S_PDE.m: celoten model (simulacijska datoteka), ki uporablja parametre definirane v datoteki params.mat in za simulacijo uporablja model ODE definiran v repressilator_S_ODE.m;
- repressilator_S_ODE.m: reakcijski del modela (brez difuzije), ki ga kli�e glavna simulacijsa datoteka (repressilator_S_PDE.m). 
 