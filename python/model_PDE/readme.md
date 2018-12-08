# MODEL v programskem jeziku Python

## Datoteke

- `parameters.py`: ima enako vlogo kot matlabova datoteka `params.m`, le da je **NI potrebno** posebej zaganjati
- `repressilator_S_PDE.py`: celoten model (simulacijska datoteka), ki uporablja parametre definirane v datoteki `parameters.py` in BO za simulacijo uporabljal model ODE definiran v `repressilator_S_ODE.py`
- `repressilator_S_ODE.py`: reakcijski del modela (brez difuzije), ki ga kliče glavna simulacijska datoteka (`repressilator_S_PDE.py`). 

## Poganjanje

v terminalu / konzoli poženi ukaz: `python repressilator_S_PDE.py`

