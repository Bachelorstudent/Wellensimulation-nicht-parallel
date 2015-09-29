# Wellensimulation-nicht-parallel

Dies ist die Wellensimulation ohne Parallelisierung.
Es wurden keine experimentellen oder halbfertige Funktionen hochgeladen.
Dieser Code wurde in der Bachelorarbeit zum ermitteln der Simulationsgeschwindigkeiten verwendet.

Kompilierung:
```
gcc -Wall -Wextra -Werror -O3 -march=native -funroll-loops -c -o wavestate3d.o -lm wavestate3d.c
gcc -Wall -Wextra -Werror -O3 -march=native -funroll-loops -c -o wave3d.o -lm wave3d.c
gcc wave3d.o wavestate3d.o -o wave3d -lm
````
Ausführung:
```
./wave3d nx ny nz iterations
```
