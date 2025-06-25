set terminal pngcairo size 800,600
unset key
set output '../Figures/Ajust Model/Ajust Model.png'

plot "../Resultats/Solucio EDO/Solucio EDO.txt" with lines linecolor rgb "blue" lw 3, \
     "../Resultats/Dades/Dades.txt" with points pt 7 lc rgb "red"

unset output

