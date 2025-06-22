# TFG: Ajust d'equacions diferencials ordinàries a dades de sèries temporals curtes
## Introducció

En aquest projecte s'implementen en llenguatge de programació C, tres mètodes heurístics per a l'ajust de un model diferencial de la dinàmica poblacional de les gavines del Delta de l'Ebre. El model amb el que es treballa és el següent 


$\frac{d}{dt}x(t)= \varphi x(t) -\beta x(t)^2-\lambda D(x(t))$

D(x, μ, σ, δ) = {
    (1 - ℰ_dir(x, μ, σ, δ)) / (1 - ℰ_dir(0, μ, σ, δ))     si 0 ≤ x ≤ δ
    (1 - ℰ(x, σ, δ))       / (1 - ℰ_dir(0, μ, σ, δ))     si x > δ



