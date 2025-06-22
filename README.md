# TFG: Ajust d'equacions diferencials ordinàries a dades de sèries temporals curtes
## Introducció

En aquest projecte s'implementen en llenguatge de programació C, tres mètodes heurístics per a l'ajust de un model diferencial de la dinàmica poblacional de les gavines del Delta de l'Ebre. El model amb el que es treballa és el següent 


$\frac{d}{dt}x(t)= \varphi x(t) -\beta x(t)^2-\lambda D(x(t))$

$D(x,\mu,\sigma,\delta) :=
	\begin{cases} 
		\frac{1-\mathcal{E}_{dir}(x,\mu,\sigma,\delta)}{1-\mathcal{E}_{dir}(0,\mu,\sigma,\delta)}& \text{si } 0 \leq x \leq \delta \\
		\frac{1-\mathcal{E}(x,\sigma,\delta)}{1-\mathcal{E}_{dir}(0,\mu,\sigma,\delta)} & \text{si } x > \delta
	\end{cases}$

