# Gradientes latitudinales en la estructura y estabilidad de redes tróficas marinas del Atlántico Sudoccidental


## Methods 


### Metrics

### Model description

We consider that the compiled food webs are a representation of all posible interactions, and at any given time if a species is not present in the local site the interaction do not realize, so the number of species and interactions can fluctuate over time. To represent this network variability we used a metaweb assembly model [@Saravia2022], assuming that each compiled food web represent the metaweb for each region.

In this model species migrate from the metaweb to a local site with a probability $c$, and become extinct from the local site with probability $e$; a reminiscence of the theory of island biogeography (MacArthur & Wilson, 1967), but with the addition of network structure. Species migrate with their potential network links from the metaweb, then in the local site only the interactions with the present species are realized. Species have a probability of secondary extinction $se$ if none of its preys are present, which only applies to non-basal species. When a species goes extinct locally it may produce secondary extinctions modulated by $se$ .

The model is implemented in R using the package `meweasmo` [@Saravia2022b]. The 3 parameters of the model are fitted for each compiled food web using approximate Bayesian computation (ABC) [@Beaumont2002]. The ABC algorithm is implemented in the package `abc` [@abc] using the `abc` function with the `rejection` method. The distance function used is the euclidean distance between the observed and simulated number of species and connectance. The prior distributions for the parameters are uniform: $c \sim U(0.01,0.3)$, $e \sim U(0.001,0.5)$, and $se \sim U(0.1,0.9)$. The number of simulations is set to 10000, and the tolerance is set to 0.05.

