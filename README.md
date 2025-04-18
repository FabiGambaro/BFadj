# Navigating sampling bias in discrete phylogeographic analysis: assessing the performance of an adjusted Bayes factor

Bayesian phylogeographic inference is a popular tool in molecular epidemiological studies, enabling the reconstruction of the dispersal history of rapidly evolving pathogens. Discrete phylogeographic analysis integrates geographic information as discrete characters and infers lineage transition events among discrete locations. A discrete phylogeographic analysis is typically followed by a Bayes factor (BF) test that assesses the statistical support for inferred transition links by comparing their posterior and prior expectations. In the standard BF (BF<sub>std</sub>) test approach, the relative abundance of the involved trait states is not considered, which can be problematic in the case of unbalanced sampling among discrete locations. Although several strategies have been proposed to address sampling bias in discrete phylogeographic analyses that employ continuous-time Markov chain (CTMC) modeling, they might require additional epidemiological information. In this study, we formally assess the performance of a modification of the BF<sub>std</sub>, the adjusted Bayes factor (BF<sub>adj</sub>), which attempts to incorporate information on the relative abundance of samples by location when inferring support for transition events and root location inference without requiring external data. Using a simulation framework, we assess the performance of BF<sub>std</sub> and BF<sub>adj</sub> under varying levels of sampling bias, estimating their type I (false positive) and type II (false negative) error rates. Our results show that BF<sub>adj</sub> complements the BF<sub>std</sub> (i) by reducing type I errors, leading to fewer false positive inferred transition events, but at the cost of increasing type II errors and thus decreasing statistical power, particularly in highly biased scenarios; and (ii) by simultaneously reducing both type I and type II errors when estimating statistical support associated with the inference of the root location. Our findings provide guidelines on how to implement and use the complementary BF<sub>adj</sub> to detect and mitigate the impact of sampling bias on the outcomes of discrete phylogeographic inference using CTMC modeling.

## Simulation framework

In this study, we use a simulation framework (i) to optimise the setup to calculate the BF<sub>adj</sub> and (ii) to formally evaluate the statistical performance of the BF<sub>adj</sub>, specifically assessing how well the BF<sub>adj</sub> can identify transition events or inferred root locations for which a high BF<sub>std</sub> support likely results from  high relative abundance of the involved locations (i.e., sampling biases). 

We used 30 simulated epidemics of rabies virus (RABV) in dogs in Morocco from the study of Layan and colleagues; full study is available [here](https://doi.org/10.1093/ve/vead010). 
In this study RABV epidemics were simualted among domestic dog populations in Morocco using a stochastic metapopulation model. From each case, RABV whole-genomes associated with each case were simulated using a simple HKY model. Viral genomes were subsenquently sampled in a biased way to generate dataset of 150 or 500 sequences. Different degrees of samplig bias were implemented with weights increasing the sampling probability of genomes from particular regions by a factor of 2.5, 5, 10, 20, or 50. From each of the simulated transmission chains, we then extracted the “true” phylogeny linking the sampled genomes by pruning all but the sampled infections. 

Two types of discrete phylogeographic analyses were performed using the software package BEAST 1.10.5 on each of the simulated RABV datasets: a “standard discrete phylogeographic analysis” and a “tip-state-swap discrete phylogeographic analysis” where the location states at the tips were randomly permuted during the run.

## Tip-state-swap discrete phylogeographic analysis

XML files were generated with generate_tipSwap_xml.r which modifies the template xml file : template_tipSwap.xml

For run these analysis, 100 evenly sampled post-burnin trees from the corresponding standard discrete phylogeographic analysis were used as an empirical tree distribution. Similarly as the standard discrete phylogeographic analysis, tip-state-swap analysis were run for 20 million or 40 million Markov chain Monte Carlo (MCMC) iterations, with samples collected every 20,000 or 40,000 iterations, for the 150 and 500 sequences dataset respectively. Tip state locations were randomised during the MCMC simulation by including the tip state swap transition kernel, which swpas the location states of two tips at each permutation event, therfore controlling the total expected number of permutation events along the chain. This parameter must be manually adjusted to reflect the number of taxa and the total MCMC length, following the equations below. Let:

  - tnp = total number of expected permutation events over the MCMC chain
  - np = number of permutation events between two consecutive sampled posterior trees
  - n = number of taxa
  - x = proportion of tips swapped between two consecutive sampled posterior trees

Then:

   tip state swap transition kernel = (tnp*sum of all transition kernel weights)/MCMC length

   tnp = np*[number of posterior samples]

   np = ln(1-x)/ln((n-2)/n)

## Estimation of the adjusted Bayes Factor (BF<sub>adj</sub>)

Following [Vrancken and colleagues](https://journals.asm.org/doi/10.1128/jvi.00683-20), the BF<sub>adj</sub> support for a transition link between locations i and j was calculated by replacing the default prior expectation based on the minimal rate configuration, which only depends on the number of sampled locations, with an empirical prior expectation p<sub>emp,i→j</sub> that accounts for the relative abundance of sampled locations (i.e., sampling intensity).The p<sub>emp,i→j</sub> is the mean posterior inclusion frequency for the transition link i to j obtained from the tip-state-swap discrete phylogeographic analysis while the p<sub>i→j</sub> is the mean posterior inclusion frequency for the transition link i to j obtained from the standard discrete phylogeographic analysis. The BF<sub>adj,i→j</sub> is then calculated as follows: 


BF<sub>adj,i→j</sub> = [p<sub>i→j</sub> /(1-p<sub>i→j</sub>)]/[p<sub>emp,i→j</sub>/(1-p<sub>emp,i→j</sub>)]


We also estimated the BF<sub>std</sub> and BF<sub>adj</sub> for the root location. Let p<sub>i</sub> and p<sub>emp,i</sub> be the location probabilities obtained at the root of the tree from the standard and tip-state-swap discrete phylogeographic analysis and q<sub>k</sub> the prior probabilty assigning equal probability to all K locations, q<sub>k</sub> = 1/K. The BF<sub>std</sub> and BF<sub>adj</sub> for the root location are calculated as follows:


BF<sub>std,i</sub> = [p<sub>i</sub>/(1-p<sub>i</sub>)]/[q<sub>k</sub>/(1-q<sub>k</sub>)] 

    
BF<sub>adj,i</sub> = [p<sub>i</sub>/(1-p<sub>i</sub>)]/[p<sub>emp,i</sub>/(1-p<sub>emp,i</sub>)]


These calculations can be computed with calculateBF.r, which can be found in the R_script folder


