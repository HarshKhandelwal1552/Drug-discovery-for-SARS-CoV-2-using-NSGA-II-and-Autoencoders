# Drug-discovery-for-SARS-CoV-2-using-NSGA-II-and-Autoencoders

## Introduction: 
The problem statement was to use machine learning for drug discovery of CoV-SARS-II. We used the prexisting drugs and we also generated new valid drugs which could surpass their parents.

## Methodology:
We took the already existing drugs as the initial population and we used multiple scoring methods such as Synthetic acceccbility, Pains filter value, etc as mentioned in the NSGA-II for SARS-CoV-2 [1], for filtering the best fitting drugs for the cure in present population, the filtering was done using non-dominated sort and crowding distance sorting. Now we took the best drugs in present population and used them to train Adversarial Autoencoders (AAEs) than we generate a batch of new drugs and after testing their validitiy and other scoring restrictions we include them into the current population and apply the scoring functions and sort the population again. Now we get the new population which is better than the previous one, we keep on doing this until we reach our expectations. At the end the AAEs will be capable enough to generate better drugs and we also passed on the best previously obtained drugs which ensures that our model wouldn’t lost.

## Results and discussion:
I was able to implement the above methodology on small scale due to hardware and time constraints. This method can be applied to any do drug discovery for any virus just by changing the constraints of scoring functions and affinity scores which can be generated for every population at once using PyRx. I used ANNs but we can use a bunch of different generative netwoks such as Variational Autoencoder, etc.

## Conclusions:
We tried to improve the generation of new drugs in the paper NSGA-II for SARS-CoV-2 [1]. Instead of replacing, deleting, inseting new symbols we used Adversarial Autoencoders which will create valid coumpounds more often and we could have a track on weather our generated drugs are getting better or not. This approach is not just restricted to Covid cure, we can extend this approach further as discussed above.

## References: 
- ![NSGA-II for SARS-CoV-2](https://arxiv.org/abs/2005.02666)
- ![NSGA-II](https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf)
- ![All smiles Variational Autoencoder](https://arxiv.org/pdf/1905.13343v2.pdf)
- ![SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
