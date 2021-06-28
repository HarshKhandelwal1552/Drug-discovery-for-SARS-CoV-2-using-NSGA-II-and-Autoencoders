
from nsga2 import NSGA2Utils
from scripts.sample_generate import mml
from scripts.run import main
from pop import Population
import pandas as pd
import time
class Evolution:

    def __init__(self, problem, num_of_generations=4, num_of_individuals=10000, num_of_tour_particips=1, tournament_prob=0.9, crossover_param=2, mutation_param=5):
        self.utils = NSGA2Utils(problem, num_of_tour_particips, tournament_prob, crossover_param, mutation_param)
        self.population = None
        self.num_of_generations = num_of_generations
        self.on_generation_finished = []
        self.num_of_individuals = num_of_individuals

    def evolve(self, fathers):
        self.population = self.utils.create_initial_population(fathers)
        print('Initializing Training.....')
        main()
        print('Initial_Perito....')
        self.utils.fast_nondominated_sort(self.population)
        
        print('Calculate Crowding....')
        for front in self.population.fronts:
            self.utils.calculate_crowding_distance(front)
        
# =============================================================================
#         print('Initial_Perito....')
#         self.utils.fast_nondominated_sort(self.population)
#         print('Calculate Crowding....')
#         for front in self.population.fronts:
#             self.utils.calculate_crowding_distance(front)
#             
#         print('Q....')
#         children = self.utils.create_children(self.population)
# =============================================================================
#        [print(x.features) for x in children]

        for i in range(self.num_of_generations):

            print('Q....')
            children = self.utils.create_children(self.population, i)
            
            print('P+Q.....')
            self.population.extend(children)
            
            #[print(x.features) for x in self.population]
            print('NON-Dominated_Sort....')
            self.utils.fast_nondominated_sort(self.population)
            new_population = Population()
            front_num = 0
            print('Crowding_Sort....')

            while len(new_population) + len(self.population.fronts[front_num]) <= self.num_of_individuals:
                
                self.utils.calculate_crowding_distance(self.population.fronts[front_num])
                len(self.population.fronts[front_num])
                
                new_population.extend(self.population.fronts[front_num])
                
                front_num += 1
            
            #[print(x.features) for x in self.population]
            self.utils.calculate_crowding_distance(self.population.fronts[front_num])
            self.population.fronts[front_num].sort(key=lambda individual: individual.crowding_distance, reverse=True)
            new_population.extend(self.population.fronts[front_num][0:self.num_of_individuals-len(new_population)])
            print(f'Top Fathers: {len(new_population)}')
            self.population = new_population
            #[print(x.features) for x in self.population]
            print('NON-Dominated_Sort - NEW....')
            self.utils.fast_nondominated_sort(self.population)
            print('Crowding_Sort - NEW....')
            for front in self.population.fronts:
                self.utils.calculate_crowding_distance(front)
            
            #[print(x.features) for x in self.population.fronts[0]]
            smiles= [x.features for x in self.population]
            samples = pd.DataFrame(smiles, columns=['SMILES'])
            base_path= 'moses/dataset/data/gen_'
            samples.to_csv(base_path+str(i+1)+ '.csv', index=False)
            print('Improving Trained Model.....')
            mml(i+1)
