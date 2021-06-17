
import random
from pop import Population
from tqdm import tqdm
from scripts.sample_generate import generator


class NSGA2Utils:

    def __init__(self, problem,
                 num_of_tour_particips=1, tournament_prob=0.9, crossover_param=2, mutation_param=5, model= 'aae'):
        self.model= model
        self.problem = problem
        self.num_of_tour_particips = num_of_tour_particips
        self.tournament_prob = tournament_prob
        self.crossover_param = crossover_param
        self.mutation_param = mutation_param

    def create_initial_population(self, l):
        population = Population()
        for x in tqdm(l):
            individual = self.problem.generate_individual(x)
            self.problem.calculate_objectives(individual)
            population.append(individual)
        return population

    def fast_nondominated_sort(self, population):
        population.fronts = [[]]
        for individual in population:
            individual.domination_count = 0
            individual.dominated_solutions = []
            for other_individual in population:
                if individual.dominates(other_individual):
                    individual.dominated_solutions.append(other_individual)
                elif other_individual.dominates(individual):
                    individual.domination_count += 1
            if individual.domination_count == 0:
                individual.rank = 0
                population.fronts[0].append(individual)
        i = 0
        while len(population.fronts[i]) > 0:
            temp = []
            for individual in population.fronts[i]:
                for other_individual in individual.dominated_solutions:
                    other_individual.domination_count -= 1
                    if other_individual.domination_count == 0:
                        other_individual.rank = i+1
                        temp.append(other_individual)
            i = i+1
            population.fronts.append(temp)
        




    def calculate_crowding_distance(self, front):
        if len(front) > 0:
            solutions_num = len(front)
            for individual in front:
                individual.crowding_distance = 0

            for m in range(len(front[0].objectives)):
                front.sort(key=lambda individual: individual.objectives[m])
                front[0].crowding_distance = 10**9
                front[solutions_num-1].crowding_distance = 10**9
                m_values = [individual.objectives[m] for individual in front]
                scale = max(m_values) - min(m_values)
                if scale == 0: scale = 1
                for i in range(1, solutions_num-1):
                    front[i].crowding_distance += (front[i+1].objectives[m] - front[i-1].objectives[m])/scale

    def crowding_operator(self, individual, other_individual):
        if (individual.rank < other_individual.rank) or \
            ((individual.rank == other_individual.rank) and (individual.crowding_distance > other_individual.crowding_distance)):
            return 1
        else:
            return -1

    def create_children(self, population,  gen_i):
        children = []
        generator(gen_i)
        with open('checkpoints/'+ self.model+ '_generated_' +str(gen_i) +'.csv') as f:
            content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
        smiles = [x.strip() for x in content] 
        children= self.create_initial_population(smiles)
# =============================================================================
#         while len(children) < len(population):
#             parent1 = self.__tournament(population)
#             child1= parent1
#             children.append(child1)
# =============================================================================
            #print(child1.features)
        return children


    def __mutate(self, child):
#         num_of_features = len(child.features)
        pass


    def __tournament(self, population):
        participants = random.sample(population.population, 2)
        #[print(x.features) for x in participants]

        best = None
        for participant in participants:
         #   print(participant.features)
            if best is None:
                best = participant
            elif best is not None:
                #print('y')
                if self.crowding_operator(best,participant) ==1:
                    best= participant
        return best

    def __choose_with_prob(self, prob):
        if random.random() <= prob:
            return True
        return False